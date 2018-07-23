//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% Code implementing the paper "Large-Scale Bounded Distortion Mappings".
//% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
//%             Please contact the author to report any bugs.
//% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//#define EIGEN_USE_MKL_ALL
#include "mex.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "mexHelpers.cpp"
//#include <chrono>

using namespace Eigen;

#define NOW std::chrono::system_clock::now

void projBlockBD2x2(VectorXd &pA, int dim, double K, double lb, double ub, VectorXd& distortions, VectorXd& flips, VectorXd& minsv, VectorXd& maxsv)
{
	int block_size = dim*dim;
	int num_blocks = pA.size() / block_size;
	double K2plus1 = K*K + 1;
	JacobiSVD<Matrix2d> svdA(dim, dim, (ComputeFullU | ComputeFullV));
	Vector2d s, snew;
	double dist, t;
	bool flipped;
	Matrix2d U, V;
	Map<Matrix2d> currA(pA.data());
	double a, b, c, d;
	double n_ab, n_cd;
	bool wasProjected;


	// project
	for (int ii = 0; ii < num_blocks; ii++)
	{
		// get current block
		new (&currA) Map<Matrix2d>(pA.data() + ii*block_size);
		
		// 2x2 sv's and sign of det
		a = (currA(0, 0) + currA(1, 1)) / 2;
		b = (currA(0, 1) - currA(1, 0)) / 2;
		c = (currA(0, 0) - currA(1, 1)) / 2;
		d = (currA(0, 1) + currA(1, 0)) / 2;
		n_ab = sqrt(a*a + b*b);
		n_cd = sqrt(c*c + d*d);
		flipped = (n_cd > n_ab);
		s(0) = n_ab + n_cd;
		s(1) = (flipped ? (n_cd - n_ab) : (n_ab - n_cd));
		// calc distortion
		dist = s(0) / s(1);
		// update sign of last singular value
		if (flipped)
			s(1) = -s(1);

		// -- projections --
		wasProjected = false;
		snew = s;
		// project BD
		if ((K > 0) & (flipped | (dist > K)))
		{
			t = (K*s(0) + s(1)) / K2plus1;
			snew << t*K, t;
			wasProjected = true;
		}
		// project LB
		if ((lb > 0) & (s(1) < lb))
		{
			snew = snew.array().max(lb);
			wasProjected = true;
		}
		// project LB
		if ((ub > 0) & (s(0) > ub))
		{
			snew = snew.array().min(ub);
			wasProjected = true;
		}
		// actually project the matrix
		if (wasProjected)
		{
			// full svd
			svdA.compute(currA);
			U = svdA.matrixU();
			V = svdA.matrixV();
			// ssvd
			if (flipped)
				U.col(1) = -U.col(1);
			// update block
			currA = U*snew.asDiagonal()*V.transpose();
		}
		// report back things
		distortions(ii) = dist;
		flips(ii) = flipped;
		maxsv(ii) = s(0);
		minsv(ii) = s(1);
	}
}

void projBlockBD3x3(VectorXd &pA, int dim, double K, double lb, double ub, VectorXd& distortions, VectorXd& flips, VectorXd& minsv, VectorXd& maxsv)
{
	int block_size = dim*dim;
	int num_blocks = pA.size() / block_size;
	double K2plus1 = K*K + 1;
	double K2plus2 = K*K + 2;
	double TwoK2plus1 = 2 * K*K + 1;
	JacobiSVD<Matrix3d> svdA(dim, dim, (ComputeFullU | ComputeFullV));
	Vector3d s, snew;
	double dist, t;
	bool flipped;
	Matrix3d U, V;
	Map<Matrix3d> currA(pA.data());
	bool wasProjected;

	int c = 0;
	// project
	for (int ii = 0; ii < num_blocks; ii++)
	{
		// get current block
		new (&currA) Map<Matrix3d>(pA.data() + ii*block_size);
		// sign of determinant 
		flipped = (currA.determinant() < 0);
		flips(ii) = flipped;
		// svd
		svdA.compute(currA);
		s = svdA.singularValues();
		// calc distortion
		dist = s(0) / s(2);
		// update sign of last singular value
		if (flipped)
			s(2) = -s(2);

		// -- projections --
		wasProjected = false;
		snew = s;
		// project BD
		if ((K > 0) & (flipped | (dist > K)))
		{
			// project
			// first, ignore central singular value (assume it remains unchanged)
			t = (K*s(0) + s(2)) / K2plus1;
			snew << t*K, s(1), t;
			// now, check if need to fix central singular value
			if (snew(1)<snew(2))
			{
				t = (K*s(0) + s(1) + s(2)) / K2plus2;
				snew << t*K, t, t;
			}
			else if (snew(1)>snew(0))
			{
				t = (K*(s(0) + s(1)) + s(2)) / TwoK2plus1;
				snew << t*K, t*K, t;
			}
			wasProjected = true;
		}
		// project LB
		if ((lb > 0) & (snew(2) < lb))
		{
			snew = snew.array().max(lb);
			wasProjected = true;
		}
		// project LB
		if ((ub > 0) & (snew(0) > ub))
		{
			snew = snew.array().min(ub);
			wasProjected = true;
		}
		// actually project the matrix
		if (wasProjected)
		{
			// full svd
			svdA.compute(currA);
			U = svdA.matrixU();
			V = svdA.matrixV();
			// ssvd
			if (flipped)
				U.col(2) = -U.col(2);
			// update block
			currA = U*snew.asDiagonal()*V.transpose();
		}
		// report back things
		distortions(ii) = dist;
		flips(ii) = flipped;
		maxsv(ii) = s(0);
		minsv(ii) = s(2);
	}
}


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])
{
	// init timer
	//std::chrono::time_point<std::chrono::system_clock> start, end;
	//std::chrono::duration<double> elapsed_seconds;

	if (nrhs!=5)
		mexErrMsgIdAndTxt("MATLAB:wrong_number_of_inputs", "must have 5 inputs");

	// assign input
	int A_rows = mxGetM(prhs[0]); // # rows of A
	int A_cols = mxGetN(prhs[0]); // # cols of A
	double *dim, *K, *lb, *ub;
	const Map<VectorXd> A(mxGetPr(prhs[0]), A_rows, A_cols);
	dim = mxGetPr(prhs[1]);
	K = mxGetPr(prhs[2]); 
	lb = mxGetPr(prhs[3]);
	ub = mxGetPr(prhs[4]);
	
	// init output
	VectorXd pA(A_rows);
	VectorXd distortions(A_rows / ((*dim)*(*dim)));
	VectorXd flips(A_rows / ((*dim)*(*dim)));
	VectorXd minsv(A_rows / ((*dim)*(*dim)));
	VectorXd maxsv(A_rows / ((*dim)*(*dim)));
	pA = A;

	// project
	//start = NOW();
	if (*dim == 2)
		projBlockBD2x2(pA, *dim, *K, *lb, *ub, distortions, flips, minsv, maxsv);
	else if (*dim == 3)
		projBlockBD3x3(pA, *dim, *K, *lb, *ub, distortions, flips, minsv, maxsv);
	else
		mexErrMsgIdAndTxt("MATLAB:wrong_dimension", "dim must be either 2 or 3");
	//end = NOW();
	//elapsed_seconds = end - start;
	//mexPrintf("projection took %f secs\n", elapsed_seconds);
	
	// assign outputs
	mapDenseMatrixToMex(pA, &(plhs[0]));
	mapDenseMatrixToMex(distortions, &(plhs[1]));
	mapDenseMatrixToMex(flips, &(plhs[2]));
	mapDenseMatrixToMex(minsv, &(plhs[3]));
	mapDenseMatrixToMex(maxsv, &(plhs[4]));
}