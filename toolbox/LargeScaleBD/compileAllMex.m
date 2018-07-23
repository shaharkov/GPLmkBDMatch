%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Large-Scale Bounded Distortion Mappings".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
cd mex


%% compile all
eigenFolder = 'eigen-eigen-b30b87236a1b';

mex('-largeArrayDims',['-I',eigenFolder],'computeMeshTranformationCoeffsMex.cpp');
mex('-largeArrayDims',['-I',eigenFolder],'projectBDMex.cpp');


%%
cd ..