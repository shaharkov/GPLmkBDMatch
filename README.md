Bounded distortion matching of Gaussian Process Landmarks
====

A Matlab implementation accompanying the papers ["Gaussian Process Landmarking for Three-Dimensional Geometric Morphometrics"](https://arxiv.org/abs/1807.11887) and ["Gaussian Process Landmarking on Manifolds"](https://arxiv.org/abs/1802.03479).

----
The functions `code/wrapperLmkMatch.m` and `code/matchSurfaceLmksBD.m` implement the matching algorithm described in section 4.2 of ["Gaussian Process Landmarking for Three-Dimensional Geometric Morphometrics"](https://arxiv.org/abs/1807.11887).

This package includes two demo pipelines:
- `demo_GPLmkMatching.m` -- implements a full pipeline for running the GP-BD matching algorithm described in the paper.
- `demo_LmkMatching_versus_CPM.m` -- a complete pipeline that includes comparison to Continuous Procrustes inter-surface Maps (CPM) as well as the alternative registration methods discussed in the paper (see Table 1).

The demo output includes the following figures, visualizing the various steps of the algorithm:
![Example Figures](./.images/figures.png?raw=true)






**Compatibility and dependencies:**
- Developed and tested with Matlab 2018a on a Windows (x64) Machine.
- Included with the code are (modified versions) of
	- [Feature Matching with Bounded Distortion](http://www.wisdom.weizmann.ac.il/~ylipman/bd_feature_match/BDFeatureMatch.zip)
	- [Accelerated Quadratic Proxy for Geometric Optimization](https://services.math.duke.edu/~shaharko//AcceleratedQuadraticProxy.html)
	- [Large-Scale Bounded Distortion Mappings](https://services.math.duke.edu/~shaharko//LargeScaleBD.html)
	- [YALMIP](https://yalmip.github.io/)
- Requires a working installation of [MOSEK](https://www.mosek.com/).
- Depends on several MEX functions. Source code is provided; compilation requires [Eigen](http://eigen.tuxfamily.org/); fast implementation of As-Rigid-As-Possible also relies on [libigl](https://github.com/libigl/libigl).

**Disclaimer:**
The code is provided as-is for academic use only and without any guarantees. Please contact the authors to report any bugs. 
Written by [Shahar Kovalsky](https://services.math.duke.edu/~shaharko/) and [Tingran Gao](https://gaotingran.com/).
