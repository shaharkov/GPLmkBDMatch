
Bounded Distortion Feature Matching


This MATLAB code, written by Eli Libman, advised by Meirav Galun, 
is an implementation of the algorithm described in the paper 
"Feature Matching with Bounded Distortion"
by Yaron Lipman, Stav Yagev, Roi Poranne, David W. Jacobs, Ronen Basri
ACM Trans. Graph., Volume 33 Issue 3, May 2014 

Start by running demo.m to see how to filter candidate SIFT matches of two images.

Dependencies:
VLFeat http://www.vlfeat.org/ - used for candidate SIFT calculaitons.
YALMIP http://users.isy.liu.se/johanl/yalmip/ - the optimization framework
MOSEK https://mosek.com/ - used for the SOCP solver

The example images are taken from: 
"Photo tourism: Exploring photo collections in 3d"
Snavely N., Seitz S., Szeliski R. 2006. 
ACM Trans. Graph., Volume 25 Issue 3, July 2006 

Disclaimer: The code is provided as-is and without any guarantees. 
