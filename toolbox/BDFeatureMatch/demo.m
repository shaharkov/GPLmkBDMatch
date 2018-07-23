%--------------------------------------------------------------------------
% Bounded Distortion Feature Matching
%
% demo.m
% run this for example on how to use this code for filtering SIFT candidate
% matches.
%
% This code, written by Eli Libman, advised by Meirav Galun, 
% is an implementation of the algorithm described in the paper 
% Feature Matching with Bounded Distortion
% by Yaron Lipman, Stav Yagev, Roi Poranne, David W. Jacobs, Ronen Basri
% ACM Trans. Graph., Volume 33 Issue 3, May 2014 
%
% Disclaimer: The code is provided as-is and without any guarantees. 
%--------------------------------------------------------------------------




%initialize
%run('F:\Tools\Matlab_tools\vlfeat-0.9.16\toolbox\vl_setup');
display('Bounded Distortion Feature matching demo. libraries vlfeat,yalmip and mosek are expected to be in your path.')
clear;close all;yalmip('clear');
IMAGES_DIR='./NotreDame';
K=3;
pnorm=0.001;

%--------------------------------------------------------------------------
% SIFT parameters
%--------------------------------------------------------------------------
dense = 0;%used: 5; % 0 for sparse sift, >0 for dense (# is step size)
dense_target = 0; %used: 5
thresh = 1.5;% used: 1.35; %the threshold for choosing a sift (default = 1.5)
sift_type = 'sift';% 'sift','dsift'
%--------------------------------------------------------------------------
% choose an image pair
%--------------------------------------------------------------------------
% NotreDame
Ia = imread(fullfile(IMAGES_DIR,'NotreDame-a.jpg'));
Ib = imread(fullfile(IMAGES_DIR,'NotreDame-b.jpg'));
% --------------------------------------------------------------------
%                                           Extract features and match
% --------------------------------------------------------------------
[matches, scores, fa, fb] = SiftMatch(Ia,Ib,dense,dense_target, thresh, sift_type);
PlotMatches(Ia,Ib,fa,fb,matches);

P=fa(1:2,matches(1,:));
Q=fb(1:2,matches(2,:));


[umatch, non2uniq, uniq2non] =unique(P','rows','first');
P=fa(1:2,matches(1,non2uniq));
Q=fb(1:2,matches(2,non2uniq));


delta=sqrt( ((max(P(1,:))-min(P(1,:)))^2  + (max(P(2,:))-min(P(2,:)) )^2)/4 );
tic;
% [ Pf, Qf, W, deformation] = BDMatch(P, Q, K, pnorm, true,true,false);
[ Pf, Qf, W, deformation] = BDMatch(P, Q, K, pnorm, false,true,false);

toc

%figure
PlotMatches( Ia, Ib, Pf,Qf,[1:size(Pf,2);1:size(Pf,2)])







