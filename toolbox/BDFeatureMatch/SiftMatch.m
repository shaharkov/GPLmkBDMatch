%--------------------------------------------------------------------------
% Bounded Distortion Feature Matching
%
% SiftMatch.m 
%
% This code, written by Eli Libman, advised by Meirav Galun, 
% is an implementation of the algorithm described in the paper 
% Feature Matching with Bounded Distortion
% by Yaron Lipman, Stav Yagev, Roi Poranne, David W. Jacobs, Ronen Basri
% ACM Trans. Graph., Volume 33 Issue 3, May 2014 
%
% Disclaimer: The code is provided as-is and without any guarantees. 
%--------------------------------------------------------------------------

function [m12, s12, F1, F2] = SiftMatch(img1,img2,dense,dense_target, thresh, sift_type)
m = 1000; % number of correspondences to plot
sk = 1;
im1 = GraySingle(img1);
im2 = GraySingle(img2);
% Sift them
if dense > 0
   
    switch sift_type
        case 'dsift'
            [F1,D1] = vl_dsift(im1,'step',dense);
            [F2,D2] = vl_dsift(im2,'step',dense);
           
        case 'sift'
           
            scales = [4 6 8 10 12]/3;
           
            % build SIFT for image 1
            if(dense==0)
                [F1,D1] = vl_sift(im1);
            else
                siz1=size(img1);
                [xx1 yy1] = meshgrid(1:dense:siz1(2),1:dense:siz1(1));
                tn1 = length(xx1(:));
                fc1 = [];
                for sc = scales
                    temp = [xx1(:)' ; yy1(:)' ; sc*ones(1,tn1) ; zeros(1,tn1)];
                    fc1 = [fc1 temp];
                end
            end
            % build SIFT for image 2
            siz2=size(img2);
            dense = dense_target;
            %     dense=3;
            [xx2 yy2] = meshgrid(1:dense:siz2(2),1:dense:siz2(1));
            tn2 = length(xx2(:));
            fc2 = [];
            for sc = mean(scales)%(4)
                temp = [xx2(:)' ; yy2(:)' ; sc*ones(1,tn2) ; zeros(1,tn2)];
                fc2 = [fc2 temp];
            end
           
            [F1,D1] = vl_sift(im1,'frames',fc1,'orientations') ;
            [F2,D2] = vl_sift(im2,'frames',fc2,'orientations') ;
                    
    end
   
else
   
    
    [F1,D1] = vl_sift(im1);
    [F2,D2] = vl_sift(im2);
end
% and match
[m12,s12] = vl_ubcmatch(D1,D2, thresh);

end
function I = GraySingle(Ia)
I= im2single(rgb2gray(Ia));
end