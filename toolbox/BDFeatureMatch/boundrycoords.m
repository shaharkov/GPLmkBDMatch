%--------------------------------------------------------------------------
% Bounded Distortion Feature Matching
%
% boundarycoords.m
%
% This code, written by Eli Libman, advised by Meirav Galun, 
% is an implementation of the algorithm described in the paper 
% Feature Matching with Bounded Distortion
% by Yaron Lipman, Stav Yagev, Roi Poranne, David W. Jacobs, Ronen Basri
% ACM Trans. Graph., Volume 33 Issue 3, May 2014 
%
% Disclaimer: The code is provided as-is and without any guarantees. 
%--------------------------------------------------------------------------



function boundry = boundrycoords(bb,PAD,N)
orgx=bb(2)-bb(1);
orgy=bb(4)-bb(3);
padx=(orgx)*PAD;
pady=(orgy)*PAD;
ngrid=sqrt(N);
dx=(orgx + 2*padx)/ngrid;
dy=(orgy + 2*pady)/ngrid;
[gx gy]=meshgrid((bb(1)-padx):dx:(bb(2)+padx),(bb(3)- pady):dy:(bb(4)+pady));
gridp=[gx(:) gy(:)];
% in_boundry=gx(:)<bb(1)-padx+dx/2 | ...
%             gx(:)>bb(2)+padx-dx/2 | ...
%             gy(:)<bb(3)-pady+dy/2 | ...
%             gy(:)>bb(4)+pady-dy/2;
        
        in_boundry=gx(:)==min(gx(:)) | ...
            gx(:)==max(gx(:)) | ...
            gy(:)==min(gy(:))| ...
            gy(:)==max(gy(:));

boundry=gridp(in_boundry,:);

end