%--------------------------------------------------------------------------
% Bounded Distortion Feature Matching
%
% PlotMatches.m 
%
% This code, written by Eli Libman, advised by Meirav Galun, 
% is an implementation of the algorithm described in the paper 
% Feature Matching with Bounded Distortion
% by Yaron Lipman, Stav Yagev, Roi Poranne, David W. Jacobs, Ronen Basri
% ACM Trans. Graph., Volume 33 Issue 3, May 2014 
%
% Disclaimer: The code is provided as-is and without any guarantees. 
%--------------------------------------------------------------------------


function fig = PlotMatches( Ia, Ib, fa,fb,matches )
[m n c] =size(Ia); % assume images have same size
P=fa(1:2,matches(1,:));
Q=fb(1:2,matches(2,:));

[umatch non2uniq uniq2non] =unique(P','rows','first');
P=fa(1:2,matches(1,non2uniq));
Q=fb(1:2,matches(2,non2uniq));

matches=matches(:,non2uniq);
ia=matches(1,:);
ib=matches(2,:);
fig=figure;
imshow([Ia Ib]);
hold on
scatter(fa(1,ia),fa(2,ia),100,fa(1,ia),'fill');
scatter(n+fb(1,ib),fb(2,ib),100,fa(1,ia),'fill');
scatter(fa(1,ia),fa(2,ia),100,fa(1,ia),'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
scatter(n+fb(1,ib),fb(2,ib),100,'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
%
%scatter(fa(1,ia),fa(2,ia),3+(fa(2,ia)/m)*100,fa(1,ia),'fill');
%scatter(n+fb(1,ib),fb(2,ib),3+(fa(2,ia)/m)*100,fa(1,ia),'fill');
%
%scatter(fa(1,ia),fa(2,ia),3+(fa(2,ia)/m)*50,fa(1,ia),'marker','o');
%scatter(n+fb(1,ib),fb(2,ib),3+(fa(2,ia)/m)*50,fa(1,ia),'marker','o');
hold off;
end

