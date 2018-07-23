function [newarr ind] = CORR_sort_and_delete_repeats(arr)
%sort arr and delete repeating elements

[newarr I] = sort(arr); %delete repeating vertices
d=diff(newarr); 
ind = find(d == 0)+1;
newarr(ind)=[];
ind=I(ind);

% 
% arr = sort(arr); %delete repeating vertices
% d=diff(arr); 
% ind = find(d == 0)+1;
% arr(ind)=[];
% 
% newarr =arr;