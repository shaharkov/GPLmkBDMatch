function [not_common] = CORR_find_not_common(a,b)%find common element of number vectors a and b


not_common=[];
a=sort(a);
b=sort(b);
while (~isempty(a) && ~isempty(b) )
    if(a(1)==b(1))
        a(1)=[];
        b(1)=[];
    elseif (a(1) > b(1))
        not_common = [not_common b(1)];
        b(1)=[];        
    else
        not_common = [not_common a(1)];
        a(1)=[];
    end   
end



not_common = [not_common a b];