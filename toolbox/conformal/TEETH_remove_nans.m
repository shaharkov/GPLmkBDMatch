function [A] = TEETH_remove_nans(A)

temp = A(1);
while (isnan(temp))
    A(1)=[];
    temp = A(1);
    
end