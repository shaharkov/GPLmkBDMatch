function [tp] = DISCtoPLANE(p, type)
%DISCTOPLANE: map disk to the (complex) plane
%   p: nx2 matrix

p = p(:,1)+1i*p(:,2);
rads = abs(p);
if strcmpi(type, 'd2p') % forward map
    tp = atanh(rads).*p./rads;
elseif strcmpi(type, 'p2d') % backward map
    tp = tanh(rads).*p./rads;
end
%the zeros are fixed
tp(p.*conj(p)<10^-20) = 0;

tp = [real(tp), imag(tp)];

end


