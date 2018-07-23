function [ca] = CORR_evaluate_disc_moebius_from_tet(ctet,z,w)

ca = zeros(size(ctet));
for k=1:length(ctet)
    tet=ctet(k);

    %solve on linear equation for a and tet
    A = z*w; A1 = real(A); A2 = imag(A);
    B = -exp(1i*tet); B1 = real(B); B2 = imag(B);
    C = -exp(1i*tet)*z+w; C1=real(C); C2=imag(C);
    Q = [(A1+B1) (A2-B2); (A2+B2) (B1-A1)];
    T = [C1 C2]';
    tca = Q\T; %the candidate a matching the given tet
    ca(k)=tca(1)+1i*tca(2);
end