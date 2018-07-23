function K = ComputeUniformization(G,options)
switch options.boundary_conditions
    case 'disc'
        bIdx=G.FindBoundaries;
        Zv=G.V(1,:)+1i*G.V(2,:);
        Zb=Zv(bIdx);
        Z0=mean(Zv);
        bValues=-log(abs(Zb-Z0));
        switch options.method
            case 'SNCCM'
                [~,K] = G.SNCCM(bIdx,bValues,[],[],2);
            case 'LSCM'
                K = G.LSCM(bIdx,bValues);
            case 'NCCM'
                K = G.NCCM(bIdx,bValues,2);
        end
        F=(Zv-Z0).*exp(K.V(1,:)+1i*K.V(2,:));
        K.V=[real(F);imag(F);zeros(size(F))];
    case 'halfplane'
        [~,bIdx]=G.FindOrientedBoundaries; bIdx=bIdx{1};
        bIdx_y=bIdx;
        bValues_y=zeros(1,numel(bIdx));
        infIdx=7;
        bIdx_x=bIdx([infIdx-1,infIdx,infIdx+1]);
        bValues_x=[-1,0,1];
        F=G.F;
        Finf_ind=any(F==bIdx(1),1); %infinity faces indices
        F(:,Finf_ind)=[];
        Ginf=Mesh('VF',G.V,F);
        switch options.method
            case 'SNCCM'
                [~,K]=G.SNCCM2d(bIdx_x,bValues_x,bIdx_y,bValues_y,2);
            case 'LSCM'
                K = Ginf.LSCM2d(bIdx_x,bValues_x,bIdx_y,bValues_y,2);
        end
    case 'triangle'
        
end