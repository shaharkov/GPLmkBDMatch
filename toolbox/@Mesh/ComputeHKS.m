 function HKS=ComputeHKS(G,Nt,Neig)
% compute the Heat Kernel Signature according to Sun,Osjanikov,Guibas
% t: Time scale as [t1,t2]
% Nt: number of samples
[L,M] = ComputeLaplacian(G);
[v,d] = eigs(M*L,Neig,'SM',struct('v0',[0,1,zeros(1,size(L,2)-2)]'));
d = abs(diag(d));
% if strcmp(t,'auto')
%     t=0;
%     t(1)=4*log(10)/d(end);
%     %t(2)=4*log(10)/d(floor(Neig/2));
%     t(2)=4*log(10)/d(2);
% end
% t=logspace(log10(t(1)),log10(t(2)),Nt);
% 
% expdt=abs(exp(-d*t)); % this is exp([d1*t1,d1*t2,d1*t3,...;d2*t1,d2*t2,d2*t3,...])
% expdtv=abs(v).^2*expdt;   % matrix where each column is the hks for time t
% HKS=expdtv./repmat(mean(expdtv,1),size(expdtv,1),1); % normalze HKS for aech t

%% test
evals = d;
evecs = v;
tmin = abs(4*log(10) / evals(end));
tmax = abs(4*log(10) / evals(2));
nstep = Nt-1;

stepsize = (log(tmax) - log(tmin)) / nstep;
logts = log(tmin):stepsize:log(tmax);
ts = exp(logts);

hks = abs( evecs(:, 2:end) ).^2 * exp( ( abs(evals(2)) - abs(evals(2:end)) )  * ts);
Am = M;
colsum = sum(Am*hks);
scale = 1.0./ colsum; 
scalem = sparse([1:length(scale)], [1:length(scale)], scale);
HKS = hks * scalem;
