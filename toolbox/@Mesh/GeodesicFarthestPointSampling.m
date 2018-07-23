function [VertSampInd] = GeodesicFarthestPointSampling(G,SampNumb,InitialSamples)
%FARTHESTPOINTSAMPLING Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    SampNumb = G.nV;
end
if nargin<3
    rng('shuffle');
    InitialSamples = randi(G.nV,1,1);
else
    if size(InitialSamples,1)>size(InitialSamples,2)
        InitialSamples = InitialSamples';
    end
end

ProcessedSampNumb = length(InitialSamples);
VertSampInd = [InitialSamples, zeros(1,SampNumb-ProcessedSampNumb)];
for k=(ProcessedSampNumb+1):SampNumb
    progressbar(k,SampNumb,20);
    [~,VertSampInd(k)] = max(G.PerformFastMarching(VertSampInd(1:(k-1))));
end

end

