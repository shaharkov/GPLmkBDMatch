function Gs = loadAndPreprocess(folder,filename,params)


%%
dataFileName = fullfile(folder,[filename '.' params.dataFileType]);
preprocessFileName = fullfile(folder,[filename params.suffix '.mat']);
GTLmkFileName = fullfile(folder,[filename '_GTLmk.mat']);


%% try to load existing preprocessed data
if ~params.force 
    try
        data = load(preprocessFileName,'params');
        if isequal(data.params,params)
            data = load(preprocessFileName);
            Gs = data.Gs;
            fprintf('loaded preprocessed data for %s\n', dataFileName);
            return
        end
    end
end
        

%% common
fprintf('loading and preprocessing %s\n', dataFileName);
% load
Gs = Mesh('off',dataFileName);
% normalize
Gs.Normalize();
% reorient (outward facing normals)
[~,~,flip] = Gs.ComputeNormal();
if flip
    Gs.F = Gs.F([1 3 2],:);
end


%% GPLmk preprocess
% compute WKS
Gs.Aux.WKS = Gs.ComputeWKS([]);


%% CPM preprocess
% preprocess for CPM
Gs.ComputeMidEdgeUniformization(params.CPM);


%% try to load GT Lmk data
try
    temp = load(GTLmkFileName);
    Gs.Aux.GTLmks = temp.GTLmks;
catch
    Gs.Aux.GTLmks = nan;
end


%% save preprocessed data
save(preprocessFileName,'Gs','params');