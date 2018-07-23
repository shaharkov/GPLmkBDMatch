function params = updateParamsPerMethod(params,typeLmkMatching)

switch typeLmkMatching
    case 'GP_BD'
        params.typeLmk = 'GP';
    case 'GP_Euc_BD'
        params.typeLmk = 'GP_Euc';
    case'GP_nW_BD'
        params.typeLmk = 'GP_nW';
    case 'GT_BD'
        params.typeLmk = 'GT';
    case 'GT2_BD'
        params.typeLmk = 'GT';
        params.computePutativeMatches.forceIdentity = true;
    case 'GT3'
        params.typeLmk = 'GT';
        params.computePutativeMatches.forceIdentity = true;
        params.matchSurfaceLmksBD.forceIdentity = true;
    otherwise
        error('Invalid params.typeLmkMatching')
end