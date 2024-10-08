% Extract and calculate delft3D parameters for colonisation, mortality and
% post-processing from second ets on - more description in technical overview (PDF)
%----------------Let's say, ets=1 here stores the results of last run -----------------------%

%% Calculate sedimentation/erosion for mortality
if mor == 1 % if morphology is included
    depth        = vs_get(NFS,'map-sed-series','DPS','quiet'); % bed topography with morphological changes
else % in case of only hydrology
    depth        = vs_get(NFS,'map-const','DPS0','quiet'); % bed topography (without morphology) at zeta points trimmed to fit grid
end

%% Calculate flood/dry times for seed colonization
% Extract water levels
WL               = vs_get(NFS,'map-series','S1','quiet'); % Water level data at zeta points for all time-steps
% dh
waterdepth       = cellfun(@plus,depth,WL,'UniformOutput',false);
%%>>dh 2019-12-04
% Calculate the flooding frequency from waterdepth
flood_temp       = cellfun(@(x) x>0.1, waterdepth,'UniformOutput',false);
flood            = sum(cat(3,flood_temp{:}),3);
watdep_marsh     = flood_temp;
clear waterdepth

%%  Elevation for warer stagnation
if Veg_WStag == 1 && ets == 1 && year > year_ini+1 % if veg in water stagnation is protected
    ElePro   = mean(cat(3,depth{:}),3); % mean pro ele
    MWL_temp = mean(cat(3,WL{:}),3);    
    MWL      = round(MWL_temp(Ndim-1,2),2); % MWL
end
clear depth WL MWL_temp

%% Calculate 90% cumulative bed shear stress for colonization 2019-12-04
Taumax           = vs_get(NFS,'map-series','TAUMAX','quiet');
Taumax           = cellfun(@(x,y) x.*y, Taumax ,flood_temp,'UniformOutput',false); % exclude invalid bed shear stress
Taumax_temp1     = cat(3,Taumax{:}); % Transform to 3d matrix
Taumax_90        = zeros(Ndim,Mdim);
for i = 1:Ndim % 
    for j = 1:Mdim
        Taumax_90(i,j) = prctile(abs(Taumax_temp1(i,j,:)),90); % 90% max. value % prctile function should be active
    end
end
clear Taumax flood_temp Taumax_temp1 i j NFS

%%
if year == 1 && ets==1 && Restart == 0 % Cold start from pristine
    Relative_flood0                                                = flood./max(max(flood)); % relative value to seek colonization location
    Tau0                                                           = Taumax_90; % Bed shear stress to evaluate seed colonization
    Relative_flood_marsh                                           = flood./max(max(flood));
    Tau0_marsh                                                     = Taumax_90;
elseif year == year_ini && ets == 1 && Restart == 1 % Hot start
    d3dparameters.Flood(year-1).PerYear(t_eco_year,1)              = {flood./max(max(flood))}; % relative value to seek colonization location
    d3dparameters.Tau(year-1).PerYear(t_eco_year,1)                = {Taumax_90}; % Bed shear stress to evaluate seed colonization
elseif ets==1
    d3dparameters.Flood(year-1).PerYear(t_eco_year,1)              = {flood./max(max(flood))}; % relative value to seek colonization location
    d3dparameters.Tau(year-1).PerYear(t_eco_year,1)                = {Taumax_90}; % Bed shear stress to evaluate seed colonization
else
    d3dparameters.Flood(year).PerYear(ets-1,1)                     = {flood./max(max(flood))}; % relative value to seek colonization location
    d3dparameters.Tau(year).PerYear(ets-1,1)                       = {Taumax_90}; % Bed shear stress to evaluate seed colonization
end
clear flood NFS Taumax_90