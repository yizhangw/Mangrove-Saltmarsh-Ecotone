%% To store model results
%>> When timesteps is large or model results is rough,
% use the original Delft3D results
%>> When timestep is small or model results are required to be elaborated,
% use the results method storing specified parameters


% 28/12/2018 Store the results in 2 different ways
if Storage == 0 % Copy the original Delft3D results to results folder
    try
        % Save results to result/folder for analysis
        copyfile(strcat(directory, 'work\trim-', ID1, '.def'), strcat(directory, 'results_', num2str(year), '\trim-', ID1, '_', num2str(ets),'.def'));
        copyfile(strcat(directory, 'work\trim-', ID1, '.dat'), strcat(directory, 'results_', num2str(year), '\trim-', ID1, '_', num2str(ets),'.dat'));
    catch
    end
else % Storage =1, convert the original results to my own results with specified parameters
    if mod(year,1)==0 && ets == 12 && marsh==1 % Save the last ets to results folder in order to the restart process
        try
            % Save results to result/folder for analysis
            copyfile(strcat(directory, 'work\trim-', ID1, '.def'), strcat(directory, 'results_', num2str(year), '\trim-', ID1, '_', num2str(ets),'.def'));
            copyfile(strcat(directory, 'work\trim-', ID1, '.dat'), strcat(directory, 'results_', num2str(year), '\trim-', ID1, '_', num2str(ets),'.dat'));
            save(strcat(directory, 'results_', num2str(year), '\d3dparameters.mat'),'d3dparameters');
        catch
        end
    elseif mod(year,1)==0 && ets == 12
        try
            % Save results to result/folder for analysis
            copyfile(strcat(directory, 'work\trim-', ID1, '.def'), strcat(directory, 'results_', num2str(year), '\trim-', ID1, '_', num2str(ets),'.def'));
            copyfile(strcat(directory, 'work\trim-', ID1, '.dat'), strcat(directory, 'results_', num2str(year), '\trim-', ID1, '_', num2str(ets),'.dat'));
            save(strcat(directory, 'results_', num2str(year), '\d3dparameters.mat'),'d3dparameters');
        catch
        end
    end
    % Function to save the specified parameters
    %>>Read NFS
    NFS      = vs_use(strcat(directory, 'work\', 'trim-', ID1,'.dat'),'quiet'); % read last trim file from previous year fine domain
    % water depth/DPS
    DPS      = vs_get(NFS,'map-sed-series','DPS','quiet'); % Bathymetry
    % water level
    S1       = vs_get(NFS,'map-series','S1','quiet'); % Water level
    % depth-averaged velocity
    U1       = vs_get(NFS,'map-series','U1','quiet'); % X-Flow velocity
    V1       = vs_get(NFS,'map-series','V1','quiet'); % X-Flow velocity
    % sediment concentration
    R1       = vs_get(NFS,'map-series','R1','quiet');  % Concentration per layer
    % bed shear stress
    Tauksi   = vs_get(NFS,'map-series','TAUKSI','quiet'); % Bed shear stress in U-point
    Taueta   = vs_get(NFS,'map-series','TAUETA','quiet'); % Bed shear stress in V-point
    %>>Save data
    MyStorage.DPS    = DPS;
    MyStorage.S1     = S1;
    MyStorage.U1     = U1;
    MyStorage.V1     = V1;
    MyStorage.R1     = R1;
    MyStorage.Tauksi = Tauksi;
    MyStorage.Taueta = Taueta;
    if Wave > 0
        % maximum bed shear stress
        Taumax   = vs_get(NFS,'map-series','TAUMAX','quiet');
        MyStorage.Taumax = Taumax; % max. bed shear stress
        % significant wave height
        Hs       = vs_get(NFS,'map-rol-series','HS','quiet');
        MyStorage.Hs = Hs;
    end
    
    save('MyStorage.mat','MyStorage');
    copyfile(strcat(directory, 'work\MyStorage.mat'), strcat(directory, 'results_', num2str(year), '\MyStorage_', num2str(ets),'.mat'));
    clear NFS DPS S1 U1 V1 R1 Tauksi Taumax MyStorage
end