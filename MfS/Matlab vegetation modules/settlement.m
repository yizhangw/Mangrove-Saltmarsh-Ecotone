%% Settlement module. Vegetation is assigned to grid cells with certain area fraction.
% Find the cells with space for colonization and the fraction that can be colonized. Fills up the space with
% the new fraction from current ETS==1. Adds the trachytope ID to new fractions for each vegetation type and saves it in matrix output.
% Opens existing trv-file and adds the new fractions, sorted after the cell numbers.
%----------------Man, say sth--------------------**--------
%>>dh: I am a referee, so allocate area fraction when coexist
%>>dh: I am a donkey, area fraction always represents no. of plants by timing total no.

%% construct root information based on the trv_trd with only stems
if Root == 1 % Include mangrove roots in the model 23/10/2018
    % duplicate root matrix from trv_trd
    trv_trd_root = trv_trd;
    % Find mangrove species and calculate the roots number
    for nv = 1:num_veg_types
        Loc_sp    = find(trv_trd(:,11)==nv); % Red = 1; Black = 2; White = 3
        % 12/09/2018 new roots number increase method-exponential increase
        trv_trd_root(Loc_sp,12) = max_root(nv)*(1./(1+exp(f*(Dmax(nv)/2-trv_trd_root(Loc_sp,13))))).*trv_trd(Loc_sp,12);
        clear Loc_sp
    end
    % 29/08/2018 Sum the different roots at different rows in particular cell into one row
    root_rc   = unique(trv_trd_root(:,[1,2]),'rows'); % cells with root
    root_temp = zeros(size(root_rc,1),size(trv_trd_root,2)); % build a temperary matrix
    for i = 1:size(root_rc,1)
        Loc             = find(trv_trd_root(:,1)==root_rc(i,1) & trv_trd_root(:,2)==root_rc(i,2)); % Target the cell number
        root_temp(i,:)  = trv_trd_root(Loc(1),:); % Attribute the first row to root matrix
        root_temp(i,11) = 900; % represent the root number has already been accumulated
        root_temp(i,12) = sum(trv_trd_root(Loc,12)); % Sum the number of roots
    end
    clear trv_trd_root root_rc i Loc
    trv_trd_root = root_temp;
    % Refresh the area fraction 4th, height 7th, density 8th, diameter 13th
    %         trv_trd_root(:,4)  = trv_trd;
    trv_trd_root(:,7)  = 0.15; % root height is a constant value: 15 cm
    trv_trd_root(:,8)  = trv_trd_root(:,12)./S_cell.*trv_trd_root(:,13)/100; % root diameter is a constant value, 1 cm
    trv_trd_root(:,13) = 1; % root diameter is a constant value, 1 cm
    % 11/10/2018  Drag coefficient of roots
    for nv = 1:num_veg_types
        trv_trd_root(trv_trd_root(:,11)==900,9) = Cd_root(nv);
    end
    % Combine stems with roots
    trv_trd                      = [trv_trd; trv_trd_root];
else % Exclude roots 23/10/2018
    num_all = num0; % Use a smaller maximum objexts number
end
trv_trd(trv_trd(:,1)==0,:)  = [];
trv_trd(trv_trd(:,2)==0,:)  = [];
%% Mangroves die early under the shadow of saltmarshes
if marsh==1 && ets==12
    if ets==1 && year==1 && Restart==0
        PH_mort = zeros(Ndim,Mdim);
    else
    PH_mort = PH_new;
    end

    IC_rc_temp = unique(trv_trd(:,[1,2]),'rows');
    for i = 1:size(IC_rc_temp,1)
        kkk = find(trv_trd(:,1)==IC_rc_temp(i,1) & trv_trd(:,2)==IC_rc_temp(i,2));
        Marsh_height = PH_mort(IC_rc_temp(i,1),IC_rc_temp(i,2));
        if max(trv_trd(kkk,7)) < Marsh_height && max(trv_trd(kkk,19)) >= 1 % Kills mangroves that have been growth-stunted for 2 years and are shorter than salt marshes
            trv_trd(kkk,12) = 0;
        end
    end
    trv_trd(trv_trd(:,12)==0,:) = [];
end
clear IC_rc_temp PH_mort i kkk Marsh_height

%% Predation probability linked to flooding
if Predation==1
    if ets==9 || ets==10
        if year==1 && ets==1 && Restart==0
            Relative_flood     = Relative_flood_marsh ;
            clear Relative_flood_marsh
        elseif ets==1
            Relative_flood     = cell2mat(d3dparameters.Flood(year-1).PerYear(t_eco_year,1));
        else
            Relative_flood     = cell2mat(d3dparameters.Flood(year).PerYear(ets-1,1));
        end

        Ppred_temp = Relative_flood;
        Ppred_temp(Ppred_temp(:)<0.8) = 0;
        Ppred_temp(Ppred_temp(:)>0) = 0.11;

        loc_temp = unique(trv_trd(:,[1,2]),'rows');

        for i = 1:size(loc_temp,1)
            kkk = find(trv_trd(:,1)==loc_temp(i,1) & trv_trd(:,2)==loc_temp(i,2));
            mloc = trv_trd(kkk(1),1);
            nloc = trv_trd(kkk(1),2);
            Ppred_P_temp = rand(1)<Ppred_temp(mloc,nloc);

            if Ppred_P_temp == 1 && max(trv_trd(kkk,7)) < 1 % Kill short seedlings
                trv_trd(kkk,12) = trv_trd(kkk,12).*0;
            end
        end
        clear i Ppred_temp mloc nloc Ppred_P_temp Flood_temp loc_temp kkk
    end
end
trv_trd(:,12) = round(trv_trd(:,12),0);
trv_trd(trv_trd(:,12)==0,:)  = []; % delete the invalid rows
trv_trd(trv_trd(:,19)>=5,:)  = []; % delete the invalid rows

%%
% Add the roots and Output the results on the basis of trv_trd
f_output

% format
trvtrd_dh(year, ets)          = {trv_trd}; % Save for next mortality

% save file trv_trd
savefile = strcat('trv_trd',num2str(ets));
savefile = strcat(directory, 'results_', num2str(year),'/', savefile);
save(savefile, 'trv_trd');
