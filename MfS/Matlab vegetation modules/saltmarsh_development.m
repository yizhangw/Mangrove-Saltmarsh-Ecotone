if marsh==1
    % Load stress parameter    
    if year==1 && ets==1 && Restart==0
        Relative_flood     = Relative_flood_marsh ;
        Tau_marsh          = Tau0_marsh;
    elseif ets==1
        Relative_flood     = cell2mat(d3dparameters.Flood(year-1).PerYear(t_eco_year,1));
        Tau_marsh          = cell2mat(d3dparameters.Tau(year-1).PerYear(t_eco_year,1));
    else
        Relative_flood     = cell2mat(d3dparameters.Flood(year).PerYear(ets-1,1));
        Tau_marsh          = cell2mat(d3dparameters.Tau(year).PerYear(ets-1,1));
    end

    % Load marsh density
    if ets==1 && year==1 && Restart==0
        PD = zeros(Ndim,Mdim);			
    else
        PD = PD_new;
        clear PD_new
    end
    
    %% Marsh fitness
    Fitness_marsh = a_marsh.*(Relative_flood.^2)+b_marsh.*Relative_flood+c_marsh;
    Fitness_marsh(Fitness_marsh>1) = 1;
    Fitness_marsh(Fitness_marsh<0) = 0;
    %% Establishment
    dt = 1;
    if ets==1 || ets==2 || ets==3 || ets==4 || ets==5 || ets==6 || ets==7    % Seed establishment during the growing season
        dPseed = (rand(Ndim,Mdim)<Seed).*P0.*dt.*(Tau_marsh<tau_maxcrp);
    else
        dPseed = zeros(Ndim,Mdim);
    end
    %% Growth and diffusion
    Diff = (-0.17.*Relative_flood.*24+5.91)./6;   % Diffusion rate links to relative hydroperiod (Zhu et al., 2019)

    if ets==8 || ets==9 || ets==10 || ets==11 || ets==12 % Saltmarshes enter a senescence phase and stop colonizing 
        dPgrowth = zeros(Ndim,Mdim);
        dPdiffx = zeros(Ndim,Mdim);
        dPdiffy = zeros(Ndim,Mdim);
    else
        dPgrowth = (r1.*(1-((PD)./K1)).*(PD).*dt).*Fitness_marsh; % logistic growth of marsh density within a grid cell in one ets
        
        dPdiffx = zeros(Ndim,Mdim);
        for j = 2:1:(Mdim-1)
            dPdiffx(:,j) = ((Diff(:,j+1).*(PD(:,j+1)-2.*PD(:,j)+PD(:,j-1)))./(dx^2)).*dt; % diffusion of plant density in x direction in one ets
        end
        clear j
        dPdiffy = zeros(Ndim,Mdim);
        for h = 2:1:(Ndim-1)
            dPdiffy(h,:) = ((Diff(h+1,:).*(PD(h+1,:)-2.*PD(h,:)+PD(h-1,:)))./(dy^2)).*dt; % diffusion of plant density in y direction in one ets
        end
        clear h
    end

    %% Height growth    
    if ets==1 && year==1 && Restart==0           % Initial height of seeds
        Pnet_growth = dPseed + dPgrowth + dPdiffx + dPdiffy;
        PH_temp = Pnet_growth;
        PH_temp(PH_temp>0) = H0;
        PH = PH_temp.*Fitness_marsh;
        clear PH_temp Pnet_growth
    else                                        % Load marsh height from previous ets
        PH = PH_new;
        clear PH_new
    end

    if ets==8 || ets==9 || ets==10 || ets==11 || ets==12 % Saltmarshes enter a senescence phase and stop growing 
        dPgrowth_height = zeros(Ndim,Mdim);
    else
        dPgrowth_height = (r2.*(1-((PH)./K2)).*(PH).*dt).*Fitness_marsh; % logistic growth of plant height within a grid cell in one ets
    end

    %% Mortality due to flow stress
    Tau_marsh(Tau_marsh<tau_mincrp) = -9999;
    Tau_marsh(Tau_marsh>=tau_maxcrp) = NaN;
    erostau = 1./(tau_maxcrp - tau_mincrp).*(Tau_marsh - tau_mincrp);
    erostau(erostau<0) = 0;
    erostau(isnan(erostau)) = 1;
    dPerostau = erostau.*PD.*dt;
    clear erostau
    
    %% Interspecific competition with mangroves
    if Competition == 1
        Height_mangrove2 = zeros(Ndim,Mdim);
        Biomass_mangrove = zeros(Ndim,Mdim);
        Diameter_mangrove = zeros(Ndim,Mdim);
        trv_trd_temp = unique(trv_trd(:,[1,2]),'rows');
        for i = 1:size(trv_trd_temp,1) % Extract mangrove height, biomass and stem diameter
            kkk  = find(trv_trd(:,1)==trv_trd_temp(i,1) & trv_trd(:,2)==trv_trd_temp(i,2) & trv_trd(:,11)==1);
            mloc = trv_trd(kkk(1),1);
            nloc = trv_trd(kkk(1),2);
            Height_mangrove2(mloc,nloc) = max(trv_trd(kkk,7));
            Biomass_mangrove(mloc,nloc) = sum(trv_trd(kkk,16));
            Diameter_mangrove(mloc,nloc)= max(trv_trd(kkk,13));
        end
        Height_relative = Height_mangrove2 - PH; % relative height between mangrove and saltmarsh
        Height_relative = Height_relative.*(PH>0);
        Height_relative(Height_relative <= 0) = NaN;         

        LW = 38.9.*Diameter_mangrove.^1.62; % mangrove leaf weight (Chen et al., 1998)
        LA = 0.01.*LW; % mangrove leaf area (Chen et al., 1998)
        AL = exp(-0.5.*(Height_relative.*LA)); % saltmarsh available light (Chen et al., 1998)
        rs = 1 - exp(-4.64.*(AL-0.05));
        ri = 2.24*(1-exp(-1.136*(AL-0.08)));
        r_final = (rs + ri)./2; % light competion outcome (Chen et al., 1998)
        r_final(isnan(r_final)) = 1;
        dPCompete_light = r_final;

        fraction_m = Biomass_mangrove./B_half(nv); % mangrove fraction
        fraction_m(fraction_m==0) = NaN;
        fraction_m(fraction_m>1) = 1;
        fraction_s = (4.6.*PH.^2 - 1.17.*PH + 0.73)./1000.*PD./6.2123; % saltmarsh fraction (Zhu et al., 2019)
        fraction_r = fraction_m - fraction_s; % relative fraction between mangrove and saltmarsh
        fraction_r(fraction_r<0) = NaN;
        dPCompete_resource = (1./(1+exp(-12.*(0.5-fraction_r)))); % resource competion outcome
        dPCompete_resource(isnan(dPCompete_resource)) = 1;

        % compare light competion and resource competion and choose the outcome with greater pressure
        Compete_temp = dPCompete_light < dPCompete_resource;
        Compete_light = Compete_temp.*dPCompete_light;
        Compete_resource = (Compete_temp==0).*dPCompete_resource;
        dPCompete = (Compete_light + Compete_resource).*(PH>0);

    else
        dPCompete_light    = ones(Ndim,Mdim);
        dPCompete_resource = ones(Ndim,Mdim);
        dPCompete          = ones(Ndim,Mdim);
    end
    clear Height_relative LW LA AL rs ri r_final Biomass_mangrove Diameter_mangrove...
        fraction_m fraction_s fraction_r Compete_temp Compete_light...
        Compete_resource i mloc nloc trv_trd_temp kkk

    %% Predtion
    if Predation == 1
        if ets==9 || ets==10
            PD_pred = (PH<1).*PD; % Kill short seedlings
            Ppred_temp = 0.83.*(Relative_flood.^2) - 1.42.*Relative_flood + 1;
            Ppred_temp(Ppred_temp<0.4) = 0.4;
            Ppred_temp(Ppred_temp>0.75) = 0.75;
            Ppred_P = zeros(Ndim,Mdim);
            for i = 1:Ndim
                for j = 1:Mdim
                    Ppred_P_temp = rand(1)<Ppred_temp(i,j);
                    Ppred_P(i,j) = Ppred_P_temp;
                    clear Ppred_P_temp
                end
            end
            dPpred = Ppred_P.*0.8.*(PD_pred>0).*PD.*dt;
        else
            dPpred = zeros(Ndim,Mdim);
        end
        clear PD_pred Ppred_temp Ppred_P i j 
    else
        dPpred = zeros(Ndim,Mdim);
    end
    %% Net development of density and height
    dPnet = PD + dPseed + dPgrowth.*dPCompete + dPdiffx + dPdiffy;
    dPnet_temp = dPnet - (dPerostau.*(dPerostau>0)) - (dPpred.*(dPpred>0));
    dPnet_temp = round(dPnet_temp,0);
    dPnet_temp(dPnet_temp<0) = 0;

    dHnet = PH + (dPgrowth_height).*dPCompete;
    dHnet_temp = dHnet.*(dPnet_temp>0);
    dheight_diff_temp = dPnet_temp;
    dheight_diff_temp(dheight_diff_temp>0) = 999;
    dheight_diff_temp1 = dHnet_temp;
    dheight_diff_temp1(dheight_diff_temp1>0) = 998;
    dheight_diff_temp2 = dheight_diff_temp - dheight_diff_temp1;
    dheight_diff_temp2(dheight_diff_temp2==1) = 0;
    dheight_diff_temp2(dheight_diff_temp2==999) = 1;
    diff_height = dheight_diff_temp2.*H0;    
    dPnet_height = dHnet_temp + diff_height;
    
    clear dHnet_temp dheight_diff_temp dheight_diff_temp1 dheight_diff_temp2 
    %% Decearse of density and hight during the non-growing season
    if ets==8 ||ets==9 || ets==10 || ets==11 || ets==12
        Ptemp = dPnet_temp.*0.88;
    else
        Ptemp = dPnet_temp;
    end

    if ets==8 ||ets==9 || ets==10 || ets==11 || ets==12
        Ptemp2 = dPnet_height.*0.83;
    else
        Ptemp2 = dPnet_height;
    end

    clear dPnet_temp dPnet_height
    %% Convert to a format that d3d can read
    
    PD_new = [Ptemp(2,1:Mdim);[Ptemp(2:Ndim-1,2),Ptemp(2:Ndim-1,2:Mdim-1),Ptemp(2:Ndim-1,Mdim-1)];Ptemp(Ndim-1,1:Mdim)];
    PH_new = [Ptemp2(2,1:Mdim);[Ptemp2(2:Ndim-1,2),Ptemp2(2:Ndim-1,2:Mdim-1),Ptemp2(2:Ndim-1,Mdim-1)];Ptemp2(Ndim-1,1:Mdim)];

    clear Ptemp Ptemp2    
    %% Update ecoengineering effects   
    
    vege_2_trachy                  = (4.6.*PH_new.^2 - 1.17.*PH_new + 0.73)./1000.*PD_new./6.2123;
    vege_2_trachy(vege_2_trachy>1) = 1;
    [row_n,col_m,fract_s]          = find(vege_2_trachy);
    
    [~,~,mD]    = find(PD_new);

    [row_n2,col_m2,mH]    = find(PH_new);
    H_vege                = [row_n2 col_m2 mH];  
    [Row2,~]              = size(H_vege);
    if Row2>0
        h_vege = ones(size(row_n,1),1).*mH;  
    else
        h_vege = zeros(size(row_n,1),1);
    end
    d_veg   = 0.005;                                         % Stem diameter
    m_veg   = d_veg.*mD;
    Cd      = ones(size(row_n,1),1);                         % Drag coefficient
    Cb      = ones(size(row_n,1),1).*65;                     % Roughness= 65
    veg_no  = 1:1:size(row_n,1);
    Form_no = ones(size(row_n,1),1).*154;
    veg_all = [row_n col_m veg_no' fract_s Form_no h_vege m_veg Cd Cb];
    veg_all(:,6)      = round(veg_all(:,6),2);
    veg_all(:,7)      = round(veg_all(:,7),4);
    trach_inp_aruv    = veg_all(:,1:4);
    veg_all_temp      = veg_all(:,4:9);
    veg_all_temp(:,1) = 1:1:size(veg_all_temp,1);
    TrNo_mark = 1;    
    while ~isempty(veg_all_temp)
        veg_temp                        = veg_all_temp(1,3:4);
        TrNo_Loc                        = find(veg_all_temp(:,3)==veg_temp(1) & veg_all_temp(:,4)==veg_temp(2));
        trach_inp_def(TrNo_mark,1:6)    = veg_all_temp(TrNo_Loc(1),1:6);
        trach_inp_def(TrNo_mark,1)      = TrNo_mark + 10000; % Differentiate from mangroves
        trach_inp_aruv(veg_all_temp(TrNo_Loc,1),3) = TrNo_mark + 10000;
        TrNo_mark                       = TrNo_mark+1;
        veg_all_temp(TrNo_Loc,:)        = [];
    end 
    %% Save marsh results
    Saltmarsh.PD            = PD_new;
    Saltmarsh.PH            = PH_new;
    Saltmarsh.D             = Diff;
    Saltmarsh.Tau           = Tau_marsh;
    Saltmarsh.Seed          = dPseed;
    Saltmarsh.Growth        = dPgrowth;
    Saltmarsh.Growth_height = dPgrowth_height;
    Saltmarsh.Diffx         = dPdiffx;
    Saltmarsh.Diffy         = dPdiffy;
    Saltmarsh.Diff_height   = diff_height;
    Saltmarsh.Net           = dPnet;
    Saltmarsh.Re_f          = Relative_flood;
    Saltmarsh.Fitness       = Fitness_marsh;
    Saltmarsh.Erostau       = dPerostau;
    Saltmarsh.Pred          = dPpred;
    Saltmarsh.Light         = dPCompete_light;
    Saltmarsh.Resource      = dPCompete_resource;
    Saltmarsh.Compete       = dPCompete;

    save(strcat(directory,'work/Saltmarsh.mat'),'Saltmarsh')
    copyfile(strcat(directory, 'work/Saltmarsh.mat'), strcat(directory, 'results_', num2str(year), '/Saltmarsh_', num2str(ets),'.mat'));

    clear length Row Row2 Column Column2 row_n row_n2 col_m col_m2 mD mH fract_s Vege_no Area_f Form_no h_vege Cd Cb vege_2_trachy vege_2_trachy2...
        botdep1 Tau_marsh Relative_flood PD PH Diff botdep_next_run Fitness_marsh dPdiffx dPdiffy dPerostau...
        dPpred dPgrowth dPseed dPnet Saltmarsh dPgrowth_height diff_height P_vege...
        H_vege dPCompete dPpred dPCompete_light dPCompete_resource
end