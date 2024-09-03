%% On the basis of {1N| 2M| 3trachNo| 4Areafractioin| 5trachids| 6rougheq| 7h(m)| 8dens(1/m)| 9Cd| 10Cz| 11vegtype| 12vegnum| 13vegdia(cm)}
% and the relative hydroperiod, P (could be every ets OR yearly)
% we can calculate the values in the 14th-20th col: 14IndS| 15SingleW| 16MultW| 17ComS| 18I*C| 19MortMark| 20MatrixNo

%%
for nv = 1:num_veg_types
    if ismember(nv,trv_trd_addIC(:,11))
        % Inundation stress
        nv_tr_Loc       = find(trv_trd_addIC(:,11)==nv); % vegetation type location
        nv_P_Loc        = sub2ind(size(P),trv_trd_addIC(nv_tr_Loc,1),trv_trd_addIC(nv_tr_Loc,2)); % convert to 1d coordinate
        trv_trd_addIC(nv_tr_Loc,14) = a(nv)*P(nv_P_Loc).^2+b(nv)*P(nv_P_Loc)+c(nv); % add inundation stress in the 14th col based on last ets results
        out_node        = nv_tr_Loc(P(nv_P_Loc)<xL(nv) | P(nv_P_Loc)>xR(nv)); % find the relative hydroperiod value out of the vegetation scale
        trv_trd_addIC(out_node,14)  = 0; % when relative inundation exceeds the scale, set to 0
        
        % Competition stress
        W_tree_a        = bio_a(nv)*(trv_trd_addIC(nv_tr_Loc,13).^2.*trv_trd_addIC(nv_tr_Loc,7)).^ind_a(nv); % aboveground tree weight, Unit-diamater is cm
        W_tree_b        = bio_b(nv)*(trv_trd_addIC(nv_tr_Loc,13).^2.*trv_trd_addIC(nv_tr_Loc,7)).^ind_b(nv); % belowground tree weight, Unit-weight(kg/tree)
        trv_trd_addIC(nv_tr_Loc,15) = (W_tree_a + W_tree_b); % Single tree weight, Unit-weight(kg/tree)
        trv_trd_addIC(nv_tr_Loc,16) = trv_trd_addIC(nv_tr_Loc,15).*trv_trd_addIC(nv_tr_Loc,12); % add plant weight in the 16th col
    end
end
clear nv_tr_Loc nv_P_Loc out_node W_tree_a W_tree_b
trv_trd_addIC_temp = trv_trd_addIC(:,14);
trv_trd_addIC_temp(trv_trd_addIC_temp>1) = 1;
trv_trd_addIC(:,14) = trv_trd_addIC_temp;
clear trv_trd_addIC_temp

IC_rc_temp = unique(trv_trd_addIC(:,[1,2]),'rows'); % Extract rows and columns and appear only once
if marsh==1 && Competition==1
    if ets==1 && year==1 && Restart==0
        PH_temp = zeros(Ndim,Mdim);
        PD_temp = zeros(Ndim,Mdim);
    else
        PD_temp = PD_new; 
        PH_temp = PH_new;
    end
   
    ecotone_marsh_height=[];
    ecotone_marsh_density=[];
    for mmm=1:size(IC_rc_temp,1)
        ecotone_marsh_height=[ecotone_marsh_height;PH_temp(IC_rc_temp(mmm,1),IC_rc_temp(mmm,2))];
        ecotone_marsh_density=[ecotone_marsh_density;PD_temp(IC_rc_temp(mmm,1),IC_rc_temp(mmm,2))];
    end
    clear mmm
    ecotone_marsh_density(ecotone_marsh_density>0 & ecotone_marsh_density<1) = 1;
 
    % C_inter for light
    for iii = 1:size(trv_trd_addIC,1) % Loop over vegetation cells
        if trv_trd_addIC(iii,11)==1
            mloc = trv_trd_addIC(iii,1);
            nloc = trv_trd_addIC(iii,2);
            Height_marsh = PH_temp(mloc,nloc);
            Density_marsh = PD_temp(mloc,nloc);
            Height_mangrove = trv_trd_addIC(iii,7);
            if Height_mangrove > Height_marsh
                Competition_inter_light = 1;
            else
                Height_relative = Height_marsh - Height_mangrove;
                LA = (0.1263.*((Height_marsh.*100).^1.5899))./10000.*Density_marsh; 
                AL = exp(-0.1.*(Height_relative.*LA));
                rs = 1 - exp(-4.64.*(AL-0.05));
                ri = 2.24*(1-exp(-1.136*(AL-0.08)));
                r_final = (rs + ri)./2;
                Competition_inter_light = r_final;
            end
            trv_trd_addIC(iii,21) = Competition_inter_light;
            clear mloc nloc Height_marsh Density_marsh Height_mangrove...
                rs ri r_final
        end
    end
    clear iii PH_temp PD_temp

 
    for iii = 1:size(IC_rc_temp,1) 
        kkk                  = find(trv_trd_addIC(:,1)==IC_rc_temp(iii,1) & trv_trd_addIC(:,2)==IC_rc_temp(iii,2));
        Bio_total            = sum(trv_trd_addIC(kkk,16));
        Competition_intra    = 1./(1+exp(d(nv)*(B_half(nv)-Bio_total))); % C_intra
        trv_trd_addIC(kkk,23) = Competition_intra;
		Height_marsh         = ecotone_marsh_height(iii);  
        Density_marsh        = ecotone_marsh_density(iii); 
    % C_inter for resource
        if ets==8 ||ets==9 || ets==10 || ets==11 || ets==12
            Competition_inter_resource = 1; % Ignoring the impact of dead saltmarshes
        else
            fraction_m = Bio_total./B_half(nv);
            fraction_m(fraction_m>1) = 1;
            fraction_s = (4.6.*Height_marsh.^2 - 1.17.*Height_marsh + 0.73)./1000.*Density_marsh./6.2123;
            if fraction_m > fraction_s
                Competition_inter_resource = 1;
            else
                fraction_r = fraction_s - fraction_m;
                C_resource = (1./(1+exp(-12.*(0.5-fraction_r))));
                Competition_inter_resource = C_resource;
                clear fraction_r C_resource
            end
        end
        mloc = IC_rc_temp(iii,1);
        nloc = IC_rc_temp(iii,2);
        jjj = trv_trd_addIC(:,1)==mloc & trv_trd_addIC(:,2)==nloc;
        trv_trd_addIC(jjj,22) = Competition_inter_resource;
    end
    clear iii kkk Bio_total Competition_intra Density_marsh...
        fraction_m fraction_s Competition_inter_resource mloc nloc jjj...
        ecotone_marsh_height ecotone_marsh_density

    trv_trd_addIC(:,17) = min(trv_trd_addIC(:,21),trv_trd_addIC(:,22)).*trv_trd_addIC(:,23);

else
    % No interspecific competition, only intraspecific competition between mangroves
    IC_rc_temp   = unique(trv_trd_addIC(:,[1,2]),'rows'); % Extract rows and columns and appear only once
    for iii = 1:size(IC_rc_temp,1) % Loop over vegetation cells
        kkk                   = find(trv_trd_addIC(:,1)==IC_rc_temp(iii,1) & trv_trd_addIC(:,2)==IC_rc_temp(iii,2));
        Bio_total             = sum(trv_trd_addIC(kkk,16));
        Competition_stress    = 1./(1+exp(d(nv)*(B_half(nv)-Bio_total))); % B_half here is a mean value
        trv_trd_addIC(kkk,17) = Competition_stress;
        trv_trd_addIC(kkk,23) = Competition_stress;
        clear kkk Bio_total Competition_stress
    end
end
trv_trd_addIC(:,18) = trv_trd_addIC(:,14).*trv_trd_addIC(:,17); % last I* present C for particular plant at particular cell, mortality will further use this value
cell_empty          = find(trv_trd_addIC(:,12)==0);
trv_trd_addIC(trv_trd_addIC(:,12)==0,:) = [];% remove the rows where their num. of veg = 0
trv_trd_addIC(:,20) = 1:1:size(trv_trd_addIC,1); % set sequence for each row
clear iii cell_empty IC_rc_temp