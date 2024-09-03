%% To create trv and trd file
%>> When trv_trd is finally set up at the end of each ets
%>> Every column will be checked again here
%:: and then
%:: trv and trd will be further created for the vegetation simulation


% 09/10/2018 exclude null trv_trd
if isempty(trv_trd)
    display(['None Vegetation appears at the Year ' num2str(year) ' ETS ' num2str(ets) ]);
else
    trv_trd      = sortrows(trv_trd,1:2); % sortrows  
    
    % 13/08/2018 consider roots
    fract_m           = trv_trd(:,16)./B_half;
    fract_m(fract_m>1) = 1;
    trv_trd(:,4)      = fract_m; % update the Area Fraction after mortality and colonization

    trv_trd(:,8)      = trv_trd(:,12)./S_cell.*trv_trd(:,13)/100; % update the density (1/m)
    
    % Sortrows
    trv_trd           = sortrows(trv_trd,1:2); % sortrows
    trv_trd(:,20)     = 1:1:size(trv_trd,1); % update matrix sequence    
    
    %% Output and format file
    % Re-number the col-'TrachytopeNr' both in trd and trv
    % remove the same value in TRD
    trv_trd(:,7)      = round(trv_trd(:,7),2); % height, reserve a decimal fraction
    trv_trd(:,8)      = round(trv_trd(:,8),4); % density, reserve 4 decimal fraction2
    trv_txt           = trv_trd(:,1:4);
    trd_txt_temp      = trv_trd(:,5:10);
    trd_txt_temp(:,7) = 1:1:size(trd_txt_temp,1); % Sequence mark
    TrNo_mark         = 1;
    while ~isempty(trd_txt_temp)
        a_temp                    = trd_txt_temp(1,3:5); % store h,n,cd to a new matrix
        TrNo_Loc                  = find(trd_txt_temp(:,3)==a_temp(1) & trd_txt_temp(:,4)==a_temp(2) & trd_txt_temp(:,5)==a_temp(3)); % find the same trd parameters
        trd_txt(TrNo_mark,1:6)    = trd_txt_temp(TrNo_Loc(1),1:6); % choose one as a representative of TRD
        trd_txt(TrNo_mark,1)      = TrNo_mark; % modify the TrNo in TRD
        trv_txt(trd_txt_temp(TrNo_Loc,7),3) = TrNo_mark; % change the TrNo in TRV as well
        TrNo_mark                 = TrNo_mark+1;
        trd_txt_temp(TrNo_Loc,:)  = [];
    end
    clear a_temp TrNo_Loc TrNo_mark
    % Write TRD
    if marsh==1
        trd_txt2=[trd_txt;trach_inp_def];
    else
        trd_txt2=trd_txt;
    end
    dlmwrite(strcat(directory, 'work\veg','.trd'),trd_txt2, '\t'); % write trd file to folder
    
    % Write TRV
    try % if there is no vegetation (error is created) just copy file otherwise sort data
        if marsh==1
            trv_txt2=[trv_txt;trach_inp_aruv];
        else
            trv_txt2=trv_txt;
        end
        % sort the vegetation in trv-file
        trv = sortrows(trv_txt2,1:2);
        % update fraction areas of all vegetation types based on changed trv file
        dlmwrite(strcat(directory, 'work\veg','.trv'),trv, '\t');
    catch
        disp('veg.trv is empty');
    end
end