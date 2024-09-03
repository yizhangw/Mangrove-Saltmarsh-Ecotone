% Colonisation: determine the available cells that can be colonized according to strategies

clear SeedLoc
% check colonisation strategy of vegetation type
if general_veg_char(1,3,1) == 1  % Colonisation strategy 1; 
    if ets==1 || ets==2 || ets==3 || ets==4 || ets==5 || ets==6 || ets==7
        ColonisationStrategy_dh % check habitat area and evaluate the growth possibility by checking I and C
    else
    end
    
elseif general_veg_char(1,3,1) == 2 % Colonistion strategy 2;
    %% unpublised version
end