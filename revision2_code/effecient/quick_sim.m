function [data_panel] = quick_sim(data_panel, state_panel, vfun, muc, cons_policy)

consumption = 1; live_rural = 2; work_urban = 3; move = 4;
move_season  = 5; movingcosts = 6; season = 7; welfare = 8; experince = 9; production = 10;
maringal_utility = 11;

%[location, season, shock_states']

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cons_policy_rural_not = reshape([cons_policy(:).rural_not],10,24);
vfun_rural_not = reshape([vfun(:).rural_not],10,24);
muc_rural_not= reshape([muc(:).rural_not],10,24);


cons_policy_rural_exp = reshape([cons_policy(:).rural_exp],10,24);
vfun_rural_exp = reshape([vfun(:).rural_exp],10,24);
muc_rural_exp= reshape([muc(:).rural_exp],10,24);

cons_policy_seasn_not = reshape([cons_policy(:).seasn_not],10,24);
vfun_seasn_not = reshape([vfun(:).seasn_not],10,24);
muc_seasn_not= reshape([muc(:).seasn_not],10,24);

cons_policy_seasn_exp = reshape([cons_policy(:).seasn_exp],10,24);
vfun_seasn_exp = reshape([vfun(:).seasn_exp],10,24);
muc_seasn_exp= reshape([muc(:).seasn_exp],10,24);

cons_policy_urban_new= reshape([cons_policy(:).urban_new],10,24);
vfun_urban_new = reshape([vfun(:).urban_new],10,24);
muc_urban_new= reshape([muc(:).urban_new],10,24);

cons_policy_urban_old = reshape([cons_policy(:).urban_old],10,24);
vfun_urban_old = reshape([vfun(:).urban_old],10,24);
muc_urban_old= reshape([muc(:).urban_old],10,24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for zzz = 1:length(data_panel)

    if state_panel(zzz,1) == 1
                
        data_panel(zzz,consumption) = cons_policy_rural_not(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,welfare) = vfun_rural_not(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,maringal_utility) = muc_rural_not(state_panel(zzz,3),state_panel(zzz,4));
        
        continue
    end
    
    if state_panel(zzz,1) == 6

        data_panel(zzz,consumption) = cons_policy_urban_old(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,welfare) = vfun_urban_old(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,maringal_utility) = muc_urban_old(state_panel(zzz,3),state_panel(zzz,4));
        
        continue
    end
    
    if state_panel(zzz,1) == 3
        data_panel(zzz,consumption) = cons_policy_rural_exp(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,welfare) = vfun_rural_exp(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,maringal_utility) = muc_rural_exp(state_panel(zzz,3),state_panel(zzz,4));
        
        continue
    end
    
    if state_panel(zzz,1) == 2

        data_panel(zzz,consumption) = cons_policy_seasn_not(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,welfare) = vfun_seasn_not(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,maringal_utility) = muc_seasn_not(state_panel(zzz,3),state_panel(zzz,4));
        
        continue
    end

    if state_panel(zzz,1) == 4
   
        data_panel(zzz,consumption) = cons_policy_seasn_exp(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,welfare) = vfun_seasn_exp(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,maringal_utility) = muc_seasn_exp(state_panel(zzz,3),state_panel(zzz,4));
        
        continue
    end
    
    if state_panel(zzz,1) == 5
   
        data_panel(zzz,consumption) = cons_policy_urban_new(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,welfare) = vfun_urban_new(state_panel(zzz,3),state_panel(zzz,4));
        
        data_panel(zzz,maringal_utility) = muc_urban_new(state_panel(zzz,3),state_panel(zzz,4));
        
        continue
    end
    

    
end

        