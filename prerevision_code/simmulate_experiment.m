function [panel_expr] = simmulate_experiment...
    (assets_policy, move_policy, assets_temp, move_temp, cons_eqiv, params,...
    state_at_expr, trans_shocks, shock_states, pref_shock)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up grid for asset holdings. 

grid = params(1:3);
params = params(4:end);

means_test = params(7);

n_asset_states = grid(1);

asset_space = linspace(grid(2),grid(3),grid(1));
% asset_space = [0, logspace(log10(grid(2)),log10(grid(3)),n_asset_states-1)];

R = params(1); 
z_rural = params(2); z_urban = params(3); W = params(4);
% These are the permanent shocks. 

m = params(4);
m_seasn = params(5);
lambda = params(6);
pi_prob = params(8);

% These control unemployment risk. 

r_shocks = trans_shocks(:,1); 
u_shocks = trans_shocks(:,2);
expr_shock = pref_shock;

n_shocks = length(r_shocks);

policy_assets_rural_nxpr = assets_policy(:,:,1);
policy_assets_rural_expr = assets_policy(:,:,2);

policy_assets_seasn_nxpr = assets_policy(:,:,3);
policy_assets_seasn_expr = assets_policy(:,:,4);

policy_assets_urban_new = assets_policy(:,:,5);
policy_assets_urban_old = assets_policy(:,:,6);

policy_move_rural_nxpr = move_policy(:,:,1);
policy_move_rural_expr = move_policy(:,:,2);
policy_move_urban_new = move_policy(:,:,3);
policy_move_urban_old = move_policy(:,:,4);

assets_rural_temp_not = assets_temp(:,:,1);
assets_rural_temp_exp = assets_temp(:,:,2);

move_rural_temp_not = move_temp(:,:,1);
move_rural_temp_exp = move_temp(:,:,2);

cons_eqiv_not = cons_eqiv(:,:,1);
cons_eqiv_exp = cons_eqiv(:,:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_series = length(shock_states);

% The initial location state...
location = zeros(time_series+1,1); 

move = zeros(time_series,1);
move_seasn = zeros(time_series,1);
move_cost = zeros(time_series,1);

experiment_flag = zeros(time_series,1);
welfare = zeros(time_series,1);

labor_income = zeros(time_series,1);

consumption = zeros(time_series,1);
assets = zeros(time_series+1,1);
season = zeros(time_series,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the way it will work, we will run this out some time say 1000
% periods. The implement the experiment. Then run it out 10 more periods to
% see what happens. 

asset_state_at_expr = state_at_expr(1);
asset_state_p = asset_state_at_expr;


location(1) = state_at_expr(2);
experiment_flag(1) = 1;
 
assets(1,1) = asset_space(asset_state_at_expr);

xxx = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Urban located guys are not in the sample....
if location(xxx) == 2 || location(xxx) == 4 || location(xxx) == 5 || location(xxx) == 6
    
location = location(1:end-1,1);
assets = assets(1:end-1,1); 
% This is important. So this now generates assets
% % held at date t (not chosen). So one can construce the budget constraint.

live_rural = location == 1 | location == 2 | location == 3 | location == 4;
work_urban = location == 2 | location == 4 | location == 5 | location == 6;
experince = location == 3;

experiment_flag(1) = 0;
welfare(1) = 0;

    season(xxx,1) = mod(shock_states(xxx),2); 
    urban_skill = z_urban.*ones(length(labor_income),1);
    panel_expr = [labor_income, consumption, assets, live_rural, work_urban,...
        move, move_seasn, move_cost, season, welfare, experiment_flag];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Rural guys that have too much wealth, not in the sample.
if assets(1,1) > means_test

location = location(1:end-1,1);
assets = assets(1:end-1,1); 
% This is important. So this now generates assets
% % held at date t (not chosen). So one can construce the budget constraint.

live_rural = location == 1 | location == 2 | location == 3 | location == 4;
work_urban = location == 2 | location == 4 | location == 5 | location == 6;
experince = location == 3;

experiment_flag(1) = 0;
welfare(1) = 0;
    
    season(xxx,1) = mod(shock_states(xxx),2); 
    urban_skill = z_urban.*ones(length(labor_income),1);
    panel_expr = [labor_income, consumption, assets, live_rural, work_urban,...
        move, move_seasn, move_cost, season, welfare, experiment_flag];
    return
end    


experiment_flag(1) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This means perform the experiment....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if location(xxx) == 1 % in rural area
    
        move(xxx,1) = move_rural_temp_not(asset_state_at_expr,shock_states(xxx))== 3;
        move_seasn(xxx,1) = move_rural_temp_not(asset_state_at_expr,shock_states(xxx))== 2;
                
        welfare(xxx,1) = cons_eqiv_not(asset_state_at_expr,shock_states(xxx));
        
        labor_income(xxx,1) = z_rural.*r_shocks(shock_states(xxx));        
        asset_state_p = assets_rural_temp_not(asset_state_at_expr,shock_states(xxx));
        
        location(xxx+1) = location(xxx);
                    
        if move_seasn(xxx,1) == 1
            location(xxx+1) = 2;
            move_cost(xxx,1) = 0 ;
        elseif move(xxx,1) == 1
            location(xxx+1) = 3;
            move_cost(xxx,1) = m;            
        end

end

if location(xxx) == 3 % in rural area
    
        move(xxx,1) = move_rural_temp_exp(asset_state_at_expr,shock_states(xxx))== 3;
        move_seasn(xxx,1) = move_rural_temp_exp(asset_state_at_expr,shock_states(xxx))== 2;
        
        welfare(xxx,1) = cons_eqiv_exp(asset_state_at_expr,shock_states(xxx));
                
        labor_income(xxx,1) = z_rural.*r_shocks(shock_states(xxx));        
        asset_state_p = assets_rural_temp_exp(asset_state_at_expr,shock_states(xxx));
        
        location(xxx+1) = location(xxx);
                    
        if expr_shock(xxx) < (1-pi_prob);
            if move_seasn(xxx,1) == 1 
                location(xxx+1,1) = 2;  
                move_cost(xxx,1) = 0;
            elseif move(xxx,1) == 1
                location(xxx+1) = 5;
                move_cost(xxx,1) = m;
            end
        else 
            if move_seasn(xxx,1) == 1 
                location(xxx+1,1) = 4; % 
                move_cost(xxx,1) = 0;
            elseif move(xxx,1) == 1
                location(xxx+1) = 6;
                move_cost(xxx,1) = m;
            end
        end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asset_state = asset_state_p;
    
assets(xxx+1,1) = asset_space(asset_state_p);
        
consumption(xxx,1) = labor_income(xxx,1) + R.*assets(xxx,1) - assets(xxx+1,1) - move_cost(xxx,1);

season(xxx,1) = mod(shock_states(xxx),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for xxx = 2:length(shock_states)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Work through stuff conditional on a location....
    %
    %                       NO-EXPERIENCE
    %
    if location(xxx) == 1 % rural area
        
        move(xxx,1) = policy_move_rural_nxpr(asset_state,shock_states(xxx))== 3;
        move_seasn(xxx,1) = policy_move_rural_nxpr(asset_state,shock_states(xxx))== 2;
                
        labor_income(xxx,1) = z_rural.*r_shocks(shock_states(xxx));        
        asset_state_p = policy_assets_rural_nxpr(asset_state,shock_states(xxx));
        
        location(xxx+1) = location(xxx);
                    
        if move_seasn(xxx,1) == 1
            location(xxx+1) = 2;
            move_cost(xxx,1) = m_seasn;
        elseif move(xxx,1) == 1
            location(xxx+1) = 5;
            move_cost(xxx,1) = m;            
        end
        
    elseif location(xxx) == 2 % seasonal movers....
       
        labor_income(xxx,1) = z_urban.*u_shocks(shock_states(xxx));
 
        asset_state_p = policy_assets_seasn_nxpr(asset_state,shock_states(xxx));
                
       if pref_shock(xxx) < (1-lambda);
            location(xxx+1,1) = 3; % get experince
       else 
            location(xxx+1,1) = 1; 
       end
               
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                       EXPERIENCE
    elseif location(xxx) == 3
        
        move(xxx,1) = policy_move_rural_expr(asset_state,shock_states(xxx))== 3;
        move_seasn(xxx,1) = policy_move_rural_expr(asset_state,shock_states(xxx))== 2;
                
        labor_income(xxx,1) = z_rural.*r_shocks(shock_states(xxx));        
        asset_state_p = policy_assets_rural_expr(asset_state,shock_states(xxx));
        
        location(xxx+1) = location(xxx);
        
        if expr_shock(xxx) < (1-pi_prob);
            if move_seasn(xxx,1) == 1 
                location(xxx+1,1) = 2;  
                move_cost(xxx,1) = m_seasn;
            elseif move(xxx,1) == 1
                location(xxx+1) = 5;
                move_cost(xxx,1) = m;
            end
        else 
            if move_seasn(xxx,1) == 1 
                location(xxx+1,1) = 4; % 
                move_cost(xxx,1) = m_seasn;
            elseif move(xxx,1) == 1
                location(xxx+1) = 6;
                move_cost(xxx,1) = m;
            end
        end
        
    elseif location(xxx) == 4
        
        labor_income(xxx,1) = z_urban.*u_shocks(shock_states(xxx));
 
        asset_state_p = policy_assets_seasn_expr(asset_state,shock_states(xxx));
        
        location(xxx+1) = 3;
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    elseif location(xxx) == 5
        
        move(xxx,1) = policy_move_urban_new(asset_state,shock_states(xxx))== 2;
        
        labor_income(xxx,1) = z_urban.*u_shocks(shock_states(xxx));
 
        asset_state_p = policy_assets_urban_new(asset_state,shock_states(xxx));
        
        location(xxx+1) = location(xxx);
        
        if move(xxx,1) == 1 
            location(xxx+1,1) = 1; % Return to being rural...
            move_cost(xxx,1) = m;
        elseif pref_shock(xxx) < (1-lambda);
            location(xxx+1,1) = 6; % Lose the aversion to urban area
        end
        
    elseif location(xxx) == 6
        
        move(xxx,1) = policy_move_urban_old(asset_state,shock_states(xxx))== 2;
        
        labor_income(xxx,1) = z_urban.*u_shocks(shock_states(xxx));
 
        asset_state_p = policy_assets_urban_old(asset_state,shock_states(xxx));
        
        location(xxx+1) = location(xxx);
        
        if move(xxx,1) == 1 
            location(xxx+1,1) = 3; % Return to being rural...
            move_cost(xxx,1) = m;
        end      
    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    asset_state = asset_state_p;
    
    assets(xxx+1,1) = asset_space(asset_state_p);
            
    consumption(xxx,1) = labor_income(xxx,1) + R.*assets(xxx,1) - assets(xxx+1,1) - move_cost(xxx,1);
        
    season(xxx,1) = mod(shock_states(xxx),2);     
end

location = location(1:end-1,1);
assets = assets(1:end-1,1); 
% This is important. So this now generates assets
% % held at date t (not chosen). So one can construce the budget constraint.

live_rural = location == 1 | location == 2 | location == 3 | location == 4;
work_urban = location == 2 | location == 4 | location == 5 | location == 6;
experince = location == 3;

urban_skill = z_urban.*ones(length(labor_income),1);

% panel_expr = [labor_income, consumption, assets, live_rural, work_urban, move,...
%     move_seasn, move_cost, season, welfare, urban_skill, experince, experiment_flag];

panel_expr = [labor_income, consumption, assets, live_rural, work_urban, move,...
    move_seasn, move_cost, season, welfare, experiment_flag];
