function [panel, states] = rural_urban_simmulate(assets_policy, move_policy, grid,...
                    params, N_obs, trans_shocks, shock_states, pref_shock)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simmulates a time series/cross section of variables that we can map
% to data. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up grid for asset holdings. 

n_asset_states = grid(1);

asset_space = linspace(grid(2),grid(3),grid(1));
% asset_space = [0, logspace(log10(grid(2)),log10(grid(3)),n_asset_states-1)];

R = params(1); 
z_rural = params(2); z_urban = params(3); 
% These are the permanent shocks. 

m = params(4);
m_seasn = params(5);
lambda = params(6);
pi = params(7);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_series = length(shock_states);

shock_states = shock_states(1:end);

asset_state = 10; 
% The initial asset state....this should not matter as long as everything
% is ergodic.
asset_state_p = asset_state;

% The initial location state...
location = zeros(time_series+1,1); 
location(1) = 1;
% This will start the person in rural area, so 2 = urban, 1 = rural. This
% should not matter either. 

move = zeros(time_series,1);
move_seasn = zeros(time_series,1);
move_cost = zeros(time_series,1);

labor_income = zeros(time_series,1);
compr_advntg = zeros(time_series,1);

consumption = zeros(time_series,1);
assets = zeros(time_series+1,1);
season = zeros(time_series,1);

assets(1,1) = asset_space(asset_state);

rec_asset_states = zeros(time_series+1,1); 
rec_asset_states(1,1) = asset_state;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin simmulation...

for xxx = 1:time_series
    
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
        
        if expr_shock(xxx) < (1-pi);
            if move_seasn(xxx,1) == 1 
                location(xxx+1,1) = 2; % 
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
    
    rec_asset_states(xxx+1,1) = asset_state_p;
        
    consumption(xxx,1) = labor_income(xxx,1) + R.*assets(xxx,1) - assets(xxx+1,1) - move_cost(xxx,1);
    
    compr_advntg(xxx,1) = (z_urban.*u_shocks(shock_states(xxx)))./(z_rural.*r_shocks(shock_states(xxx)));
        
    season(xxx,1) = mod(shock_states(xxx),2);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Record the stuff....
% trans_shocks(shock_states,:); This is a way to see how transitory shocks
% matter. 
location = location(1:end-1,1);
assets = assets(1:end-1,1); 
rec_asset_states = rec_asset_states(1:end-1,1);

% This is important. So this now generates assets
% held at date t (not chosen). So one can construce the budget constraint.

live_rural = location == 1 | location == 2 | location == 3 | location == 4;
work_urban = location == 2 | location == 4 | location == 5 | location == 6;

panel = [labor_income, consumption, assets, live_rural, work_urban, move, move_seasn, move_cost, season];

states = [rec_asset_states, location, season, shock_states'];

panel = panel(time_series-(N_obs-1):end,:);
states = states(time_series-(N_obs-1):end,:);


