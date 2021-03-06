function [data_panel, params] = just_simmulate(params, move, solve_types, assets, vfun, cft_params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_sims = 5000;
time_series = 100000;
N_obs = 25000;

params.N_obs = N_obs;

rng(03281978)

[~, shock_states_p] = hmmgenerate(time_series,params.trans_mat,ones(params.n_shocks));

pref_shocks = rand(time_series,params.n_perm_shocks);
move_shocks = rand(time_series,params.n_perm_shocks);

states_panel = zeros(N_obs,4,params.n_types);

[~,type_weights] = pareto_approx(params.n_perm_shocks, 1./params.perm_shock_u_std);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%params.means_test = 0;

mtest = params.asset_space < params.means_test;
mtest_move = params.m_season.*(~mtest)';
params.m_season = mtest_move;
params.m_fiscal = params.m_season(end) - params.m_season;



sim_panel = zeros(N_obs,15,params.n_types);   
    
parfor xxx = 1:params.n_types 

% Interestingly, this is not a good part of the code to use parfor... it
% runs much faster with just a for loop.
       
    [sim_panel(:,:,xxx), ~] = simmulate_prefshock_welfare(...
        assets(xxx), move(xxx),params, solve_types(xxx,:), params.trans_shocks, shock_states_p, pref_shocks(:,xxx), move_shocks(:,xxx), vfun(xxx));
    % This is the same one as in baseline model
end
    
    
n_draws = floor(N_obs/max(N_obs*type_weights)); % this computes the number of draws.
sample = min(n_draws.*round(N_obs*type_weights),N_obs); % Then the number of guys to pull.
s_count = 1;

for xxx = 1:params.n_types

    e_count = s_count + sample(xxx)-1;
        
    data_panel(s_count:e_count,:) = sim_panel(N_obs-(sample(xxx)-1):end,:,xxx);
    
    s_count = e_count+1;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
income = 1; consumption = 2; assets = 3; live_rural = 4; work_urban = 5;
move = 6; move_season = 7; movingcosts = 8; season = 9; net_asset = 10;
welfare = 11; experince = 12; fiscalcost = 13; tax = 14; production = 15;

rural_not_monga = (data_panel(:,live_rural)==1 & data_panel(:,season)~=1);

if ( isempty(cft_params) || cft_params == 0)
    % if we have not specify, do this. Or if it's zero, then we are again 
    % just trying to replicate the original economy.
    
    params.means_test = (prctile(data_panel(rural_not_monga,3),55) + prctile(data_panel(rural_not_monga,3),45))./2;
    % This is so we can just replicate the old stuff...
else
    params.means_test = cft_params;
    % Here if we are doing the counterfactuall, we want the "same exact
    % guys", policies may change asset distribtuion, so this holds the
    % asset threshold at whatever it was chosen to be. 
end



