function [rural_data, urban_data, params] = compute_outcomes_prefshock_GE(cal_params, wages, cft_params, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver file for the code which is consistent with RR paper at
% Econometrica (late 2017-on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.sigma_nu_not = cal_params(end);
params.sigma_nu_exp = cal_params(end);

params.R = 0.95; % Storage technology that looses value over time. We are thinking currency. Citation for the 0.92 number?

params.beta = 0.95;  params.abar = 0; % Discount factor

params.ubar = cal_params(5);   params.lambda = cal_params(6); params.pi_prob = cal_params(7);

params.rural_options = 3;
params.urban_options = 2;

gamma = 2;
params.A = (1-gamma).^-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shocks...
shock_rho = cal_params(4); % Persistance of shocks
shock_std = cal_params(1).*sqrt((1-shock_rho).^2); % Standard Deviation of shocks

perm_shock_u_std = cal_params(2); % Permenant ability differs in the urban area.

urban_tfp = cal_params(3); rural_tfp = 1./urban_tfp; % Urban TFP

seasonal_factor = 0.55; % The seasonal fluctuation part. 


params.m_season = 0.08; % This is the bus ticket
params.m = 2*params.m_season; % This is the moving cost. 

gamma_urban = cal_params(8); % Gamma parameter (set to 1?)

% m_error_national_survey = 0; % Mesurment error. Set to zero, then expost pick to high variances. 
% m_error_expr = 0;
 
n_perm_shocks = 48; %48
n_tran_shocks = 30; %30
% Number of permenant and transitory types. 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This first sets up the transitory shocks. 
% 
m_adjust_rural = -1./(1-shock_rho.^2).*(shock_std.^2).*(1/2);

m_adjust_urban = -1./(1-shock_rho.^2).*(shock_std.^2).*(1/2).*(gamma_urban);
% This here should set the shocks so that in levels, the mean value is one.
% thus changes in std or gamma only affect variances. 

% [shocks_trans,trans_mat] = tauchen(n_tran_shocks,0,shock_rho,shock_std,2.5);

[shocks_trans,trans_mat] = rouwenhorst(n_tran_shocks,shock_rho,shock_std);

seasonal_trans_mat = [0 , 1 ; 1, 0]; 

if isempty(wages) 
    seasonal_shocks = [log(seasonal_factor); log(1./seasonal_factor)];
else
    seasonal_shocks = [log(wages(1)); log(wages(2))];
end

trans_mat = kron(trans_mat, seasonal_trans_mat);

shocks_trans_temp = repmat(shocks_trans',2,1);

seasonal_shocks_temp = repmat(seasonal_shocks,1,n_tran_shocks);

shocks_trans_r = shocks_trans_temp(:) + seasonal_shocks_temp(:) + m_adjust_rural ;
shocks_trans_u = gamma_urban.*(shocks_trans_temp(:) + m_adjust_urban);

trans_shocks = [exp(shocks_trans_r),exp(shocks_trans_u)];

n_shocks = length(shocks_trans_u);

% Note: Seasonality will only show up in odd periods. Even periods are the
% good season...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This sets up the permenant type of shocks...for now, I'm just going to
% use this tauchen method which with no persistance will correspond with a
% normal distribution and then the transition matrix will determine the
% relative weights of the guys in the population. 

[zurban , zurban_prob] = pareto_approx_alt_v2(n_perm_shocks, 1./perm_shock_u_std);

types = [ones(n_perm_shocks,1), zurban];

type_weights = zurban_prob;

[n_types , ~] = size(types);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up asset space and parameters to pass to the value function
% itteration.
    
params.grid = [100, 0, 3];

asset_space = linspace(params.grid(2),params.grid(3),params.grid(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre generate the shocks
n_sims = 5000;
time_series = 100000;
N_obs = 25000;

rng(03281978)

params.N_obs = N_obs;

rng(03281978)

[~, shock_states_p] = hmmgenerate(time_series,trans_mat,ones(n_shocks));

pref_shocks = rand(time_series,n_perm_shocks);
move_shocks = rand(time_series,n_perm_shocks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the policy functions and then simmulate the time paths. See the
% routines below for details.

% Note depending on the computer you have (and toolboxes with Matlab) using
% the parfor command here does help. It distributes the instructions within
% the for loop across different cores. It this case it leads to a big speed
% up. 

%assets = struct();
% move = zeros(n_types);
% vguess = zeros(n_types);

solve_types = [rural_tfp.*types(:,1), types(:,2)];

if isempty(cft_params) 
    params.means_test = 0;
else
    params.means_test = cft_params;
end

parfor xxx = 1:n_types 
        
    %params = [R, solve_types(xxx,:), beta, m, gamma, abar, ubar, lambda, pi_prob, m_temp];
    [assets(xxx), move(xxx), vguess(xxx)] = ...
        rural_urban_value_prefshock_GE(params, solve_types(xxx,:), trans_shocks, trans_mat);

%     [assets(:,:,:,xxx), move(:,:,:,xxx), vguess(:,:,:,xxx)] = ...
%         rural_urban_value_addit(params, trans_shocks, trans_mat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now simulate the model...

sim_panel = zeros(N_obs,9,n_types);
states_panel = zeros(N_obs,4,n_types);


for xxx = 1:n_types 
% Interestingly, this is not a good part of the code to use parfor... it
% runs much faster with just a for loop.
       
    [sim_panel(:,:,xxx), states_panel(:,:,xxx)] = rural_urban_simmulate_prefshock(...
        assets(xxx), move(xxx),params, solve_types(xxx,:), trans_shocks, shock_states_p, pref_shocks(:,xxx),move_shocks(:,xxx));

end 
% this is how the panel is organized...
% panel = [labor_income, consumption, assets, live_rural, work_urban, move, move_seasn, move_cost, season];

% Now record the data. What we are doing here is creating a
% cross-section/pannel of guys that are taken in porportion to their
% distributed weights. 

n_draws = floor(N_obs/max(N_obs*type_weights)); % this computes the number of draws.
sample = min(n_draws.*round(N_obs*type_weights),N_obs); % Then the number of guys to pull.
s_count = 1;

for xxx = 1:n_types 

    e_count = s_count + sample(xxx)-1;
        
    data_panel(s_count:e_count,:) = sim_panel(N_obs-(sample(xxx)-1):end,:,xxx);
    
    s_count = e_count+1;
   
end

rural_not_monga = data_panel(:,4)==1 & data_panel(:,end)~=1;
%params.means_test = median(data_panel(rural_not_monga,3));


params.means_test = (prctile(data_panel(rural_not_monga,3),55) + prctile(data_panel(rural_not_monga,3),45))./2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This section of the code now performs the expirements. 

% assets_temp = zeros(n_asset_states,n_shocks,2,n_types);
% move_temp = zeros(n_asset_states,n_shocks,2,n_types);
% cons_eqiv = zeros(n_asset_states,n_shocks,2,n_types);

% assets_surv = zeros(n_asset_states,n_shocks,n_types);
% move_surv = zeros(n_asset_states,n_shocks,n_types);
 
sim_expr_panel = zeros(n_sims,13,11,n_types);
sim_cntr_panel = zeros(n_sims,9,11,n_types);
% sim_surv_panel = zeros(n_sims,10,3,n_types);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

periods = 1:length(states_panel(:,:,1))-20;
monga = periods(rem(periods,2)==0)-1;
pref_shocks = pref_shocks((N_obs+1):end,:);
move_shocks = move_shocks((N_obs+1):end,:);

for xxx = 1:n_types     
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, perform the field experiment...

    [assets_temp(xxx), move_temp(xxx), cons_eqiv(xxx)] = field_experiment_welfare_prefshock(params, solve_types(xxx,:), trans_shocks, trans_mat, vguess(xxx));
    
    %[assets_temp(:,:,:,xxx), move_temp(:,:,:,xxx)] = field_experiment(params, trans_shocks, trans_mat, vguess(:,:,:,xxx));

    % This generates an alternative policy function for rural households associated with a
    % the field experiment of paying for a temporary move. The asset_temp
    % provides the asset policy conditional on a temporary move. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    rng(02071983+xxx)
    
    monga_index = monga(randi(length(monga),1,n_sims))';

    [sim_expr_panel(:,:,:,xxx), sim_cntr_panel(:,:,:,xxx)]...
        = experiment_driver_prefshock(assets(xxx), move(xxx), assets_temp(xxx), move_temp(xxx), cons_eqiv(xxx),...
          params, solve_types(xxx,:), trans_shocks, monga_index, states_panel(:,:,xxx), pref_shocks(:,xxx), move_shocks(:,xxx), sim_panel(:,:,xxx));
         
    % This then takes the policy functions, simmulates the model, then
    % after a period of time, implements the experirment.     
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Now the code below constructs a panel so the approriate types are where
% % they should be....

n_draws = floor(n_sims/max(n_sims*type_weights));
sample_expr = min(n_draws.*round(n_sims*type_weights),n_sims);
s_expr_count = 1;

for xxx = 1:n_types
        
    e_expr_count = s_expr_count + sample_expr(xxx)-1;
    
    data_panel_expr(s_expr_count:e_expr_count,:,1) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,1,xxx);

    data_panel_cntr(s_expr_count:e_expr_count,:,1) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,1,xxx);
    
    data_panel_expr(s_expr_count:e_expr_count,:,2) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,2,xxx);

    data_panel_cntr(s_expr_count:e_expr_count,:,2) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,2,xxx);
        
    data_panel_expr(s_expr_count:e_expr_count,:,3) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,3,xxx);
    data_panel_cntr(s_expr_count:e_expr_count,:,3) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,3,xxx);
    
    data_panel_expr(s_expr_count:e_expr_count,:,4) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,4,xxx);
    data_panel_cntr(s_expr_count:e_expr_count,:,4) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,4,xxx);
    
    data_panel_expr(s_expr_count:e_expr_count,:,5) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,5,xxx);
    data_panel_cntr(s_expr_count:e_expr_count,:,5) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,5,xxx);
    
    data_panel_expr(s_expr_count:e_expr_count,:,7) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,7,xxx);
    data_panel_cntr(s_expr_count:e_expr_count,:,7) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,7,xxx);
    
    data_panel_expr(s_expr_count:e_expr_count,:,11) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,11,xxx);
    data_panel_cntr(s_expr_count:e_expr_count,:,11) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,11,xxx);
        
% if flag == 1
%     
%         data_panel_surv(s_expr_count:e_expr_count,:,1) = sim_surv_panel(n_sims-(sample_expr(xxx)-1):end,:,1,xxx);
%         data_panel_surv(s_expr_count:e_expr_count,:,2) = sim_surv_panel(n_sims-(sample_expr(xxx)-1):end,:,2,xxx);
% end
                    
    s_expr_count = e_expr_count + 1;
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rural = data_panel(:,4)==1;
rural_monga = data_panel(:,4)==1 & data_panel(:,end)==1;
rural_not_monga = data_panel(:,4)==1 & data_panel(:,end)~=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part just focuses on the entire sample...

% Income earned by rural residents relative to urban...
m_income = [mean((data_panel(rural,1))), mean((data_panel(~rural,1)))];

% Income in monga relative to non-monga...
m_income_season = [mean((data_panel(rural_monga,1))), mean((data_panel(rural_not_monga,1)))];

% Consumption...
m_consumption = [mean((data_panel(rural,2))), mean((data_panel(~rural,2)))];

% Fraction of residents residing in the rural area...
avg_rural = sum(rural)./length(data_panel);

var_income = [var(log(data_panel(rural,1))), var(log(data_panel(~rural,1)))];
% No emasurment error here, we add it on expost. 

var_consumption = [var(log(data_panel(rural,2))), var(log(data_panel(~rural,2)))];

perm_moves = sum(data_panel(rural,6)==1)./length(data_panel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now use the control and expirement stuff...
% [labor_income, consumption, assets, live_rural, work_urban, move, move_seasn, move_cost, season, experiment_flag];

% First drop people that did not have the experiment performed on them....
rural_cntr = data_panel_cntr(:,4,1)==1 & data_panel_expr(:,end,1)==1;

control_data = data_panel_cntr(rural_cntr,:,:);
expermt_data = data_panel_expr(rural_cntr,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Migration Elasticity

temp_migrate_cntr = control_data(:,7,1) == 1;
temp_migrate_expr = expermt_data(:,7,1) == 1;

temp_migration = sum(temp_migrate_cntr)./sum(rural_cntr);

temp_expr_migration = sum(temp_migrate_expr)./sum(rural_cntr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number_workers_cntr = length(control_data(:,:,1));

number_workers_expr = length(expermt_data(:,:,1));

rural_data.cntr_labor_units = sum((1./rural_tfp).*control_data(~temp_migrate_cntr,1,1))./number_workers_cntr;
% Number of labor units that remain in the control vilage. Note this is
% selected on the means test...need to think about this.

rural_data.expr_labor_units = sum((1./rural_tfp).*expermt_data(~temp_migrate_expr,1,1))./number_workers_expr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This stuff is to figure out the correct scaling parameter.
% panel = [labor_income, consumption, assets, live_rural, work_urban, move, move_seasn, move_cost, season];

labor_units_rural_monga = (data_panel(:,5)~=1 & data_panel(:,end)==1);
% this says, in monga, I'm working in rural area

labor_units_urban_monga = (data_panel(:,5)==1 & data_panel(:,end)==1);
% this says, in monga, I'm working in urban area

labor_units_rural_not_monga = (data_panel(:,5)~=1 & data_panel(:,end)==0);
% this says, NOT monga, I'm working in rural area

labor_units_urban_not_monga = (data_panel(:,5)==1 & data_panel(:,end)==0);
% this says, NOT monga, I'm working in urban area

number_workers_monga = sum(labor_units_rural_monga) + sum(labor_units_urban_monga);
% In the monga, this is the number of guys in total...

number_workers_not_monga = sum(labor_units_rural_not_monga) + sum(labor_units_urban_not_monga);
% outside of monga, number of guys in total...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rural_data.labor_units_monga = sum((1./rural_tfp).*data_panel(labor_units_rural_monga,1))./number_workers_monga;
% the want here is how many effective labor units are in the rural area,
% during the monga. So first sum over all income payments in rural area and
% then remove TFP term (we will put that back in). 
% Then divide through by the total mass of works...just a normalization(?)

% Same idea below...
rural_data.labor_units_not_monga = sum((1./rural_tfp).*data_panel(labor_units_rural_not_monga,1))./number_workers_not_monga;

rural_data.labor_monga = sum(labor_units_rural_monga)./number_workers_monga;

rural_data.labor_not_monga = sum(labor_units_rural_not_monga)./number_workers_not_monga;

rural_data.seasonal_factor = seasonal_factor;


rural_data.labor_units_monga = sum((1./rural_tfp).*data_panel(labor_units_rural_monga,1))./number_workers_monga;
% the want here is how many effective labor units are in the rural area,
% during the monga. So first sum over all income payments in rural area and
% then remove TFP term (we will put that back in). 
% Then divide through by the total mass of works...just a normalization(?)

% Same idea below...
rural_data.labor_units_not_monga = sum((1./rural_tfp).*data_panel(labor_units_rural_not_monga,1))./number_workers_not_monga;

rural_data.labor_monga = sum(labor_units_rural_monga)./number_workers_monga;

rural_data.labor_not_monga = sum(labor_units_rural_not_monga)./number_workers_not_monga;

rural_data.seasonal_factor = seasonal_factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

urban_data.labor_units_monga = sum(data_panel(labor_units_urban_monga,1))./number_workers_monga;

urban_data.labor_units_not_monga = sum(data_panel(labor_units_urban_not_monga,1))./number_workers_not_monga;



if flag == 1
    
    disp('Average Rural Population')
    disp(avg_rural)
    
    disp('Fraction of Rural Who are Migrants')
    disp(temp_migration)
    
end
    
    







