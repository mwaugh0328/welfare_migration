function [targets] = compute_outcomes(params, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver file for the code which is consistent with RR paper at
% Econometrica (late 2017-on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 2; 

R = 0.95; % Storage technology that looses value over time. We are thinking currency. Citation for the 0.92 number?

beta = 0.95;  abar = 0; % Discount factor

ubar = params(5);   lambda = params(6); pi_prob = params(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shocks...
shock_rho = params(4); % Persistance of shocks
shock_std = params(1).*sqrt((1-shock_rho).^2); % Standard Deviation of shocks

perm_shock_u_std = params(2); % Permenant ability differs in the urban area.

urban_tfp = params(3); rural_tfp = 1./urban_tfp; % Urban TFP

seasonal_factor = 0.55; % The seasonal fluctuation part. 

m = 0.16; % This is the moving cost. 
m_temp = 0.08; % This is the bus ticket

gamma_urban = params(8); % Gamma parameter (set to 1?)

% m_error_national_survey = 0; % Mesurment error. Set to zero, then expost pick to high variances. 
% m_error_expr = 0;
 
n_perm_shocks = 48; 
n_tran_shocks = 30;
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

seasonal_shocks = [log(seasonal_factor); log(1./seasonal_factor)];

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
    
grid = [100, 0, 6];
n_asset_states = grid(1);

asset_space = linspace(grid(2),grid(3),grid(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre generate the shocks
n_sims = 3000;
time_series = 50000;
N_obs = 25000;

rng(03281978)

[~, shock_states_p] = hmmgenerate(time_series,trans_mat,ones(n_shocks));

pref_shocks = rand(time_series,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the policy functions and then simmulate the time paths. See the
% routines below for details.

% Note depending on the computer you have (and toolboxes with Matlab) using
% the parfor command here does help. It distributes the instructions within
% the for loop across different cores. It this case it leads to a big speed
% up. 

assets = zeros(n_asset_states,n_shocks,6,n_types);
move = zeros(n_asset_states,n_shocks,4,n_types);

solve_types = [rural_tfp.*types(:,1), types(:,2)];

parfor xxx = 1:n_types 
    
    params = [R, solve_types(xxx,:), beta, m, gamma, abar, ubar, lambda, pi_prob, m_temp];

    [assets(:,:,:,xxx), move(:,:,:,xxx), vguess(:,:,:,xxx)] = ...
        rural_urban_value_mex(params, trans_shocks, trans_mat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now simulate the model...

sim_panel = zeros(N_obs,9,n_types);
states_panel = zeros(N_obs,4,n_types);

for xxx = 1:n_types 
% Interestingly, this is not a good part of the code to use parfor... it
% runs much faster with just a for loop.
    
    params = [R, solve_types(xxx,:),  m, m_temp, lambda, pi_prob];
    
    [sim_panel(:,:,xxx), states_panel(:,:,xxx)] = rural_urban_simmulate(assets(:,:,:,xxx), move(:,:,:,xxx),...
        grid, params, N_obs, trans_shocks, shock_states_p, pref_shocks');
    
%     [sim_panel(:,:,xxx), states_panel(:,:,xxx)] = rural_urban_simmulate_mex_p(assets(:,:,:,xxx), move(:,:,:,xxx),...
%         grid, params, N_obs, trans_shocks, shock_states_p, pref_shocks',trans_mat);
% 
%     [sim_panel(:,:,xxx), states_panel(:,:,xxx)] = rural_urban_simmulate_plot(assets(:,:,:,xxx), move(:,:,:,xxx),...
%         grid, params, N_obs, trans_shocks, shock_states_p, pref_shocks, trans_mat);
%   
end 

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
means_test = median(data_panel(rural_not_monga,3));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This section of the code now performs the expirements. 

assets_temp = zeros(n_asset_states,n_shocks,2,n_types);
move_temp = zeros(n_asset_states,n_shocks,2,n_types);
cons_eqiv = zeros(n_asset_states,n_shocks,2,n_types);

% assets_surv = zeros(n_asset_states,n_shocks,n_types);
% move_surv = zeros(n_asset_states,n_shocks,n_types);
 
sim_expr_panel = zeros(n_sims,11,11,n_types);
sim_cntr_panel = zeros(n_sims,9,11,n_types);
% sim_surv_panel = zeros(n_sims,10,3,n_types);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

periods = 1:length(states_panel(:,:,1))-20;
monga = periods(rem(periods,2)==0)-1;
pref_shocks = pref_shocks((N_obs+1):end,1);

for xxx = 1:n_types     
    
    params = [R, solve_types(xxx,:), beta, m, gamma, abar, ubar, lambda, m_temp, pi_prob];
    
    params_sim = [grid, R, solve_types(xxx,:),  m, m_temp, lambda, means_test, pi_prob];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, perform the field experiment...

    [assets_temp(:,:,:,xxx), move_temp(:,:,:,xxx)] = field_experiment_welfare(grid, params, trans_shocks, trans_mat, vguess(:,:,:,xxx));
    
    %[assets_temp(:,:,:,xxx), move_temp(:,:,:,xxx)] = field_experiment(params, trans_shocks, trans_mat, vguess(:,:,:,xxx));

    % This generates an alternative policy function for rural households associated with a
    % the field experiment of paying for a temporary move. The asset_temp
    % provides the asset policy conditional on a temporary move. 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the survey... tic

%     [assets_surv(:,:,xxx), move_surv(:,:,xxx)] = survey_results(grid, params, trans_shocks, trans_mat, vguess(:,:,:,xxx));
   
     rng(02071983+xxx)
    monga_index = monga(randi(length(monga),1,n_sims))';

    [sim_expr_panel(:,:,:,xxx), sim_cntr_panel(:,:,:,xxx)]...
        = experiment_driver(assets(:,:,:,xxx), move(:,:,:,xxx), ...
          assets_temp(:,:,:,xxx), move_temp(:,:,:,xxx), cons_eqiv(:,:,:,xxx),...
          params_sim, trans_shocks, monga_index, states_panel(:,:,xxx), pref_shocks, sim_panel(:,:,xxx));
         
    % This then takes the policy functions, simmulates the model, then
    % after a period of time, implements the experirment. 
    
% if flag == 1
%     [assets_surv(:,:,xxx), move_surv(:,:,xxx)] = survey_results(grid, params, trans_shocks, trans_mat, vguess(:,:,:,xxx));
%     
%     [sim_surv_panel(:,:,:,xxx), sim_cntr_panel(:,:,:,xxx)]...
%         = experiment_driver_mex(assets(:,:,:,xxx), move(:,:,:,xxx), ...
%           assets_surv(:,:,xxx), move_surv(:,:,xxx),...
%           params_sim, trans_shocks, monga_index, states_panel(:,:,xxx), pref_shocks, sim_panel(:,:,xxx));
% end
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data_panel(:,1) = exp(log(data_panel(:,1)) + m_error_national_survey.*randn(n_obs_panel,1));
% data_panel(:,2) = exp(log(data_panel(:,2)) + m_error_national_survey.*randn(n_obs_panel,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% panel = [labor_income, consumption, assets, live_rural, work_urban, move, move_seasn, move_cost, expected_urban, season];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First devine some indicator variables...

rural = data_panel(:,4)==1;
rural_monga = data_panel(:,4)==1 & data_panel(:,end)==1;
rural_not_monga = data_panel(:,4)==1 & data_panel(:,end)~=1;

% urban_monga = data_panel(:,4)~=1 & data_panel(:,end)==1;
% urban_not_monga = data_panel(:,4)~=1 & data_panel(:,end)~=1;

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

migration_elasticity = temp_expr_migration - temp_migration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Migration Elasticity Year 2

temp_migrate_cntr_y2 = control_data(:,7,3) == 1;
temp_migrate_expr_y2 = expermt_data(:,7,3) == 1;

temp_migration_y2 = sum(temp_migrate_cntr_y2)./sum(rural_cntr);

temp_expr_migration_y2 = sum(temp_migrate_expr_y2)./sum(rural_cntr);

migration_elasticity_y2 = temp_expr_migration_y2 - temp_migration_y2;

cont_y2 = control_data(:,7,1) == 1 & control_data(:,7,3) == 1;
control_migration_cont_y2 = sum(cont_y2)./sum(rural_cntr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Migration Elasticity Year 4

temp_migrate_cntr_y3 = control_data(:,7,7) == 1;
temp_migrate_expr_y3 = expermt_data(:,7,7) == 1;

temp_migration_y3 = sum(temp_migrate_cntr_y3)./sum(rural_cntr);

temp_expr_migration_y3 = sum(temp_migrate_expr_y3)./sum(rural_cntr);

migration_elasticity_y3 = temp_expr_migration_y3 - temp_migration_y3;

cont_y3 = control_data(:,7,1) == 1 & control_data(:,7,3) == 1 & control_data(:,7,7) == 1;
control_migration_cont_y3 = sum(cont_y3)./sum(rural_cntr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Migration Elasticity Year 5

temp_migrate_cntr_y5 = control_data(:,7,11) == 1;
temp_migrate_expr_y5 = expermt_data(:,7,11) == 1;

temp_migration_y5 = sum(temp_migrate_cntr_y5)./sum(rural_cntr);

temp_expr_migration_y5 = sum(temp_migrate_expr_y5)./sum(rural_cntr);

migration_elasticity_y5 = temp_expr_migration_y5 - temp_migration_y5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the LATE estimate in the model, the same way as in the
% data...so first stack stuff the way we want it....

all_migration = [temp_migrate_cntr; temp_migrate_expr];

not_control = [zeros(length(temp_migrate_cntr),1); ones(length(temp_migrate_expr),1)];

first_stage_b = regress(all_migration, [ones(length(not_control),1), not_control]);

predic_migration = first_stage_b(1) + first_stage_b(2).*not_control;

consumption_noerror = [control_data(:,2,2); expermt_data(:,2,2)];

c_noerror_no_migrate = control_data(~temp_migrate_cntr,2,2);

AVG_C = mean(consumption_noerror);

OLS_beta = regress(consumption_noerror, [ones(length(predic_migration),1), all_migration]);
OLS = OLS_beta(2)./AVG_C;

LATE_beta = regress(consumption_noerror, [ones(length(predic_migration),1), predic_migration]);
LATE = LATE_beta(2)./AVG_C ;

var_consumption_no_migrate_control = var(log(c_noerror_no_migrate));

cons_drop = mean(log(control_data(~temp_migrate_cntr,2,1))-log(control_data(~temp_migrate_cntr,2,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cons_data_no_error_r1 = [control_data(:,2,1); expermt_data(:,2,1)];
cons_data_no_error_r2 = [control_data(:,2,2); expermt_data(:,2,2)];

% cons_data_r1 = exp(log(cons_data_no_error_r1) + m_error.*randn(length(cons_data_no_error_r1),1)); 
% cons_data_r2 = exp(log(cons_data_no_error_r2) + m_error.*randn(length(cons_data_no_error_r2),1));
% 
% cons_model = [ [temp_migrate_cntr; temp_migrate_expr], [zeros(length(temp_migrate_cntr),1); ...
%                 ones(length(temp_migrate_expr),1)], log(cons_data_r1), log(cons_data_r2)];
            
cons_model_growth = log(cons_data_no_error_r1) - log(cons_data_no_error_r2);
var_cons_growth = std(cons_model_growth);
% Again, no measurment error here, we can add it on expost. Need
% consistency in language. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assets...
frac_no_assets = sum(control_data(:,3,1) < asset_space(2))./sum(rural_cntr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aggregate_moments = [m_income(2)./m_income(1), avg_rural, var_income(2), frac_no_assets];

experiment_moments = [migration_elasticity, migration_elasticity_y2, LATE];

control_moments = [temp_migration, control_migration_cont_y2, control_migration_cont_y3, OLS, var_consumption_no_migrate_control];

experiment_hybrid = [temp_migration, migration_elasticity, migration_elasticity_y2, LATE, OLS, var_cons_growth];


targets = [aggregate_moments, experiment_hybrid] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
if flag == 1
    
disp('')
disp('')
disp('Average Rural Population')
disp(avg_rural)
disp('Temporary Moving Cost Relative to Mean Consumption')
disp(m_temp./mean(AVG_C))
disp('Fraction of Rural Who are Migrants')
disp(temp_migration)
disp('Expr Elasticity: Year One, Two, Four')
disp([migration_elasticity, migration_elasticity_y2, migration_elasticity_y3, migration_elasticity_y5])
disp('Control: Year One, Repeat Two, Four')
disp([temp_migration, control_migration_cont_y2 , control_migration_cont_y3])
disp('OLS Estimate')
disp(OLS)
disp('LATE Estimate')
disp(LATE)
disp('Wage Gap')
disp(m_income(2)./m_income(1))
disp('Mean Consumption Rural and Urban')
disp(m_consumption)
disp('Variance of Consumption Rural and Urban')
disp(var_consumption)
disp('Variance of Consumption Growth')
disp(var_cons_growth)
disp('Variance of Log Income Rural and Urban')
disp(var_income)
disp('Fraction of Rural with No Assets')
disp(frac_no_assets)
disp('Permenant Moves')
disp(perm_moves)
disp('Ratio of Income in Monga vs Non-Monga')
disp(m_income_season(1)/m_income_season(2))
disp('Consumption Drop')
disp(cons_drop)

cd('.\Analysis')

m_rates = [migration_elasticity, migration_elasticity_y2, NaN, migration_elasticity_y3, NaN, migration_elasticity_y5];
m_rates = 100.*m_rates';

plot_migration

cd('..')

figure

subplot(3,2,1), hist(log(data_panel(rural,1)),50)
 
subplot(3,2,2), hist(log(data_panel(~rural,1)),50)

subplot(3,2,3), hist(log(data_panel(rural,2)),50)
 
subplot(3,2,4), hist(log(data_panel(~rural,2)),50)
 
subplot(3,2,5), hist((data_panel(rural,3)),50)
 
subplot(3,2,6), hist((data_panel(~rural,3)),50)
    


end









