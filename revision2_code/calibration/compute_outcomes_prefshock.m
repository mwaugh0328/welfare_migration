function [targets] = compute_outcomes_prefshock(cal_params, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver file for the code which is consistent with RR2 paper at
% Econometrica (late 2020-on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.rural_options = 3;
params.urban_options = 2;

%Preferences
params.sigma_nu_not = cal_params(9); %These are the logit shocks
params.sigma_nu_exp = cal_params(9);

params.R = cal_params(12);  %The storage technology

params.beta = cal_params(13);  % The discount factor

params.abar = cal_params(14); % The discount factor

params.ubar = cal_params(5);   % ubar, disutility of being in urban area

params.lambda = cal_params(6); % getting experince and losing it

params.pi_prob = cal_params(7);

params.pref_gamma = 2; % Riskaversion (need to look at paper and change name)

params.A = (1-params.pref_gamma).^-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shocks...

shock_rho = cal_params(4); % Persistance of shocks

shock_std = cal_params(1).*sqrt((1-shock_rho).^2); % Standard Deviation of shocks

perm_shock_u_std = cal_params(2); % Permenant ability differs in the urban area.

urban_tfp = cal_params(3); 

rural_tfp = 1./urban_tfp; % Urban TFP

params.seasonal_factor = cal_params(11); % The seasonal fluctuation part. 

params.m_season = cal_params(10); % This is the bus ticket

params.m = 2*params.m_season; % This is the moving cost. 

gamma_urban = cal_params(8); % Gamma parameter (set to 1?)

% Number of permenant and transitory types. 
n_perm_shocks = 36; %48
n_tran_shocks = 15; %30


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This first sets up the transitory shocks. 

m_adjust_rural = -1./(1-shock_rho.^2).*(shock_std.^2).*(1/2);

m_adjust_urban = -1./(1-shock_rho.^2).*(shock_std.^2).*(1/2).*(gamma_urban);
% This here should set the shocks so that in levels, the mean value is one.
% thus changes in std or gamma only affect variances. 

%[shocks_trans,trans_mat] = tauchen(n_tran_shocks,0,shock_rho,shock_std,2.5);

[shocks_trans,trans_mat_temp] = rouwenhorst(n_tran_shocks,shock_rho,shock_std);

seasonal_trans_mat = [0 , 1 ; 1, 0]; 

seasonal_shocks = [log(params.seasonal_factor); log(1./params.seasonal_factor)];

params.trans_mat = kron(trans_mat_temp, seasonal_trans_mat);

shocks_trans_temp = repmat(shocks_trans',2,1);

seasonal_shocks_temp = repmat(seasonal_shocks,1,n_tran_shocks);

shocks_trans_r = shocks_trans_temp(:) + seasonal_shocks_temp(:) + m_adjust_rural ;
shocks_trans_u = gamma_urban.*(shocks_trans_temp(:) + m_adjust_urban);

params.trans_shocks = [exp(shocks_trans_r),exp(shocks_trans_u)];

params.n_shocks = length(shocks_trans_u);

% Note: Seasonality will only show up in odd periods. Even periods are the
% good season...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This sets up the permenant type of shocks...for now, I'm just going to
% use this tauchen method which with no persistance will correspond with a
% normal distribution and then the transition matrix will determine the
% relative weights of the guys in the population. 

[zurban , zurban_prob] = pareto_approx(n_perm_shocks, 1./perm_shock_u_std);

types = [ones(n_perm_shocks,1), zurban];

type_weights = zurban_prob;

[n_types , ~] = size(types);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up asset space and parameters to pass to the value function itteration.
% This is the new grid setup. It places a finer grid near the constraint
% and moving cost...

%params.grid = [100, 0, 2];

grid1 = [30, 0.0, 0.29];
grid2 = [70, 0.30, 2];

params.asset_space = [linspace(grid1(2),grid1(3),grid1(1)), linspace(grid2(2),grid2(3),grid2(1))];

%params.asset_space = linspace(params.grid(2),params.grid(3),params.grid(1));
asset_space = params.asset_space;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre generate the shocks
n_sims = 10000;
time_series = 100000;
N_obs = 25000;

params.N_obs = N_obs;

rng(03281978)

[~, shock_states_p] = hmmgenerate(time_series,params.trans_mat,ones(params.n_shocks));

pref_shocks = rand(time_series,n_perm_shocks);
move_shocks = rand(time_series,n_perm_shocks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the policy functions and then simmulate the time paths. See the
% routines below for details.

% Note depending on the computer you have (and toolboxes with Matlab) using
% the parfor command here does help. It distributes the instructions within
% the for loop across different cores. It this case it leads to a big speed
% up. 

solve_types = [rural_tfp.*types(:,1), types(:,2)];

parfor xxx = 1:n_types 

    [assets(xxx), move(xxx), vguess(xxx)] = ...
        rural_urban_value_prefshock(params, solve_types(xxx,:));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now simulate the model...

sim_panel = zeros(N_obs,9,n_types);
states_panel = zeros(N_obs,4,n_types);

parfor xxx = 1:n_types 
% Interestingly, this is not a good part of the code to use parfor... it
% runs much faster with just a for loop.
       
    [sim_panel(:,:,xxx), states_panel(:,:,xxx)] = rural_urban_simmulate_prefshock(...
        assets(xxx), move(xxx), params, solve_types(xxx,:), shock_states_p, pref_shocks(:,xxx),move_shocks(:,xxx));

end 

% Now record the data. What we are doing here is creating a
% cross-section/pannel of guys that are taken in porportion to their
% distributed weights. 

n_draws = floor(N_obs/max(N_obs*type_weights)); % this computes the number of draws.
sample = min(n_draws.*round(N_obs*type_weights),N_obs); % Then the number of guys to pull.
s_count = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This then pulls from the pannel
for xxx = 1:n_types 

    e_count = s_count + sample(xxx)-1;
        
    data_panel(s_count:e_count,:) = sim_panel(N_obs-(sample(xxx)-1):end,:,xxx);
    
    s_count = e_count+1;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEre is the means test for Mushfiq's expr...the latter is to smooth
% things

rural_not_monga = data_panel(:,4)==1 & data_panel(:,end)~=1;
%params.means_test = median(data_panel(rural_not_monga,3));


params.means_test = (prctile(data_panel(rural_not_monga,3),55) + prctile(data_panel(rural_not_monga,3),45))./2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This section of the code now performs the expirements. 
 
sim_expr_panel = zeros(n_sims,13,11,n_types);
sim_cntr_panel = zeros(n_sims,9,11,n_types);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

periods = 1:length(states_panel(:,:,1))-20;
monga = periods(rem(periods,2)==0)-1;
pref_shocks = pref_shocks((N_obs+1):end,:);
move_shocks = move_shocks((N_obs+1):end,:);

parfor xxx = 1:n_types     
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, perform the field experiment...

    [assets_temp(xxx), move_temp(xxx), cons_eqiv(xxx)] = field_experiment_welfare_prefshock(params, solve_types(xxx,:), vguess(xxx));
    
    %[assets_temp(:,:,:,xxx), move_temp(:,:,:,xxx)] = field_experiment(params, trans_shocks, trans_mat, vguess(:,:,:,xxx));

    % This generates an alternative policy function for rural households associated with a
    % the field experiment of paying for a temporary move. The asset_temp
    % provides the asset policy conditional on a temporary move. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    rng(02071983+xxx)
    
    monga_index = monga(randi(length(monga),1,n_sims))';

    [sim_expr_panel(:,:,:,xxx), sim_cntr_panel(:,:,:,xxx)]...
        = experiment_driver_prefshock(assets(xxx), move(xxx), assets_temp(xxx), move_temp(xxx), cons_eqiv(xxx),...
          params, solve_types(xxx,:), monga_index, states_panel(:,:,xxx), pref_shocks(:,xxx), move_shocks(:,xxx), sim_panel(:,:,xxx));
         
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

exp_index = [1,2,3,4,5,7,11];

for xxx = 1:n_types
        
    e_expr_count = s_expr_count + sample_expr(xxx)-1;
    
    for zzz = 1:length(exp_index)
    
        data_panel_expr(s_expr_count:e_expr_count,:,exp_index(zzz)) = sim_expr_panel(n_sims-(sample_expr(xxx)-1):end,:,exp_index(zzz),xxx);

        data_panel_cntr(s_expr_count:e_expr_count,:,exp_index(zzz)) = sim_cntr_panel(n_sims-(sample_expr(xxx)-1):end,:,exp_index(zzz),xxx);
        
    end
                            
    s_expr_count = e_expr_count + 1;
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%frac_no_assets = sum(control_data(:,3,1) < asset_space(2))./sum(rural_cntr);

frac_no_assets = 0.95*(sum(control_data(:,3,1) == asset_space(1)))/sum(rural_cntr) + 0.05*(sum(control_data(:,3,1) == asset_space(2)))/sum(rural_cntr);
% Trying to smmoth this thing out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aggregate_moments = [m_income(2)./m_income(1), avg_rural, var_income(2), frac_no_assets];

experiment_moments = [migration_elasticity, migration_elasticity_y2, LATE];

control_moments = [temp_migration, control_migration_cont_y2, control_migration_cont_y3, OLS, var_consumption_no_migrate_control];

experiment_hybrid = [temp_migration, migration_elasticity, migration_elasticity_y2, LATE, OLS, params.m_season./mean(AVG_C), var_cons_growth];

experiment_hybrid_v2 = [temp_migration, migration_elasticity, migration_elasticity_y2, LATE, OLS,...
    control_migration_cont_y2./temp_migration, var_cons_growth];


targets = [aggregate_moments, experiment_hybrid] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
if flag == 1
    
disp('')
disp('')
disp('Average Rural Population')
disp(avg_rural)
disp('Temporary Moving Cost Relative to Mean Consumption')
disp(params.m_season./mean(AVG_C))
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
disp([mean((data_panel(rural_not_monga,1))),mean((data_panel(rural_monga,1)))])

% cd('..\Analysis')
% 
% m_rates = [migration_elasticity, migration_elasticity_y2, NaN, migration_elasticity_y3, NaN, migration_elasticity_y5];
% m_rates = 100.*m_rates';
% 
% plot_migration
% 
% cd('..\calibration')

figure

subplot(3,2,1), hist(log(data_panel(rural,1)),50)
 
subplot(3,2,2), hist(log(data_panel(~rural,1)),50)

subplot(3,2,3), hist(log(data_panel(rural,2)),50)
 
subplot(3,2,4), hist(log(data_panel(~rural,2)),50)
 
subplot(3,2,5), hist((data_panel(rural,3)),50)
 
subplot(3,2,6), hist((data_panel(~rural,3)),50)
    


end









