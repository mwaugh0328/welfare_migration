function [move, solve_types, assets, params, vfun, ce] = just_policy(cal_params, wages, vft_fun, cft_params, tax, policyfun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver file for the code which is consistent with RR paper at
% Econometrica (late 2017-on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(tax) 
    
    params.tax.rate = 1;
    params.tax.prog = 0;
else
    params.tax.rate = tax(1);
    params.tax.prog = tax(2);
end

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

params.perm_shock_u_std = cal_params(2); % Permenant ability differs in the urban area.

urban_tfp = cal_params(3); rural_tfp = 1./urban_tfp; % Urban TFP
params.urban_tfp = urban_tfp;
params.rural_tfp = rural_tfp;

seasonal_factor = 0.55; % The seasonal fluctuation part. 
params.seasonal_factor = seasonal_factor;

params.m_season = 0.08; % This is the bus ticket
params.m = 2*params.m_season; % This is the moving cost. 

gamma_urban = cal_params(8); % Gamma parameter (set to 1?)

% m_error_national_survey = 0; % Mesurment error. Set to zero, then expost pick to high variances. 
% m_error_expr = 0;
 
n_perm_shocks = 24; %48
params.n_perm_shocks = n_perm_shocks;
n_tran_shocks = 5; %30
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

% This is if the thing is in wages or in productivity terms.
if isempty(wages) 
    seasonal_shocks = [log(seasonal_factor); log(1./seasonal_factor)];
else
    seasonal_shocks = [log(wages(1)); log(wages(2))];
end

trans_mat = kron(trans_mat, seasonal_trans_mat);
params.trans_mat = trans_mat;

shocks_trans_temp = repmat(shocks_trans',2,1);

seasonal_shocks_temp = repmat(seasonal_shocks,1,n_tran_shocks);

shocks_trans_r = shocks_trans_temp(:) + seasonal_shocks_temp(:) + m_adjust_rural ;
shocks_trans_u = gamma_urban.*(shocks_trans_temp(:) + m_adjust_urban);

trans_shocks = [exp(shocks_trans_r),exp(shocks_trans_u)];
params.trans_shocks = trans_shocks;


n_shocks = length(shocks_trans_u);
params.n_shocks = n_shocks;


% Note: Seasonality will only show up in odd periods. Even periods are the
% good season...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This sets up the permenant type of shocks...for now, I'm just going to
% use this tauchen method which with no persistance will correspond with a
% normal distribution and then the transition matrix will determine the
% relative weights of the guys in the population. 

[zurban , zurban_prob] = pareto_approx(n_perm_shocks, 1./params.perm_shock_u_std);

types = [ones(n_perm_shocks,1), zurban];

type_weights = zurban_prob;

[n_types , ~] = size(types);
params.n_types = n_types;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up asset space and parameters to pass to the value function
% itteration.
    
params.grid = [100, 0, 3];

asset_space = linspace(params.grid(2),params.grid(3),params.grid(1));
params.asset_space = asset_space;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the policy functions and then simmulate the time paths. See the
% routines below for details.

% Note depending on the computer you have (and toolboxes with Matlab) using
% the parfor command here does help. It distributes the instructions within
% the for loop across different cores. It this case it leads to a big speed
% up. 

solve_types = [rural_tfp.*types(:,1), types(:,2)];

% The counterfactual is a means tested moving cost removal. if it's zero,
% this is the baseline model. Otherwise, it;s the means tested value...
if isempty(cft_params) 
    
    params.means_test = 0;
else
    
    params.means_test = cft_params;

end

% Then here is the value function stuff, again depends if there is the
% means test. Also, the final imput can be a value function. This is used
% to compute welfare...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% consumption.rural_not = ones(n_shocks);
% consumption.rural_exp = ones(n_shocks);
% consumption.seasn_not = ones(n_shocks);
% consumption.seasn_exp = ones(n_shocks);
% consumption.urban_new = ones(n_shocks);
% consumption.urban_old = ones(n_shocks);
% 
% move.rural_not = cumsum((1./3).*ones(n_shocks,3),2);
% move.rural_exp = cumsum((1./3).*ones(n_shocks,3),2);
% move.urban_new = cumsum((1./2).*ones(n_shocks,2),2);
% move.urban_old = cumsum((1./2).*ones(n_shocks,2),2);
% 
% [move, vfinal, ~] = effecient_valuefun(consumption, move, params, solve_types(1,:), trans_shocks, trans_mat,[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if isempty(vft_fun) 
    
    if isempty(policyfun)
        
        parfor xxx = 1:n_types 
        
            [assets(xxx), move(xxx), vfun(xxx), ce(xxx)] = ...
            rural_urban_value_prefshock_GE(params, solve_types(xxx,:), trans_shocks, trans_mat, []);
        
        end
    else
        parfor xxx = 1:n_types 
                [assets(xxx), move(xxx), vfun(xxx), ce(xxx)] = ...
            policy_valuefun(policyfun.assets(xxx), policyfun.move(xxx), params, solve_types(xxx,:), trans_shocks, trans_mat, []);
        end
    end
        

else
    if isempty(policyfun)
        
        parfor xxx = 1:n_types 
        
            [assets(xxx), move(xxx), vfun(xxx), ce(xxx)] = ...
            rural_urban_value_prefshock_GE(params, solve_types(xxx,:), trans_shocks, trans_mat, vft_fun(xxx));
        
        end
    else
        parfor xxx = 1:n_types 
                [assets(xxx), move(xxx), vfun(xxx), ce(xxx)] = ...
            policy_valuefun(policyfun.assets(xxx), policyfun.move(xxx), params, solve_types(xxx,:), trans_shocks, trans_mat, vft_fun(xxx));
        end
    end
    
end


%test = 1