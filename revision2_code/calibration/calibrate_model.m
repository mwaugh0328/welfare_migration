function theta = calibrate_model(cal_params,specs,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cal_params should have the following order
% 1: standard Deviation of shocks (todo, veryfy its stand dev or variance)
% 2: Pareto shape parameter for permenant ability in the urban area.
% 3: Urban TFP
% 4: Persistance of transitory shocks
% 5: Ubar, disutility of being in urban area
% 6: Getting experince 
% 7: Losing it. (TO DO, veryify the 6 and 7 is correct rel. paper)
% 8: Gamma parameter in shock process
% 9: Logit shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then there are some pre-determined, hand calibrated values...

cal_params(10) = 0.08; % this is the moving cost. It was hand calibrated to 
% deliver it being approx 10% of half yearly consumption. Did try internal
% calibration...does provide better fit actually.

cal_params(11) = 0.55; % This is the seasonal shock. Again, hand calibrated
% to approx match the large drop in income and consumption in Monga from
% the Khandeker paper.

cal_params(12) = 0.95; % This is the gross real interest rate on the riskfree
%asset/storage technology. Consistent with their primary asset being cash and
% an annual 10% inflation rate

cal_params(13) = 0.95; % This is the discount factor. Would be something worth 
%trying to calibratate in the future. 

cal_params(14) = 0.0; % This is the abar...it's set to zero. We are able to 
% to have it be positive. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(specs)
    grid1 = [30, 0.0, 0.29]; % this is done to overwieght grid at constraint and
                             % moving cost.
    grid2 = [70, 0.30, 2];

    specs.asset_space = [linspace(grid1(2),grid1(3),grid1(1)), linspace(grid2(2),grid2(3),grid2(1))];
    %specs.asset_space = linspace(0,2,100); % this is the equally spaced grid

    specs.n_perm_shocks = 24;
    specs.n_trans_shocks = 15;

    specs.time_series = 100000; % length of the time series for each perm type
    specs.N_obs = 25000; % grab last number of observations
    specs.n_sims = 10000; % given the pannel above how many times to sample for experiment
end

%[yyy] = compute_outcomes(cal_params, specs,1);

[yyy] = analyze_outcomes(cal_params, specs,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aggregate Moments....
aggregate_moments = [1.89, 0.61, 0.625, 0.47];

%%% Description:
% Wage gap
% The rural share
% The urban variance... note that this is position number 3 (see below)
% Fraction with no liquid assets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Moments...
%experiment_moments = [0.22, 0.092, 0.30];

%control_moments = [0.36, 0.25, 0.16, 0.10,  0.19];

%experiment_hybrid_v2 = [0.36, 0.22, 0.092, 0.30, 0.10, 0.25/0.36, 0.40];

experiment_hybrid = [0.36, 0.22, 0.092, 0.30, 0.10, 0.40];

% The experiment hybrid is a combination of conrol and experiment...
% seasonal migration in control
% increase in r1 (22 percent)
% increase in r2 (9.2 percent)
% LATE estiamte
% OLS estimate
% Standard deviation of consumption growth. See line 400 in
% ``compute_outcomes.``

% Note there is currently an inconsistency between the numbers in the table
% and what I have here. 0.40 corresponds with a variance of 0.16, note 0.18
% which is what is in my slides. To edit: onece we have a consistent
% number, we should change code so we work in only variances or stds.

% Also note how this is working, in ``compute_outcomes'' we do not include
% measurment error in the out put moments. The idea is that since it is
% additive, we only need to check ex-post what the measurment error is.
% This simplifies the calibration since we are calibrating 2 less
% parameters. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targets = [aggregate_moments, experiment_hybrid];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_theta = (targets'-yyy')./targets';

yyy([3,end]) = [];
targets([3,end]) = [];

g_theta = zeros(length(targets),1);

g_theta(1) = log(targets(1)'./yyy(1)');

g_theta(2:end) = (targets(2:end)')-(yyy(2:end)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is how the variances are being treated. 
% const_var = false(length(g_theta),1);
% 
% % This is just saying if the variances are above targets, pennelize. If
% % not...then count as zero penelty. 
% const_var(3) = g_theta(3) > 0;
% const_var(end) = g_theta(end) > 0;
% 
% g_theta(const_var) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = eye(length(g_theta));

if flag == 1

theta = g_theta'*W*g_theta;

elseif flag == 2
    
theta = yyy;

end