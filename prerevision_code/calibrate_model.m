function theta = calibrate_model(zzz,flag)

% The parameters in order...

params = zzz;

[yyy] = compute_outcomes(params,0);

% The current moments that we are targeting....
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aggregate Moments....
aggregate_moments = [1.80, 0.63, 0.68, 0.47];

% Wage gap of 2.62
% 63 percent reside in rural area
% Aggregate variance of log income (urban)
% Fraction of rural households with no liquid assets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Moments...
experiment_moments = [0.22, 0.092, 0.30];

control_moments = [0.36, 0.25, 0.16, 0.10,  0.40];

experiment_hybrid = [0.36, 0.22, 0.092, 0.30, 0.10,  0.40];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Survey Moments...
%survey_moments = [0.17];

% Offered Average Migrant wage zero unemployment risk, 86 migrate.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targets = [aggregate_moments, experiment_hybrid];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_theta = (targets'-yyy')./targets';
% g_theta = log(targets'./yyy');

g_theta = (targets')-(yyy');

const_var = false(length(g_theta),1);
% This is just saying if the variances are above targets, pennelize. If
% not...then count as zero penelty. 
const_var(3) = g_theta(3) > 0;
const_var(end) = g_theta(end) > 0;

g_theta(const_var) = 0;

W = eye(length(g_theta));

if flag == 1

theta = g_theta'*W*g_theta;

elseif flag == 2
    
theta = g_theta';

end