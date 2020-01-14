function theta = calibrate_model_NAG(zzz,flag)

% The parameters in order...

%upper_bound = [2.00, 0.75, 2.00, 0.85, 2.00, 0.75, 0.75,2.00]; 
%lower_bound = [0.01, 0.10, 1.00, 0.10, 1.00, 0.25, 0.25,0.35]; 

must_be_positive = [1,3,5,8];
must_be_zero_one = [2,4,6,7];

params = zeros(8,1);
% params(must_be_positive) = exp(zzz(must_be_positive));
% params(must_be_zero_one) = 1./(1+exp(-zzz(must_be_zero_one)));

params(must_be_positive) = (zzz(must_be_positive));
params(must_be_zero_one) = (zzz(must_be_zero_one));

%params(2) = [];

[yyy] = compute_outcomes_prefshock(params,0);

% The current moments that we are targeting....
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aggregate_moments = [1.89, 0.61, 0.625, 0.47];

%%% Description:
% Wage gap
% The rural share
% The urban variance... note that this is position number 3 (see below)
% Fraction with no liquid assets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Moments...
experiment_moments = [0.22, 0.092, 0.30];

control_moments = [0.36, 0.25, 0.16, 0.10,  0.19];

experiment_hybrid = [0.36, 0.22, 0.092, 0.30, 0.10,  0.40];

experiment_hybrid_v2 = [0.36, 0.22, 0.092, 0.30, 0.10, 0.25,  0.40];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Survey Moments...
%survey_moments = [0.17];

% Offered Average Migrant wage zero unemployment risk, 86 migrate.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targets = [aggregate_moments, experiment_hybrid];

% yyy(3) = yyy(3) + exp(zzz(end-1));
% yyy(end) = yyy(end) + exp(zzz(end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_theta = (targets'-yyy')./targets';
% g_theta = log(targets'./yyy');

yyy([3,end]) = [];
targets([3,end]) = [];

%g_theta = -(targets') + (yyy');

g_theta = log(targets'./yyy');

% const_var = false(length(g_theta),1);
% %This is just saying if the variances are above targets, pennelize. If
% %not...then count as zero penelty. 
% const_var(3) = 1;
% const_var(end) = 1;
% 
% g_theta(const_var) = 0;
% 
% W = eye(length(g_theta));
% 
 if flag == 1

theta = g_theta;
 end
% elseif flag == 2
%     
% theta = g_theta';
% 
% end