clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guess = [0.8029    0.4860    1.4901    0.5292    1.5365    0.4229    0.4611];

upper_bound = [2.00, 0.75, 2.00, 0.85, 2.00, 2.00]; 
lower_bound = [0.01, 0.10, 0.50, 0.10, 1.00, 0.35]; 

% upper_bound_t = [1.50, 0.60, 2.00, 0.85, 1.70, 0.75, 1.00]; 
% lower_bound_t = [0.60, 0.40, 1.30, 0.40, 1.30, 0.25, 0.35]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ObjectiveFunction = @(x) calibrate_model((x),1);

options = gaoptimset('Display','iter','Generations',100,'Initialpopulation',[]);

cal = ga(ObjectiveFunction, length(upper_bound),[],[],[],[],(lower_bound),(upper_bound),[],options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options_pa = optimoptions('patternsearch','Display','iter','MaxFunEvals',2e3);

guess = cal;

new_cal = patternsearch(ObjectiveFunction,guess,[],[],[],[],(lower_bound),(upper_bound),[],options_pa) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save calibration_allmoments_new new_cal

[results] = compute_outcomes(new_cal,1);

% 0.8458    0.3412    1.1063    0.5222    1.4602    0.4227
% 0.8458    0.3412    1.1063    0.5222    1.4602    0.4227

% options = optimoptions('fsolve', 'Display','iter','MaxFunEvals',1e3,'MaxIter',1e1,...
%     'Algorithm','trust-region-reflective','FinDiffRelStep',10^-3);
% 
% guess = cal;
% 
% tic
% new_cal = fsolve(@(xxx) calibrate_model(xxx,0),guess,options);
% toc