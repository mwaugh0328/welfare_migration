clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%guess = [1.2893    0.5511    1.5691    0.7426    1.5548    0.5939    0.6258    0.5327    0.1161];
%guess(4) = []; % for ubar to be removed

guess = [1.2624    0.5637    1.5557    1.5436    0.4389    0.7698    0.8639    0.1674]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 addpath('C:\Users\mwaugh\github\welfare_migration\revision_code\calibration')
% 
 ObjectiveFunction = @(x) calibrate_model_norho(exp(x),1);
 
%upper_bound = [2.00, 2.00, 0.85, 2.00, 0.75, 2.00]; 
%lower_bound = [0.01, 1.00, 0.10, 1.00, 0.25, 0.35];

max_funs = 20;

options_pa = optimoptions('patternsearch','Display','iter','MaxFunEvals',max_funs);

new_val = patternsearch(ObjectiveFunction,log(guess) + 0.01*randn(1,length(guess)),[],[],[],[],[],[],[],options_pa) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options_fmin = optimset('Display','iter','MaxFunEvals',500,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-10);

[new_val, fval]= fminsearch(ObjectiveFunction,new_val,options_fmin);

rmpath('C:\Users\mwaugh\github\welfare_migration\revision_code\calibration')

save calibration_norho new_val fval


