clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guess = [1.1169    1.7259    0.7472    1.7128    0.2679    0.7231    1.5];

upper_bound = [2.00, 2.00, 0.85, 2.00, 0.75, 0.75, 2.50]; 
lower_bound = [0.01, 1.00, 0.10, 1.00, 0.25, 0.25, 0.35]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ObjectiveFunction = @(x) calibrate_model((x),1);

options = gaoptimset('Display','iter','Generations',20,'Initialpopulation',guess);

cal = ga(ObjectiveFunction, length(upper_bound),[],[],[],[],(lower_bound),(upper_bound),[],options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options_pa = optimoptions('patternsearch','Display','iter','MaxFunEvals',1e3);

guess = cal;

new_cal = patternsearch(ObjectiveFunction,guess,[],[],[],[],(lower_bound),(upper_bound),[],options_pa) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save calibration_noselection new_cal

[results] = compute_outcomes(new_cal,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different permutations on this...
% The main calibration...

% guess = [1.3000    0.4946    1.4758    0.7377    1.4382    0.4478    0.5322    0.6907];
% 
% %upper_bound = [2.00, 0.75, 2.00, 0.85, 2.00, 0.75, 0.75,2.00]; 
% %lower_bound = [0.01, 0.10, 1.00, 0.10, 1.00, 0.25, 0.25,0.35]; 
% 
% upper_bound = [1.50, 0.55, 1.60, 0.85, 1.70, 0.75, 0.75,1.00]; 
% lower_bound = [0.90, 0.45, 1.30, 0.40, 1.30, 0.25, 0.25,0.45]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No Selection...

guess = [1.1169    1.7259    0.7472    1.7128    0.2679    0.7231    1.9998];

upper_bound = [2.00, 2.00, 0.85, 2.00, 0.75, 0.75, 2.50]; 
lower_bound = [0.01, 1.00, 0.10, 1.00, 0.25, 0.25, 0.35]; 