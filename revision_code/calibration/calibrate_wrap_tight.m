clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.3065    0.4961    1.4559    0.7239    1.4440    0.5445    0.4921    0.6590
%1.3066    0.4961    1.4598    0.7219    1.4479    0.5745    0.5109    0.6590
guess = [1.2779    0.5376    1.5492    0.7439    1.5131    0.6706    0.6244    0.5661    0.1040];

%guess = [1.0333    0.5488    1.5184    0.6904    0.8454    0.4534    0.6359    0.7371]

upper_bound = [1.75, 0.60, 1.70, 0.90, 1.7, 0.85, 0.85, 0.85, 0.20]; 
lower_bound = [1.00, 0.40, 1.20, 0.50, 1.4, 0.20, 0.35, 0.35, 0.04]; 

% ObjectiveFunction = @(x) calibrate_model((exp(x)),1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% options = gaoptimset('Display','iter','Generations',10,'Initialpopulation',log(guess));
% % 
% cal = ga(ObjectiveFunction, length(upper_bound),[],[],[],[],log(lower_bound),log(upper_bound),[],options)
% 
% save calibration_ga cal
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ObjectiveFunction = @(x) calibrate_model(exp(x),1);

options_pa = optimoptions('patternsearch','Display','iter','MaxFunEvals',5e2);


new_cal = patternsearch(ObjectiveFunction,log(guess),[],[],[],[],[],[],[],options_pa) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save calibration_final new_cal

[results] = compute_outcomes_prefshock(exp(new_cal),1);

% Just Control...
% 1.1075    0.5033    1.3988    0.6742    1.5932    0.6849    0.3735
% Control and Experiment
% 1.2639    0.5179    1.4820    0.7380    1.4825    0.4578    0.5423
%
% Hybrid
% 1.0703    0.4909    1.4755    0.6660    1.4657    0.4700    0.5964
% 0.9870    0.4961    1.4717    0.6283    1.5053    0.5114    0.4879
% 1.0734    0.5010    1.5054    0.6682    1.4772    0.4723    0.5582

% 

% options = optimoptions('fsolve', 'Display','iter','MaxFunEvals',1e3,'MaxIter',1e1,...
%     'Algorithm','trust-region-reflective','FinDiffRelStep',10^-3);
% 
% guess = cal;
% 
% tic
% new_cal = fsolve(@(xxx) calibrate_model(xxx,0),guess,options);
% toc