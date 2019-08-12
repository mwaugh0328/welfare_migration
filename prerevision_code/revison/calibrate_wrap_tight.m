clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.3065    0.4961    1.4559    0.7239    1.4440    0.5445    0.4921    0.6590
%1.3066    0.4961    1.4598    0.7219    1.4479    0.5745    0.5109    0.6590
%guess = [1.6523    0.4652    1.3913    0.8073    1.5226    0.5389    0.5926    0.6049];

guess = [1.2497    0.5116    1.5024    0.7129    1.4788    0.6553    0.6049    0.6587]

upper_bound = [1.75, 0.65, 1.70, 0.85, 1.75, 0.80, 0.85, 0.85]; 
lower_bound = [1.10, 0.40, 1.20, 0.60, 1.4, 0.35, 0.35, 0.35]; 

ObjectiveFunction = @(x) calibrate_model((x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options = gaoptimset('Display','iter','Generations',10,'Initialpopulation',guess);
% % 
% cal = ga(ObjectiveFunction, length(upper_bound),[],[],[],[],(lower_bound),(upper_bound),[],options)
% 
% save calibration_ga cal
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ObjectiveFunction = @(x) calibrate_model((x),1);

options_pa = optimoptions('patternsearch','Display','iter','MaxFunEvals',5e3);

new_cal = patternsearch(ObjectiveFunction,(guess),[],[],[],[],(lower_bound),(upper_bound),[],options_pa) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save calibration_final new_cal

[results] = compute_outcomes(new_cal,1);


options = optimset('Display','iter','MaxFunEvals',200,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-5);

guess = [1.2261    0.5213    1.5200    0.7139    1.4759    0.6443    0.6064    0.6292];
%guess = log([4.7675    0.2327    1.3456    0.9995]);

model_params = fminsearch(@(x) calibrate_model((x),1),guess,options);

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