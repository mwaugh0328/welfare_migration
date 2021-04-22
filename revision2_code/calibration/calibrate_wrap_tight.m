clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.3065    0.4961    1.4559    0.7239    1.4440    0.5445    0.4921    0.6590
%1.3066    0.4961    1.4598    0.7219    1.4479    0.5745    0.5109    0.6590
% guess = [1.2893    0.5511    1.5691    0.7426    1.5548    0.5939    0.6258    0.5327    0.1161 0.08];
% 
% guess = [1.2879    0.5484    1.5817    0.7409    1.5645    0.5909    0.5919    0.5339    0.1172    0.0969] + .01*randn(1,length(guess));
% 
% %guess = [1.0333    0.5488    1.5184    0.6904    0.8454    0.4534    0.6359    0.7371]
% 
% upper_bound = [1.75, 0.60, 1.70, 0.95, 1.9, 0.85, 0.85, 0.95, 0.30]; 
% lower_bound = [1.00, 0.40, 1.20, 0.25, 1.1, 0.20, 0.35, 0.15, 0.01]; 

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

% ObjectiveFunction = @(x) calibrate_model(exp(x),1);
% 
% options_pa = optimoptions('patternsearch','Display','iter','MaxFunEvals',2e3);
% 
% new_cal = patternsearch(ObjectiveFunction,log(guess),[],[],[],[],[],[],[],options_pa) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('Display','iter','MaxFunEvals',100);

load('calibration_final.mat')

tic
[new_cal, fval] =fminsearch(@(xxx) calibrate_model(exp(xxx),[],1),new_val,options);
toc


