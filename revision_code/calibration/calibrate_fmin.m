clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');

options = optimset('Display','iter','MaxFunEvals',300,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-10);

guess = [1.2840    0.5343    1.5437    0.7462    1.5024    0.6541    0.6021    0.5948    0.1069];
%1.2881    0.5357    1.5406    0.7477    1.5028    0.6481    0.6038    0.5908    0.1084
% [new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),log(guess),options);
% 
% disp(new_val)
% % 
%  save calibration_best new_val fval
% 
% [new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),new_val,options);
% 
%  save calibration_best new_val fval

%[new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),new_val,options);

% compute_outcomes_prefshock(exp(new_val),1);

%rng(08182016)
%rng(08182019)
% rng(08182017)
% 
niters = 10;
new_val = zeros(niters,length(guess));
fval = zeros(niters,1);

for zzz = 1:niters
rng(0818201 + zzz)

[new_val(zzz,:),fval(zzz)]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),log(guess) + .01*randn(1,length(guess)),options);

end
% 
 disp(new_val)

 save calibration_multistart new_val fval
