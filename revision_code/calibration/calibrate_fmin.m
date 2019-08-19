clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');

options = optimset('Display','iter','MaxFunEvals',500,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-10);

guess = [1.2871    0.5397    1.5547    0.7345    1.5151    0.6703    0.6310    0.5718    0.1063];

[new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),log(guess),options);

save calibration_highgrid new_val fval

[new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),new_val,options);

save calibration_highgrid new_val fval

[new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),new_val,options);

% rng(08182016)
% 
% niters = 10;
% new_val = zeros(niters,length(guess));
% fval = zeros(niters,1);
% 
% for zzz = 1:niters
% rng(08182016 + zzz)
% 
% [new_val(zzz,:),fval(zzz)]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),log(guess) + .01*randn(1,length(guess)),options);
% 
% end

disp(new_val)
%compute_outcomes_prefshock(exp(new_val),1);

save calibration_highgrid new_val fval
