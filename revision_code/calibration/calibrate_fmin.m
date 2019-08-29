clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');

options = optimset('Display','iter','MaxFunEvals',50,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-10);

guess = [1.2906    0.5508    1.5654    0.7403    1.5519    0.5964    0.6312    0.5375    0.1198];

[new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),log(guess),options);

disp(new_val)

%save calibration_highgrid new_val fval

[new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),new_val,options);

%save calibration_highgrid new_val fval

%[new_val,fval]= fminsearch(@(xxx)calibrate_model(exp(xxx),1),new_val,options);

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

%save calibration_highgrid new_val fval
