function calibrate_wrap_NAG_fsolve
clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2170    0.5158    1.5150    0.7087    1.4996    0.7665    0.6968    0.6146
% high 48/30 grid...

%guess = [1.1975    0.5158    1.5150    0.7009    1.4996    0.7665    0.6791    0.6322];
%1.2438    0.5276    1.5266    0.7028    1.4990    0.7899    0.6936    0.5933
guess = [1.2758    0.5493    1.6050    0.7237    1.5138    0.6106    0.6157    0.6322    0.1573];


% guess = zeros(1,length(old_cal)+2);
% guess(must_be_positive) = log(old_cal(must_be_positive));
% guess(must_be_zero_one) = -log(1./old_cal(must_be_zero_one)-1);

% guess = zeros(1,length(old_cal)+2);
% guess(must_be_positive) = (old_cal(must_be_positive));
% guess(must_be_zero_one) = old_cal(must_be_zero_one);
% 
% guess(end-1) = log(2.*10^-1);
% guess(end) = log(2.*10^-1);

%guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean(abs(calibrate_model(guess,2)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nval = length(guess);
diag_adjust = round(nval-1);

tic
[new_val, fvec, diag, nfev,~,~,~,~,ifail] = c05qc(@fcn, (guess), int64(diag_adjust),...
    int64(diag_adjust), int64(1), ones(nval,1), int64(5),'epsfcn', 10^-5);
toc


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options = optimset('Display','iter','MaxFunEvals',100,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-5);
% 
% 
% %guess = log([4.7675    0.2327    1.3456    0.9995]);
% 
% model_params = fminsearch(@(xxx)calibrate_model_NAG((xxx),1),new_val,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(new_val)


disp(new_val)
compute_outcomes_prefshock(exp(new_val),1);

save calibration_test_NAG params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fvec, user, iflag] = fcn(n, x, fvec, user, iflag)
  if iflag ~=0
    fvec = zeros(n, 1);
    fvec(1:n) = calibrate_model((x),2);
  else
    fprintf('objective = %10e\n', mean(abs(fvec)))
  end