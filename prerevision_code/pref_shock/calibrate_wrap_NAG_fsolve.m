function calibrate_wrap_NAG_fsolve
clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2220    0.4842    1.4478    0.7088    1.4486    0.6276    0.5188    0.6638

guess = [1.3710    0.5350    1.5770    0.7220    1.5416    0.7667    0.3501    0.6190];

must_be_positive = [1,3,5,8];
must_be_zero_one = [2,4,6,7];

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
mean(abs(calibrate_model_NAG(guess,1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nval = length(guess);
diag_adjust = round(nval-1);

tic
[new_val, fvec, diag, nfev,~,~,~,~,ifail] = c05qc(@fcn, (guess), int64(diag_adjust),...
    int64(diag_adjust), int64(1), ones(nval,1), int64(5),'epsfcn', 10^-4);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(new_val)
params = zeros(8,1);
params(must_be_positive) = (new_val(must_be_positive));
params(must_be_zero_one) = (new_val(must_be_zero_one));

disp(params)
compute_outcomes_prefshock(params,1);

save calibration_test_NAG params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fvec, user, iflag] = fcn(n, x, fvec, user, iflag)
  if iflag ~=0
    fvec = zeros(n, 1);
    fvec(1:n) = calibrate_model_NAG((x),1);
  else
    fprintf('objective = %10e\n', mean(abs(fvec)))
  end