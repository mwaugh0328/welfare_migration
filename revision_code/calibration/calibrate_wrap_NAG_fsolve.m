function calibrate_wrap_NAG_fsolve
clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guess = [1.2876    0.5398    1.5543    0.7343    1.5153    0.6702    0.6309    0.5719    0.1066];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean(abs(calibrate_model(guess,2)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nval = length(guess);
diag_adjust = round(nval-1);

tic
[new_val, fvec, diag, nfev,~,~,~,~,ifail] = c05qc(@fcn, log(guess), int64(diag_adjust),...
    int64(diag_adjust), int64(1), ones(nval,1), int64(5),'epsfcn', 10^-5);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(new_val)


disp(new_val)
compute_outcomes_prefshock(exp(new_val),1);

save calibration_NAG new_val

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fvec, user, iflag] = fcn(n, x, fvec, user, iflag)
  if iflag ~=0
    fvec = zeros(n, 1);
    fvec(1:n) = calibrate_model(exp(x),2);
  else
    fprintf('objective = %10e\n', mean(abs(fvec)))
  end