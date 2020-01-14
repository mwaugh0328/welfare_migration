clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');

options = optimset('Display','iter','MaxFunEvals',500,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-10);

guess = [1.1903    0.4466    1.4888    0.7620   -0.3595    0.5266    0.0852    1.1220    0.1088];

[new_val,fval]= fminsearch(@(xxx)calibrate_model((xxx),1),(guess) + .01*randn(1,length(guess)),options);
 
disp(new_val)

save calibration_addit new_val

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.8968 wage gap
%0.6741 rural share
%0.4317 variance 
%0.5181 assets
%0.4661 migration rate
%0.3004 jump in migration
%-0.0071 next period migration
%0.1478 LATE
%0.0533 OLS
%0.6020 repeat migration
%0.1844 variance oncsutp
% 1.8951    0.6741    0.4317    0.5181    0.4661    0.3004   -0.0071    0.1478    0.0533    0.6020    0.1844