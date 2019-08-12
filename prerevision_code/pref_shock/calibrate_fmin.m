clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');

options = optimset('Display','iter','MaxFunEvals',500,'MaxIter',1e6,'TolFun',1e-3,'TolX',1e-10);

guess = [1.2959    0.5367    1.5555    0.7254    1.5138    0.6953    0.6437    0.5856    0.1054];


new_val = fminsearch(@(xxx)calibrate_model(exp(xxx),1),log(guess),options);

disp(new_val)
compute_outcomes_prefshock((new_val),1);
