clc; clear;
close all

warning('off','stats:regress:RankDefDesignMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.3065    0.4961    1.4559    0.7239    1.4440    0.5445    0.4921    0.6590
%1.3066    0.4961    1.4598    0.7219    1.4479    0.5745    0.5109    0.6590
guess = [1.2893    0.5511    1.5691    0.7426    1.5548    0.5939    0.6258    0.5327    0.1161 0.08];

guess = [1.2879    0.5484    1.5817    0.7409    1.5645    0.5909    0.5919    0.5339    0.1172    0.0969] + .01*randn(1,length(guess));
rns =50;

new_cal = zeros(rns,10);
fval(rns,1);

for xxx = 1:rns
 
    tic
    [new_cal(xxx,:), fval(xxx)] =fminsearch(@(xxx) calibrate_model(exp(xxx),1),log(guess),options);
    toc
    
    guess = exp(new_cal(xxx,:))+ .01*randn(1,length(guess));
    
end