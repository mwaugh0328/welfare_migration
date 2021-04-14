clear
warning('off','stats:regress:RankDefDesignMat');

addpath('../calibration')

load calibration_final
load wages

[move, solve_types, assets, params, specs, vfun, ce] = just_policy(exp(new_val), testwage, [], [], [], []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data_panel, params] = just_simmulate(params, move, solve_types, assets, specs, ce, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[labor, govbc, tfp] = aggregate(params, data_panel, testwage, [], 1);

cft = params.means_test;

taxprog = 0.0;

policyfun.move = move;
policyfun.assets = assets;

compute_eq([testwage; 1.0], exp(new_val), tfp, cft, vfun, taxprog, policyfun, 1)

options = optimoptions('fsolve', 'Display','iter','MaxFunEvals',2000,'MaxIter',20,...
'TolX',1e-3,'Algorithm','trust-region-reflective','FiniteDifferenceType','central',...
'FiniteDifferenceStepSize', 10^-3);

guess = [testwage; 1.0];

tic
[wageseq, ~, ~] = fsolve(@(xxx) compute_eq((xxx), exp(new_val), tfp, cft, [], taxprog, policyfun, 0), guess,options);
toc

disp(wageseq)

compute_eq([wageseq], exp(new_val), tfp, cft, vfun, taxprog, policyfun, 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guess = [testwage; 1.0];

tic
[wageseq, ~, ~] = fsolve(@(xxx) compute_eq((xxx), exp(new_val), tfp, cft, [], taxprog,[], 0), guess,options);
toc

disp(wageseq)

compute_eq([wageseq], exp(new_val), tfp, cft, vfun, taxprog, [], 1)