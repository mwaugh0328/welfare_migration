
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First add some paths as some code will be called from other folders, then
% add the calibration. Note the -r2 is in levels not logs and is called
% new_cal, so just comment out the line below.

addpath('../calibration')
addpath('../ge_taxation')

load calibration_final
new_cal = exp(new_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now let's run the main file, this creates the wages for the decentralized
% equillibrium

[targets, wage] = analyze_outcomes(new_cal, [], [], [], [], 1);

cd('..\effecient')

wage_de = [wage.monga, wage.notmonga];
%These are teh wages in the decentralized equillibrium, then we are going
%to pass this through one last time to get policy functions and the
%primitive TFP. The move policy function below is used to construct a good
%initial guess for the optimization and bounds on the problem.

[move_de, solve_types, assets, params, specs, vfun, ce] = just_policy(new_cal, wage_de, [], [], [], []);

[data_panel, params] = just_simmulate(params, move_de, solve_types, assets, specs, vfun, []);

[labor, govbc, tfp, ~, welfare_decentralized] = aggregate(params, data_panel, wage_de, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

best = load('move_fminsearch.mat','x1');

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('')
disp('Now Compute the Effecient Allocation...')

tic
[~, social_welfare] = compute_effecient(best.x1,  new_cal, tfp, 1);
toc

cons_eqiv.all = ((social_welfare.all ./ welfare_decentralized.all)).^(1./(1-params.pref_gamma)) - 1;
% This is just the standard thing. Think of guys behind the vale, so social
% welfare in the effecient allocation relative to decentralized. This is
% what each should recive (expost paths and outcomes may be different) but
% this is again a behind the vale calcuation.

% you could also compute, take this compared to a rural guy, what would he
% get in expectation if living in the effecient world or urban.
cons_eqiv.rural = ((social_welfare.all ./ welfare_decentralized.rural)).^(1./(1-params.pref_gamma)) - 1;
cons_eqiv.urban = ((social_welfare.all ./ welfare_decentralized.urban)).^(1./(1-params.pref_gamma)) - 1;

disp("Al, Welfare Gain in %: From Decentralized to Centralized/Effecient Allocaiton")
disp(100.*cons_eqiv.all)
disp("Rural, Welfare Gain in %: From Decentralized to Centralized/Effecient Allocaiton")
disp(100.*cons_eqiv.rural)
disp("Urban, Welfare Gain in %: From Decentralized to Centralized/Effecient Allocaiton")
disp(100.*cons_eqiv.urban)

rmpath('../calibration')
rmpath('../ge_taxation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All this stuff below is to solve for the effecient allocation, not just
% compute it given a solution which is what it does above.
% ntypes = 24;
% n_shocks = 5*2;
% 
% asset_loc = 15;
% 
% move_vec = [];
% count = 0;
% 
% A = [];
% rural_A = [eye(10), eye(10); zeros(10), zeros(10)];
% urban_A = zeros(10);
% 
% %14.0000   76.0000   12.0000    1.0000   
% %15.0000   79.0000   86.0000    2.0000
% %12.0000   92.0000   35.0000    3.0000 
% for xxx = 1:ntypes
%     
%     foo = squeeze(move_de(xxx).rural_not(12,:,:));
%   
%     foo = [foo(:,1),  diff(foo,1,2)];
%     foo = foo(:,1:2);
% 
%     move_vec = [move_vec ; foo(:)];
%     A = blkdiag(A,rural_A);
%     
%     foo = squeeze(move_de(xxx).rural_exp(92,:,:));
%     
%     foo = [foo(:,1),  diff(foo,1,2)];
%     foo = foo(:,1:2);
% 
%     move_vec = [move_vec ; foo(:)];
%     A = blkdiag(A,rural_A);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     foo = squeeze(move_de(xxx).urban_new(35,:,:));
% 
%     foo = [foo(:,1),  diff(foo,1,2)];
%     move_vec = [move_vec ; foo(:,1)];
%     A = blkdiag(A,urban_A);
%     
%     foo = squeeze(move_de(xxx).urban_old(3,:,:));
%     
%     foo = [foo(:,1),  diff(foo,1,2)];
%     move_vec = [move_vec ; foo(:,1)];
%     A = blkdiag(A,urban_A);
%     bar = 1;
%     
% end
% 
% testmove = make_movepolicy(move_vec,n_shocks,ntypes);
% 
% x0 = move_vec;
% 
% LB = (zeros(length(x0),1));
% UB = (ones(length(x0),1));
% b = ones(length(x0),1);

%tic
%[social_welfare] = compute_effecient(x0, new_cal, tfp, 1);
%toc

% opts = optimoptions(@fmincon,'Algorithm','active-set','Display','iter',...
% 'UseParallel',true,'FiniteDifferenceType','central',...
% 'FiniteDifferenceStepSize', 10^-3,'TolX',1e-3,'TolFun',10^-3,'MaxFunEvals',100000,...
% 'ConstraintTolerance',10^-3);
% 
% load('move_psearch_best5.mat');
% 
% [x1] = fmincon(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0), x1,A,b,[],[],(LB),(UB),[],opts);
% 

% two = load('move_psearch_best2.mat','x1');
% thr = load('move_psearch_best3.mat','x1');
% xinit = [x0'; one.x1';two.x1';thr.x1'];

% tic
% 
% opts = optimoptions('ga','Display','iter','UseParallel',true,'MaxGenerations',100,'TolFun',10^-3);
% 
% [x1] = ga(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0),length(x0), A,b,[],[],(LB),(UB),[],opts);
% 
% x1 = x1';
% 
% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% 
% tic
% 
% opts = optimoptions('patternsearch','Display','iter','UseParallel',true,'MaxFunEvals',10000,'TolFun',10^-3,'TolX',10^-3);
% 
% x1 = patternsearch(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0),x1, A,b,[],[],(LB),(UB),[],opts);
% 
% save move_psearch_best7 x1
% 
% toc
% tic 
% 
% opts = optimset('Display','iter','UseParallel',true,'MaxFunEvals',50000,'TolFun',10^-3,'TolX',10^-3);
% 
% x1 = fminsearchcon(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0), best.x1,LB,UB,A,b,[],opts);
% 
% save move_fminsearch_surf x1
% 
% toc

















