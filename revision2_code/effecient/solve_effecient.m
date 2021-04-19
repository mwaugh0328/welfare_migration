

addpath('../calibration')

load calibration_final
load wages


[move_de, solve_types, assets, params, vfun, ce] = just_policy(exp(new_val), wages, [], [], [], []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data_panel, params] = just_simmulate(params, move_de, solve_types, assets, vfun, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[labor, govbc, tfp] = just_aggregate(params,data_panel, wages, [], 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntypes = 24;
n_shocks = 5*2;

asset_loc = 15;

move_vec = [];
count = 0;

A = [];
rural_A = [eye(10), eye(10); zeros(10), zeros(10)];
urban_A = zeros(10);

%14.0000   76.0000   12.0000    1.0000   
%15.0000   79.0000   86.0000    2.0000
%12.0000   92.0000   35.0000    3.0000 
for xxx = 1:ntypes
    
    foo = squeeze(move_de(xxx).rural_not(12,:,:));
  
    foo = [foo(:,1),  diff(foo,1,2)];
    foo = foo(:,1:2);

    move_vec = [move_vec ; foo(:)];
    A = blkdiag(A,rural_A);
    
    foo = squeeze(move_de(xxx).rural_exp(92,:,:));
    
    foo = [foo(:,1),  diff(foo,1,2)];
    foo = foo(:,1:2);

    move_vec = [move_vec ; foo(:)];
    A = blkdiag(A,rural_A);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    foo = squeeze(move_de(xxx).urban_new(35,:,:));

    foo = [foo(:,1),  diff(foo,1,2)];
    move_vec = [move_vec ; foo(:,1)];
    A = blkdiag(A,urban_A);
    
    foo = squeeze(move_de(xxx).urban_old(3,:,:));
    
    foo = [foo(:,1),  diff(foo,1,2)];
    move_vec = [move_vec ; foo(:,1)];
    A = blkdiag(A,urban_A);
    bar = 1;
    
end

testmove = make_movepolicy(move_vec,n_shocks,ntypes);


x0 = move_vec;

tic
[social_welfare] = compute_effecient(x0, exp(new_val), tfp, 1);
toc

LB = (zeros(length(x0),1));
UB = (ones(length(x0),1));
b = ones(length(x0),1);


% opts = optimoptions(@fmincon,'Algorithm','active-set','Display','iter',...
% 'UseParallel',true,'FiniteDifferenceType','central',...
% 'FiniteDifferenceStepSize', 10^-3,'TolX',1e-3,'TolFun',10^-3,'MaxFunEvals',100000,...
% 'ConstraintTolerance',10^-3);
% 
% load('move_psearch_best5.mat');
% 
% [x1] = fmincon(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0), x1,A,b,[],[],(LB),(UB),[],opts);
% 
% two = load('move_psearch_best1.mat','x1');
% two = load('move_psearch_best2.mat','x1');
% thr = load('move_psearch_best3.mat','x1');
% xinit = [x0'; one.x1';two.x1';thr.x1'];

tic

opts = optimoptions('ga','Display','iter','UseParallel',true,'MaxGenerations',100,'TolFun',10^-3);

[x1] = ga(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0),length(x0), A,b,[],[],(LB),(UB),[],opts);

x1 = x1';

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

tic

opts = optimoptions('patternsearch','Display','iter','UseParallel',true,'MaxFunEvals',10000,'TolFun',10^-3,'TolX',10^-3);

x1 = patternsearch(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0),x1, A,b,[],[],(LB),(UB),[],opts);

save move_psearch_best7 x1

toc

% two = load('move_psearch_best2.mat','x1');
% 
% opts = optimset('Display','iter','UseParallel',true,'MaxFunEvals',2000,'TolFun',10^-3,'TolX',10^-3);
% 
% x1 = fminsearchcon(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0),two.x1,LB,UB,A,b,[],opts);
% 
% save move_fminsearch_best7 x1


tic
[social_welfare] = compute_effecient(x1,  exp(new_val), tfp, 1);
toc


















