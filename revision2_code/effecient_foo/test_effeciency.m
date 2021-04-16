

addpath('../calibration')

load calibration_final
load wages


% [move_de, solve_types, assets, params, vfun, ce] = just_policy(exp(new_val), wages, [], [], [], []);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [data_panel, params] = just_simmulate(params, move_de, solve_types, assets, vfun, []);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [labor, govbc, tfp] = just_aggregate(params,data_panel, wages, [], 0);


[move, solve_types, assets, params, specs, vfun, ce] = just_policy(exp(new_val), testwage, [], [], [], []);

[data_panel, params] = just_simmulate(params, move, solve_types, assets, specs, ce, []);

[labor, govbc, tfp] = aggregate(params, data_panel, testwage, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntypes = 24;
n_shocks = 5*2;

asset_loc = 15;

move_vec = [];
count = 0;

A = [];
rural_A = [eye(10), eye(10); zeros(10), zeros(10)];
urban_A = zeros(10);

for xxx = 1:ntypes
    
    foo = squeeze(move_de(xxx).rural_not(15,:,:));
  
    foo = [foo(:,1),  diff(foo,1,2)];
    foo = foo(:,1:2);

    move_vec = [move_vec ; foo(:)];
    A = blkdiag(A,rural_A);
    
    foo = squeeze(move_de(xxx).rural_exp(79,:,:));
    
    foo = [foo(:,1),  diff(foo,1,2)];
    foo = foo(:,1:2);

    move_vec = [move_vec ; foo(:)];
    A = blkdiag(A,rural_A);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    foo = squeeze(move_de(xxx).urban_new(86,:,:));

    foo = [foo(:,1),  diff(foo,1,2)];
    move_vec = [move_vec ; foo(:,1)];
    A = blkdiag(A,urban_A);
    
    foo = squeeze(move_de(xxx).urban_old(2,:,:));
    
    foo = [foo(:,1),  diff(foo,1,2)];
    move_vec = [move_vec ; foo(:,1)];
    A = blkdiag(A,urban_A);
    bar = 1;
    
end



move_vec = make_movevec(move_de, 12, params)

testmove = make_movepolicy(move_vec,n_shocks,ntypes);


x0 = move_vec;

tic
[social_welfare] = compute_effecient(x0, exp(new_val), tfp, 1);
toc

LB = (zeros(length(x0),1));
UB = (ones(length(x0),1));
b = ones(length(x0),1);

% opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
% 'UseParallel',true,'FiniteDifferenceType','central',...
% 'FiniteDifferenceStepSize', 10^-3,'TolX',1e-2,'MaxFunEvals',50000,...
% 'ConstraintTolerance',10^-3);
% 
% 
% [x1] = fmincon(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0), x1,A,b,[],[],(LB),(UB),[],opts);

% opts = optimoptions('ga','Display','iter','UseParallel',true,'MaxGenerations',10);
% 
% [x1] = ga(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0),length(x0), A,b,[],[],(LB),(UB),[],opts);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

opts = optimoptions('patternsearch','Display','iter','UseParallel',true,'MaxFunEvals',100000);

x1 = patternsearch(@(xxx) compute_effecient(xxx, exp(new_val), tfp, 0), x0, A,b,[],[],(LB),(UB),[],opts);

toc

% 
% 
% % tic
% % [x1,fval] = runobjconstr(x1,  exp(new_val), tfp, opts);
% % toc
% 
tic
[social_welfare] = compute_effecient(x1,  exp(new_val), tfp, 1);
toc


















