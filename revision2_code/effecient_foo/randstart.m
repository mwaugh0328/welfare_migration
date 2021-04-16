addpath('../calibration')

load calibration_final
load wages


[move_de, solve_types, assets, params, vfun, ce] = just_policy(exp(new_val), wages, [], [], [], []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data_panel, params] = just_simmulate(params, move_de, solve_types, assets, vfun, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[labor, govbc, tfp] = just_aggregate(params,data_panel, wages, [], 0);

recval = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntypes = 24;
n_shocks = 5*2;

nruns = 1000;
rng(02071983)
asset_loc = randi(100,nruns,4);

for zzz = 1:nruns

    
disp(zzz)

move_vec = [];
count = 0;

A = [];
rural_A = [eye(10), eye(10); zeros(10), zeros(10)];
urban_A = zeros(10);

for xxx = 1:ntypes
    
    foo = squeeze(move_de(xxx).rural_not(asset_loc(zzz,1),:,:));
  
    foo = [foo(:,1),  diff(foo,1,2)];
    foo = foo(:,1:2);

    move_vec = [move_vec ; foo(:)];
    A = blkdiag(A,rural_A);
    
    foo = squeeze(move_de(xxx).rural_exp(asset_loc(zzz,2),:,:));
    
    foo = [foo(:,1),  diff(foo,1,2)];
    foo = foo(:,1:2);

    move_vec = [move_vec ; foo(:)];
    A = blkdiag(A,rural_A);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    foo = squeeze(move_de(xxx).urban_new(asset_loc(zzz,3),:,:));

    foo = [foo(:,1),  diff(foo,1,2)];
    move_vec = [move_vec ; foo(:,1)];
    A = blkdiag(A,urban_A);
    
    foo = squeeze(move_de(xxx).urban_old(asset_loc(zzz,4),:,:));
    
    foo = [foo(:,1),  diff(foo,1,2)];
    move_vec = [move_vec ; foo(:,1)];
    A = blkdiag(A,urban_A);
    bar = 1;
    
end

x0 = move_vec;

tic
[social_welfare] = compute_effecient(x0, exp(new_val), tfp, 1);
toc

recval = [recval; asset_loc(zzz,:), social_welfare];

end

save diffstart recval