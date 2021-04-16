function [move] = make_migration(params,choices)

nstates = params.n_shocks;

rural_options = params.rural_options;
urban_options = params.urban_options;

n_rural_choice = 2*nstates*(rural_options-1);

rural_choice = choices(1:n_rural_choice);
urban_choice = choices(n_rural_choice+1:end);

rural_choice = reshape(rural_choice,nstates,rural_options-1,2);
urban_choice = reshape(urban_choice,nstates,urban_options-1,2);

move.rural_not = rural_choice(:,:,1);
move.rural_not = [move.rural_not, ones(nstates,1)];

move.rural_exp = rural_choice(:,:,2);
move.rural_exp = [move.rural_exp, ones(nstates,1)];

move.urban_new = urban_choice(:,:,1);
move.urban_new = [move.urban_new, ones(nstates,1)];

move.urban_old = urban_choice(:,:,2);
move.urban_old = [move.urban_old, ones(nstates,1)];


