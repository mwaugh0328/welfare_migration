function [policy] = mobility_labor_supply_smth(params,wages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the households problem
n_iterations = 2000;

tol = 10^-3; % Not sensitive to lower values

n_options = params.n_options;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_smooth = params.bsmooth;
w_smooth = params.wsmooth;

R = params.R;

beta = params.beta;

B = params.labor_disutility;

m = params.m;

trans_mat = params.trans_mat;
invar_trans_mat = params.invar_trans_mat';

%grid = params.grid;

home_production = params.home_production;
% convert into units of final good. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks = wages;
n_shocks = length(wages);

%A = (1-gamma).^-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up grid for asset holdings.

asset_space = params.asset_space;
n_asset_states = length(asset_space);

% asset_grid(:,i) == asset_space(i)
asset_grid  = meshgrid(asset_space,asset_space);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the matricies for value function itteration

v_prime = zeros(n_asset_states,n_shocks);

policy_assets = zeros(n_asset_states,n_shocks,n_options);
policy_move = zeros(n_asset_states,n_shocks,n_options);
value_fun_out = zeros(n_asset_states,n_shocks,n_options);
% v_old = zeros(n_asset_states,n_shocks); % this is the standard guess for
% asset states and for the shocks.

% we have net_assets(i, j) = R * a_i  - a_j
net_assets = R.*asset_grid' - asset_grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I'm going to pre-generate the period utility function and then the
% feasibility conditions... layout is (a, a', w)

utility_stay = zeros(n_asset_states,n_asset_states,n_shocks);
utility_move = zeros(n_asset_states,n_asset_states,n_shocks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
utility_rest = log(max(net_assets + home_production,10^-8)) ;

feasible_rest = net_assets > 0;

utility_rest = utility_rest + -1e10.*(~feasible_rest);

utility_rest_move = log(max(net_assets + home_production - m,10^-8)) ;

feasible_rest_move = net_assets - m > 0;

utility_rest_move = utility_rest_move + -1e10.*(~feasible_rest_move);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% feasible_stay = false(n_asset_states,n_asset_states,n_shocks);
% feasible_move = false(n_asset_states,n_asset_states,n_shocks);

for zzz = 1:n_shocks
    utility_stay(:,:,zzz) = log(max(net_assets + shocks(zzz), 10^-8)) - B;
    utility_move(:,:,zzz) = log(max(net_assets + shocks(zzz)- m, 10^-8)) - B;

    feasible_stay = net_assets + shocks(zzz) > 0;
    feasible_move = net_assets + shocks(zzz) - m > 0;

    utility_stay(:,:,zzz) = utility_stay(:,:,zzz) + -1e10.*(~feasible_stay);
    utility_move(:,:,zzz) = utility_move(:,:,zzz) + -1e10.*(~feasible_move);
end

v_old = repmat(diag(median(utility_stay,3)),1,n_shocks)./(1-beta);
v_old = ones(size(v_old));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commence value function itteration.

for iter = 1:n_iterations

    v_hat = v_old;

    value_newjob = beta.*(invar_trans_mat*v_hat');
    % New job is expected value across all posible locations....

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute value of Resting/not moving...note that this is independent
    % of the shock, so can be done outside loop...

    value_fun = bsxfun(@plus,utility_rest_move,value_newjob);

    [v_rest_move, ~ ] = max(value_fun,[],2);
    
    value_fun_R = value_fun - utility_rest_move;
    % Now work through different options that do depend upon the shock...    
    for zzz = 1:n_shocks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute value of styaing...

        expected_value = beta.*(trans_mat(zzz,:)*v_hat');

        value_fun = bsxfun(@plus,utility_stay(:,:,zzz),expected_value);

        [v_stay, ~ ] = max(value_fun,[],2);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute value of Resting/not moving...note the ording here is to
    % exploit that value_fun in memoroy from above...

        value_fun = value_fun + utility_rest - utility_stay(:,:,zzz);
        
        %value_fun_rtest = bsxfun(@plus,utility_rest,expected_value);
        
        % here take the value fun above, swap out the period utility
        % function since the epected value part is the same. This avoides
        % doing the bsx function again...some loss in clarity here...

        [v_rest, ~ ] = max(value_fun,[],2);
        
        %norm(value_fun-value_fun_rtest)
        
    % So if I rest I get the expected value of the next period...

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute value of working/moving
        value_fun = value_fun_R + utility_move(:,:,zzz);
        
        %value_fun_mtest = bsxfun(@plus,utility_move(:,:,zzz),value_newjob);
        
        %norm(value_fun-value_fun_mtest)

        [v_move, ~ ] = max(value_fun,[],2);

    % This says the value of moving is equall to
    % beta times the new job. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pi_denom_work_stay = exp(v_stay./w_smooth) + exp(v_rest./w_smooth) ;
    pi_denom_work_move = exp(v_move./w_smooth) + exp(v_rest_move./w_smooth) ;
    

    E_utility_stay = w_smooth.*log(pi_denom_work_stay);
    E_utility_move = w_smooth.*log(pi_denom_work_move);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the move or not nest
    
    pi_stay = exp(E_utility_stay./b_smooth);
    
    pi_move = exp(E_utility_move./b_smooth);
    
    pi_denom_move_stay = pi_stay + pi_move;
            
    %v_prime(:,zzz) = pi_stay.*v_stay + pi_move.*v_move +  pi_rest.*v_rest + pi_rest_move.*v_rest_move;
    v_prime(:,zzz) = b_smooth.*log(pi_denom_move_stay);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now do the max across these different options...
%     pi_stay = exp(v_stay./b_smooth);
%     
%     pi_move = exp(v_move./b_smooth);
%     
%     pi_rest = exp(v_rest./b_smooth);
%     
%     pi_rest_move = exp(v_rest_move./b_smooth);
%     
%     pi_denom = pi_stay + pi_move + pi_rest + pi_rest_move;
%         
%    %v_prime(:,zzz) = pi_stay.*v_stay + pi_move.*v_move +  pi_rest.*v_rest + pi_rest_move.*v_rest_move;
%    %
%    
%    %v_prime(:,zzz) = (pi_stay.*v_stay + pi_move.*v_move +  pi_rest.*v_rest + pi_rest_move.*v_rest_move)./pi_denom;
%    v_prime(:,zzz) = b_smooth.*log(pi_denom);
%    
   % norm([v_prime(:,zzz)- max([ v_stay , v_move, v_rest, v_rest_move],[],2)])
    %[v_prime(:,zzz), ~ ] = max([ v_stay , v_move, v_rest, v_rest_move],[],2) ;

    v_hat(:,zzz) = v_prime(:,zzz);
    % update the value function within the itteration
    end
    
    % disp(iter)
    % disp([norm(v_old - v_prime,Inf),iter])
    if norm(v_old - v_prime,Inf) < tol
        break
    end
    
%     if isnan(norm(v_old - v_prime,Inf))
%         break
%     end
    
    v_old = v_prime;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now pull out the policy function


value_newjob = beta.*(invar_trans_mat*v_hat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute value of Resting/not moving

value_fun = bsxfun(@plus,utility_rest_move,value_newjob);

[v_rest_move, p_asset_rest_move] = max(value_fun,[],2);

value_fun_R = value_fun - utility_rest_move;

for zzz = 1:n_shocks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute value of styaing...

        expected_value = beta.*(trans_mat(zzz,:)*v_hat');

        value_fun = bsxfun(@plus,utility_stay(:,:,zzz),expected_value);

        [v_stay, p_asset_stay] = max(value_fun,[],2);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute value of Resting/not moving
        
        value_fun = value_fun + utility_rest - utility_stay(:,:,zzz);
        
        %value_fun = bsxfun(@plus,utility_rest,expected_value);
        
        % here take the value fun above, swap out the period utility
        % function since the epected value part is the same. This avoides
        % doing the bsx function again...some loss in clarity here...

        [v_rest, p_asset_rest] = max(value_fun,[],2);        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute value of working/moving
    
    value_fun = value_fun_R + utility_move(:,:,zzz);
        
        %value_fun = bsxfun(@plus,utility_move(:,:,zzz),value_newjob);
        
        %norm(value_fun-value_fun_mtest)

    [v_move, p_asset_move] = max(value_fun,[],2);

    % So this is the new part. This says the value of moving is equall to
    % beta times the value associated with the best job. Given the
    % organization of the state space, the best job is at the end of the
    % coloumn of the matrix.

    policy_assets(:,zzz,:)= [ p_asset_stay , p_asset_move, p_asset_rest, p_asset_rest_move ];
        
%     pi_denom = exp(v_stay./b_smooth) + exp(v_move./b_smooth) + exp(v_rest./b_smooth) + exp(v_rest_move./b_smooth);
%     
%     pi_stay = exp(v_stay./b_smooth)./pi_denom;
%     
%     pi_move = exp(v_move./b_smooth)./pi_denom;
%     
%     pi_rest = exp(v_rest./b_smooth)./pi_denom;
%     
%     pi_rest_move = exp(v_rest_move./b_smooth)./pi_denom;
%     
%     policy_move(:,zzz,:) = [pi_stay, pi_move, pi_rest, pi_rest_move];
%     
%     %v_prime(:,zzz) = pi_stay.*v_stay + pi_move.*v_move +  pi_rest.*v_rest + pi_rest_move.*v_rest_move;
%     v_prime(:,zzz) = b_smooth.*log(pi_denom);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the work or not nest
    
    pi_work_stay = exp(v_stay./w_smooth);
    pi_work_move = exp(v_move./w_smooth);
    
    pi_denom_work_stay = exp(v_stay./w_smooth) + exp(v_rest./w_smooth) ;
    pi_denom_work_move = exp(v_move./w_smooth) + exp(v_rest_move./w_smooth) ;
    

    E_utility_stay = w_smooth.*log(pi_denom_work_stay);
    E_utility_move = w_smooth.*log(pi_denom_work_move);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the move or not nest
    
    pi_stay = exp(E_utility_stay./b_smooth);
    
    pi_move = exp(E_utility_move./b_smooth);
    
    pi_denom_move_stay = pi_stay + pi_move;
    
    
    pi_stay_work = (pi_stay./pi_denom_move_stay).*(pi_work_stay./pi_denom_work_stay);
    
    pi_move_work = (pi_move./pi_denom_move_stay).*(pi_work_move./pi_denom_work_move);
    
    pi_stay_rest = (pi_stay./pi_denom_move_stay).*(1 - (pi_work_stay./pi_denom_work_stay));
    
    pi_move_rest = (pi_move./pi_denom_move_stay).*(1 - (pi_work_move./pi_denom_work_move));
    
    test = [pi_stay_work, pi_move_work, pi_stay_rest, pi_move_rest];
   
    test(isnan(test)) = 0;
    
    policy_move(:,zzz,:) = test;
    
    %v_prime(:,zzz) = pi_stay.*v_stay + pi_move.*v_move +  pi_rest.*v_rest + pi_rest_move.*v_rest_move;
    v_prime(:,zzz) = b_smooth.*log(pi_denom_move_stay);
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
   value_fun_out(:,zzz,:) = [v_stay , v_move, v_rest, v_rest_move];
   
end

policy.move = policy_move;
policy.assets = policy_assets;
policy.value_fun = value_fun_out;
policy.vhat = v_prime;

