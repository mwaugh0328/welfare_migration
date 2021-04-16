clear
warning('off','stats:regress:RankDefDesignMat');

addpath('../calibration')

load calibration_final

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First compute the competitive equillibrium and infer seasonal
% productivity...

%[average_rural_labor_units_monga, average_rural_labor_units_not_monga, seasonal_factor, labor_monga, labor_not_monga];
%    urban_data = [average_urban_labor_units_monga, average_urban_labor_units_not_monga];
tic
[rural, urban, vfun, params] = compute_outcomes_prefshock_GE(exp(new_val),[], [],[], []);
toc

alpha = 0.845;

monga_productivity = (rural.seasonal_factor)./(alpha.*(rural.labor_units_monga).^(alpha-1))

not_monga_productivity = (1./rural.seasonal_factor)./(alpha.*rural.labor_units_not_monga.^(alpha-1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct wages

wage_monga = monga_productivity.*alpha.*(rural.labor_units_monga).^(alpha-1);

wage_not_monga = not_monga_productivity.*alpha.*(rural.labor_units_not_monga).^(alpha-1);

wages = [wage_monga; wage_not_monga];

mtest = params.means_test;

save wages.mat wages vfun mtest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the elasticity of wages w.r.t. labor...


wage_increase = 100*(log((rural.expr_labor_units).^(alpha-1)) - log((rural.cntr_labor_units).^(alpha-1)));
% Compare the implied change across the two groups...

% From ACM: For every extra 10% of the landless population that emigrates, wages increase by 2.2%
% With a 22 percent increase in emigration between the two...this means that 
% Agricultural wages are therefore predicted to increase by 2.2 percent *
% 2.2 = 4.84
disp('')
disp('')
disp('Wage Change from ACM, Wage Change in Model')
disp([4.48, wage_increase])
%disp([rural.labor_monga, rural.labor_not_monga])

labor_units_old = [rural.labor_units_monga, rural.labor_units_not_monga];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computes eq. gives us some baseline statistics....
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('')
disp('First Compute the Baseline: Note welfare here should be 0 as nothing happend')
compute_outcomes_prefshock_GE(exp(new_val), wages, 0, vfun,1);

wages_old = wages;

disp('')
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('')
disp('Now compute counterfactual holding wages fixed...')

compute_outcomes_prefshock_GE(exp(new_val), wages, params.means_test, vfun,1);

disp('')
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is where the counter-factual eq would be...
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('')
disp('Calculating the Equillibrium')
n_iters = 10;
relax = 0.25;

for xxx = 1:n_iters

    [ruralnew, ~] = compute_outcomes_prefshock_GE(exp(new_val), wages, params.means_test,[],[]);
    
    labor_units_new = [ruralnew.labor_units_monga, ruralnew.labor_units_not_monga];
    
        
    if norm(log(labor_units_old) - log(labor_units_new)) < 10^-2
        
        disp('Converged')
        disp(labor_units_new)
                
        break
    
    end

    disp(norm(log(labor_units_old) - log(labor_units_new)) )
        
    labor_units_old = relax.*labor_units_new + (1-relax).*labor_units_old;
    
    wage_monga = monga_productivity.*alpha.*(labor_units_old(1)).^(alpha-1); 

    wage_not_monga = not_monga_productivity.*alpha.*(labor_units_old(2)).^(alpha-1); 

    wages = [wage_monga; wage_not_monga];
    
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('')
disp('Change in Wages')
disp(wages./wages_old)
disp('')
disp('Now compute counterfactual welfare with wages changing...')
disp('')

compute_outcomes_prefshock_GE(exp(new_val), wages, params.means_test, vfun,1);

rmpath('../calibration')


