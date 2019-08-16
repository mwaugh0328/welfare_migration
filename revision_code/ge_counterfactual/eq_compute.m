clear
warning('off','stats:regress:RankDefDesignMat');

addpath('../calibration')

load calibration_highgrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First compute the competitive equillibrium and infer seasonal
% productivity...

%[average_rural_labor_units_monga, average_rural_labor_units_not_monga, seasonal_factor, labor_monga, labor_not_monga];
%    urban_data = [average_urban_labor_units_monga, average_urban_labor_units_not_monga];

[rural, urban, params] = compute_outcomes_prefshock_GE(exp(new_val),[], []);

alpha = 0.91;

monga_productivity = (rural.seasonal_factor)./(alpha.*(rural.labor_units_monga).^(alpha-1));

not_monga_productivity = (1./rural.seasonal_factor)./(alpha.*rural.labor_units_not_monga.^(alpha-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct wages

wage_monga = monga_productivity.*alpha.*(rural.labor_units_monga).^(alpha-1);

wage_not_monga = not_monga_productivity.*alpha.*(rural.labor_units_not_monga).^(alpha-1);

wages = [wage_monga; wage_not_monga];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the elasticity of wages w.r.t. labor...

wage_adj_productivity_not_monga = wage_not_monga./not_monga_productivity;

wage_adj_productivity_monga = wage_monga./monga_productivity;

change_wage = log(wage_adj_productivity_monga)-log(wage_adj_productivity_not_monga);

change_labor = log(rural.labor_monga)- log(rural.labor_not_monga);

wage_elasticity = change_wage./change_labor;

disp('Wage Elasticity')
disp(wage_elasticity)
disp([rural.labor_monga, rural.labor_not_monga])

labor_units_old = [rural.labor_units_monga, rural.labor_units_not_monga];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is where the counter-factual eq would be...

disp('Calculating the Equillibrium')
n_iters = 10;
relax = 0.25;

for xxx = 1:n_iters

    [ruralnew, ~] = compute_outcomes_prefshock_GE(exp(new_val), wages, params.means_test);
    
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

rmpath('../calibration')


