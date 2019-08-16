clear
warning('off','stats:regress:RankDefDesignMat');

load min_calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First compute the competitive equillibrium and infer seasonal
% productivity...

[rural_labor_units, urban_labor_units] = compute_outcomes_CE(cal,[], 1);

alpha = 0.91;

monga_productivity = rural_labor_units(3)./(alpha.*rural_labor_units(1).^(alpha-1));

not_monga_productivity = (1./rural_labor_units(3))./(alpha.*rural_labor_units(2).^(alpha-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct wages

wage_monga = monga_productivity.*alpha.*(rural_labor_units(1)).^(alpha-1);

wage_not_monga = not_monga_productivity.*alpha.*(rural_labor_units(2)).^(alpha-1);

wages = [wage_monga; wage_not_monga];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the elasticity of wages w.r.t. labor...

wage_adj_productivity_not_monga = wage_not_monga./not_monga_productivity;

wage_adj_productivity_monga = wage_monga./monga_productivity;

change_wage = log(wage_adj_productivity_monga)-log(wage_adj_productivity_not_monga);

change_labor = log(rural_labor_units(4))- log(rural_labor_units(5));

wage_elasticity = change_wage./change_labor;

disp('Wage Elasticity')
disp(wage_elasticity)

labor_units_old = [rural_labor_units(1), rural_labor_units(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accounting_cal = cal;
accounting_cal(2) = 0.15;


disp('Calculating the Equillibrium')
n_iters = 50;
relax = 0.25;

for xxx = 1:n_iters

    [labor_units_new, ~] = compute_outcomes_CE(accounting_cal, wages, 0);
    
        
    if norm(log(labor_units_old) - log(labor_units_new)) < 10^-2
        
        disp('Converged')
        
        [rural_labor_units, urban_labor_units] = compute_outcomes_CE(accounting_cal, wages, 1);
        
        break
    
    end

    disp(norm(log(labor_units_old) - log(labor_units_new)) )
        
    labor_units_old = relax.*labor_units_new + (1-relax).*labor_units_old;
    
    wage_monga = monga_productivity.*alpha.*(labor_units_old(1)).^(alpha-1); 

    wage_not_monga = not_monga_productivity.*alpha.*(labor_units_old(2)).^(alpha-1); 

    wages = [wage_monga; wage_not_monga];
    
end




