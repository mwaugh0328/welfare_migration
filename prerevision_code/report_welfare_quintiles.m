
income_assets(:,[1,2]) = income_assets(:,[1,2]) + .01.*randn(size(income_assets(:,[1,2])));
% Add just a bit of noise to smooth stuf out...

income_prct = 20:20:80;
edges_income = prctile(income_assets(:,1),income_prct);
edges_income = [0, edges_income, 10];


welfare_bin = zeros(length(edges_income)-1,1);
welfare_bin_cond = zeros(length(edges_income)-1,1);
migration_bin = zeros(length(edges_income)-1,1);
counts = zeros(length(edges_income)-1,1);

income_gain_bin = zeros(length(edges_income)-1,1);
cons_gain_bin = zeros(length(edges_income)-1,1);
urban_bin = zeros(length(edges_income)-1,1);
expr_bin = zeros(length(edges_income)-1,1);

for xxx = 1:length(edges_income)-1
        
        income_yes = edges_income(xxx) <= income_assets(:,1) & ... 
        income_assets(:,1) < edges_income(xxx+1) ;        
        
        test = (income_yes == 1 );
        test_migrate = (test ==1) & (income_assets(:,4) == 1);
        
        counts(xxx) = sum(test);
        
        welfare_bin(xxx) = mean(income_assets(test,3));
        migration_bin(xxx) = mean(income_assets(test,4));
        welfare_bin_cond(xxx) = mean(income_assets(test_migrate,3));
        
        income_gain_bin(xxx) = mean(income_gain(test_migrate));
        cons_gain_bin(xxx) = mean(cons_gain(test_migrate));
        urban_bin(xxx) = mean(urban_prd(test_migrate));
        expr_bin(xxx) = mean(expr_prd(test_migrate));

end

disp(sum(counts(:)))
%disp('Welfare by Income Quintile')
%disp(100.*welfare_bin)
