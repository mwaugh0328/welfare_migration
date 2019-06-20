
income_assets(:,[1,2]) = income_assets(:,[1,2]) + .001.*randn(size(income_assets(:,[1,2])));
urban_prd = urban_prd + .001.*randn(size(urban_prd));

urban_prct = 20:20:80;
edges_z = prctile(urban_prd,urban_prct);
edges_z = [0    0.5786    0.6649    0.8009    1.1103   10.0000];
% Add just a bit of noise to smooth stuf out...

income_prct = 20:20:80;
edges_income = prctile(income_assets(:,1),income_prct);
edges_income = [0    0.6152    0.8267    1.0687    1.4156   10.0000];


welfare_bin = zeros(length(edges_income)-1,length(edges_z)-1);
welfare_bin_cond = zeros(length(edges_income)-1,length(edges_z)-1);
migration_bin = zeros(length(edges_income)-1,length(edges_z)-1);
counts = zeros(length(edges_income)-1,length(edges_z)-1);

income_gain_bin = zeros(length(edges_income)-1, length(edges_z)-1);
cons_gain_bin = zeros(length(edges_income)-1, length(edges_z)-1);
urban_bin = zeros(length(edges_income)-1, length(edges_z)-1);
expr_bin = zeros(length(edges_income)-1, length(edges_z)-1);

for xxx = 1:length(edges_income)-1
        
        income_yes = edges_income(xxx) <= income_assets(:,1) & ... 
        income_assets(:,1) < edges_income(xxx+1) ;        
        
    for zzz = 1:length(edges_z)-1
        
        urban_z_yes = edges_z(zzz) <= urban_prd & ... 
        urban_prd < edges_z(zzz+1); 
    
    
        test = (income_yes == 1 ) & (urban_z_yes == 1 );
        test_migrate = (test == 1) & (income_assets(:,4)==1);
        counts(xxx,zzz) = sum(test);
        
        welfare_bin(xxx,zzz) = mean(income_assets(test,3));
        welfare_bin_cond(xxx,zzz) = mean(welfare_migrate_data(test_migrate,1));
        migration_bin(xxx,zzz) = mean(income_assets(test,4));
        expr_bin(xxx,zzz) =  mean(expr_prd(test));
    end

end

disp(sum(counts(:)))
%disp('Welfare by Income Quintile')
%disp(100.*welfare_bin)
