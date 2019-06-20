
urban_prd = urban_prd + .001.*randn(size(urban_prd));

urban_prct = 20:20:80;
edges_z = prctile(urban_prd,urban_prct);
edges_z = [0 0.5778    0.6568    0.8023    1.1269 10.0000];


welfare_bin = zeros(length(edges_z)-1,1);
welfare_bin_cond = zeros(length(edges_z)-1,1);
migration_bin = zeros(length(edges_z)-1,1);
counts = zeros(length(edges_z)-1,1);

urban_bin = zeros(length(edges_z)-1,1);
expr_bin = zeros(length(edges_z)-1,1);

income_bin = zeros(length(edges_z)-1,1);
asset_bin  = zeros(length(edges_z)-1,1);

for xxx = 1:length(edges_z)-1
        
        urban_z_yes = edges_z(xxx) <= urban_prd & ... 
        urban_prd < edges_z(xxx+1) ;        
        
        test = (urban_z_yes == 1 );
        test_migrate = (test ==1) & (welfare_migrate_data(:,2) == 1);
        
        counts(xxx) = sum(test);
        
        welfare_bin(xxx) = mean(welfare_migrate_data(test,1));
        migration_bin(xxx) = mean(welfare_migrate_data(test,2));
        welfare_bin_cond(xxx) = prctile(welfare_migrate_data(test_migrate,1),90);
        
        urban_bin(xxx) = mean(urban_prd(test_migrate));
        expr_bin(xxx) = mean(expr_prd(test_migrate));
        
        income_bin(xxx) = mean(welfare_migrate_data(test,3));
        asset_bin(xxx) = mean(welfare_migrate_data(test,4));
end

disp(sum(counts(:)))
%disp('Welfare by Income Quintile')
%disp(100.*welfare_bin)