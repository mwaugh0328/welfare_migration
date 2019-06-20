load migrate_rates
close all

years = migratesS1(:,1);

data = migratesS1(:,4);


error_one = - migratesS1(:,5).*1.96;
error_two = + migratesS1(:,5).*1.96;

data_isnan = isnan(data);


H = shadedErrorBar(years(~data_isnan), data(~data_isnan),[error_two(~data_isnan)],'b');
hold

D = plot(years(~data_isnan), m_rates(~data_isnan),'r-' );
xlim([2008 2013])
set(gca,'XTick',[2008:2013])
box off
ylabel('Migration Rate: Treatment Minus Control')

legend([H.mainLine H.patch D],'Data','95-5 Confidence Interval (Data)', 'Model Prediction')
legend('boxoff')