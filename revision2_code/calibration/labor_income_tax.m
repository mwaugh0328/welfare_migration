function [aftertax, tax, production] = labor_income_tax(laborincome, tax)

% if location == 1 %ubran guys pay the tax

production = laborincome;

aftertax = tax.rate.*(laborincome).^(1-tax.prog);

tax = laborincome - aftertax;
    
% else %rural guys DONT pay the tax
%     
%     production = laborincome;
% 
%     aftertax = laborincome;
% 
%     tax = 0;
% end

