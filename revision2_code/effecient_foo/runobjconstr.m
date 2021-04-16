function [x,f,eflag,outpt] = runobjconstr(x0, params, tfp, opts)

if nargin == 1 % No options supplied
    opts = [];
end

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint


fun = @(xxx) objfun(xxx, params, tfp); % The objective function, nested below
cfun = @(xxx) constr(xxx, params, tfp); % The constraint function, nested below

LB = [(0.5.*ones(2,1)); (zeros(length(x0)-2,1))];

UB = [(1.5*ones(2,1)); (ones(length(x0)-2,1))];

%ga(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
% Call fmincon
%[x,f,eflag,outpt] = fmincon(fun,x0,[],[],[],[],(LB),(UB),cfun,opts);
[x,f] = ga(fun,length(x0),[],[],[],[],(LB),(UB),cfun,opts);

    function y = objfun(x, params, tfp)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myceq] = compute_effecient((x), params, tfp, 0);
            myc = [];
            xLast = x;
        end
        % Now compute objective function
        y = myf;
    end

    function [c,ceq] = constr(x, params, tfp)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myceq] = compute_effecient((x), params, tfp, 0);
            myc = [];
            xLast = x;
        end
        % Now compute constraint function
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end

end