%David Lagakos -- RA: Alejandro Nakab
%Standard Errors Calculation for LMW 2019 ECMA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%% Elasticities for parameters with respect to Moments
% order of moments from compute_outcome.m is:
% cons_growth P_noassets mig_rate cons_gainols  mig_inc  mig_inc_y2 cons_gainlate wage_gap p_rural var_logwage 
% order of parameters: 
% order_table = [2,5,6,7,8,3,1,4]=[theta,ubar,lambda,pi_prob,gamma,A_u,sigma_s,rho // sigma_rc, sigma_ui]
% order_moments = [10, 4, 5, 9, 6, 7, 8, 1, 2, 3]=[cons_growth P_noassets mig_rate cons_gainols  mig_inc  mig_inc_y2 cons_gainlate wage_gap p_rural var_logwage];
% 1- Wage gap
% 2- The rural share
% 3- The urban variance... note that this is position number 3 (see below)
% 4- Fraction with no liquid assets
% 5- seasonal migration in control
% 6- increase in r1 (22 percent)
% 7- increase in r2 (9.2 percent)
% 8- LATE estiamte
% 9- OLS estimate
% 10- Standard deviation of consumption growth.

clear
T=6;%number of changes from 1.0025 to 1.0150
Matrix_M_new=zeros(10,10*(T+1));

for j=1:T
    load('calibration_allmoments_new.mat')

    params = new_cal;

    n_params = length(params); % how many paramters we need to do...

    eps=1+0.0025*(j); %(1.05 before) This is the change. One issue is that some of this stuff was
    % not changing that much, this is a question of how accuratly we are
    % solving it, here is an area of investigation.

    els_moments = zeros(10,10);

    for xxx = 1:n_params

        disp(xxx)

        cal_eps_for = params;
        cal_eps_bak = params;

        cal_eps_for(xxx) = params(xxx).*eps; % forward
        cal_eps_bak(xxx) = params(xxx)./eps; % backward

        change_cal = (cal_eps_for(xxx))-(cal_eps_bak(xxx));

        moments_for = calibrate_model(cal_eps_for,2); % moments forward
        moments_bak = calibrate_model(cal_eps_bak,2); % moments backward

        change_moments = (moments_for) - (moments_bak);

        els_moments(xxx,:) = change_moments'./change_cal; % compute change.

    % So each row is a parameter, the each column is the moment

    end

    %%%%%%%%%
    % aggregate_moments = [1.80, 0.63, 0.68, 0.47];
    % experiment_hybrid = [0.36, 0.22, 0.092, 0.30, 0.10,  0.40];
    % 

    order_table = [2,5,6,7,8,3,1,4];
    order_moments = [5, 6, 7, 9, 8, 1, 2, 4];
    test = els_moments(order_table,order_moments);
    round(test,2)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to have same order as var-covar matrix:
    order_table = [2,5,6,7,8,3,1,4];
    order_moments = [10, 4, 5, 9, 6, 7, 8, 1, 2, 3];
    test = els_moments(order_table,order_moments);
    Matrix_M = [test' zeros(10,2)];
    Matrix_M(1,9)=1;
    Matrix_M(10,10)=1;
    Matrix_M_new(:,j*10+1:j*10+10)=Matrix_M;
end

M_mean=(1/T)*(Matrix_M_new(:,11:20)+Matrix_M_new(:,21:30)+Matrix_M_new(:,31:40) ...
          +Matrix_M_new(:,41:50)+Matrix_M_new(:,51:60)+Matrix_M_new(:,61:70)) ;
Matrix_M_new(:,1:10)=M_mean;

%Example of an M matrix: below under extras.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Variance Covariance Matrix

                   %cons_growth P_noassets    mig_rate  cons_gainols    mig_inc     mig_inc_y2  cons_gainlate
Exp_mom_varcovar=  [0.00128584	0.00001115	-0.00003059  0.00019928	  0.00002116	0.00004212	0.00000031       % cons_growth
                    0.00001115	0.00012712	0.00000817	 0.00001947	  -0.00001500	-0.00000289	-0.00003335      % P_noassets
                    -0.00003059	0.00000817	0.00069678	 -0.00000131  -0.00036251	-0.00017265	0.00034881       % mig_rate
                    0.00019928	0.00001947	-0.00000131  0.00199811	  -0.00001064	0.00001570	-0.00037858      % cons_gainols
                    0.00002116	-0.00001500	-0.00036251	 -0.00001064  0.00057068	0.00020858	-0.00031222      % mig_inc
                    0.00004212	-0.00000289	-0.00017265	 0.00001570   0.00020858	0.00059693	-0.00019477      % mig_inc_y2 
                    0.00000031	-0.00003335	0.00034881	 -0.00037858  -0.00031222	-0.00019477	0.00934914];     % cons_gainlate

                 
                    % wage_gap    p_rural   var_logwage             
Surv_mom_varcovar=  [ 0.020911   0.000063   0.002406    % wage_gap
                      0.000063   0.000046   -0.000027   % p_rural
                      0.002406   -0.000027  0.001946];  % var_logwage                  

zeros_7x3=zeros(7,3); 
zeros_3x7=zeros(3,7);

V=[Exp_mom_varcovar  zeros_7x3
   zeros_3x7         Surv_mom_varcovar]
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Weighting Matrix

W= eye( 10 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Gamma Matrix

g=1; %gamma parameter
g=[ones(1,7)*1/g ones(1,3)];
GAM=eye(10).*g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SD Calculations 
%%% For the mean matrix (levels off due to the 1.5% simulation)
M=M_mean;
SD_mat= inv(M'*W*M) * M' * W * GAM * V * GAM * W * M * inv(M'*W*M)
SD=diag(SD_mat)';
mean_sd_mean=[params(order_table) 0.37 0.31; SD];

%%% For the 0.5% matrix (smallest so far)
M=Matrix_M_new(:,11:20);
SD_mat= inv(M'*W*M) * M' * W * GAM * V * GAM * W * M * inv(M'*W*M);
SD=diag(SD_mat)';
mean_sd_05=[params(order_table) 0.37 0.31; SD];

%%% for all:
sd_all=zeros(7,10);
for j=0:6
M=Matrix_M_new(:,j*10+1:j*10+10);
SD_mat= inv(M'*W*M) * M' * W * GAM * V * GAM * W * M * inv(M'*W*M)
SD=diag(SD_mat)';
sd_all(j+1,:)=SD;
end
mean_sd_all=[params(order_table) 0.37 0.31; sd_all];


















%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% EXTRAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


             %theta               ubar               lambda               pi_prob                gamma                A_u             sigma_s                 rho         sigma_rc  sigma_ui
M=   [ -0.0206323396970548	0.0446789192079261	0.0120378273753058    -0.00583267092339692	0.0742013811631713	0.00617415639497008	0.121882698796568	-0.460028166581215  	1       0       % cons_growth
        -0.140158624408723	-0.0953100790789674	0.0754185672147121	  -0.0369400707530855	-0.0935474485497274	-0.256173770085839	-1.02216559983323	3.24247567273733        0       0       % P_noassets
        -2.42342106251463	-1.90256998771693	-0.320080752961127	  0.909710809079324     -0.679108422009364	0.615712548286370	-0.126189199358721	0.356216831131724       0       0       % mig_rate
        0.680503262035387	0.167092971298870	-0.115897516656599	  0.180841117676576	    0.511786122199301	-0.121435435280994	-0.0453488232386904	0.159009346391686       0       0       % cons_gainols
        -0.603650359144762	-0.169376364238634	0.107350472553348	  -0.0506694637163164	-0.102107239654800	0.360497797115580	0.212489370617591	-0.411684057893901      0       0       % mig_inc
        -0.255017897413134	-0.0502351821695374	-0.108377835176062	  0.0504231965779612	0.0395828205933962	0.130129535277801	0.0803376904122831	-0.101359609620140      0       0       % mig_inc_y2
        -0.164837252539745	0.230882318919817	0.000846655105998828  -0.0589016204621054	-0.0208491822826699	0.105056867950595	-0.0725723423850904	0.273083248832135       0       0       % cons_gainlate
        4.57728285993360	-0.145263405501969	0.196838551366081	  0.0850838982147289	0.752882421867080	-0.403360249557969	-0.529640348548826	0.625868349941994       0       0       % wage_gap 
        0.865842828296977	-0.329482685784942	-0.0694189610817254	  0.209939567641733     0.229589800446487	-0.603941372403899	-0.174079031439711	0.299797211519787       0       0       % p_rural
        1.28504185646498	-0.0127722118971912	-0.0692925240820247	  0.0191042584192962	0.299834821319474	0.0388075095573367	0.261051683487046	-0.556282376059827      0       1];     % var_logwage

%check=Matrix_M-M;
%M=Matrix_M;

M2=M';
SD_mat2= inv(M2'*W*M2) * M2' * W * GAM * V * GAM * W * M2 * inv(M2'*W*M2)
SD2=diag(SD_mat2)'

% Means to compare:
% theta  ubar  lambda  pi_prob  gamma  A_u sigma_s  rho  sigma_rc  sigma_ui

mean_sd=[params(order_table) 1 1; SD];



%%%%%%%%%%%%%%%%%%%%%%%

elast_bar_u = [ones(1,10)*0.5]; % cons_growth                                           
elast_lambda = [ones(1,10)*0.2]; % P_noassets
elast_pi = [ones(1,10)*0.53]; % mig_rate
elast_theta = [ones(1,10)*0.54]; % cons_gainols
elast_gamma = [ones(1,10)*0.55]; % mig_inc
elast_A_u = [ones(1,10)*0.56]; % mig_inc_y2
elast_sigma_s = [ones(1,10)*4]; % cons_gainlate
elast_rho = [ones(1,10)*0.58]; % wage_gap 
elast_sigma_rc = [0 0 0 0 0 0 0 0 0 0]; %p_rural
elast_sigma_ui = [ones(1,10)*1]; % var_logwage
%}

