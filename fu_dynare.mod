close all;

%----------------------------------
%       ENDOGENOUS VARIABLES
%----------------------------------

var C, N, U, V_Cf, i, V_YPh, V_P, C_f W, P, Z, C_h, P_h, Y, Pi_hat_W   V_W V_W_tilde W_tilde DW B ; 

%----------------------------------
%         EXOGENOUS SHOCKS
%----------------------------------

varexo eps_z;

%----------------------------------
%           PARAMETERS
%----------------------------------

parameters alpha eta transferCoef beta epsilon_W sigma rho varphi  chi  gamma P_f omega_W k_W C_f_param C_mean_param iferg i_star p tau Z_volatility ; 

%----------------------------------
%           CALIBRATION
%----------------------------------
    beta = 0.99;
    chi=1;
    sigma = 2;
    varphi=3;
    epsilon_W = 6;
    rho=0.95;
    eta=1.01;
    omega_W=.87;
    P_f=0;
    p = -0.001;
    i_star = 1/beta-1;
    load parameterfile;
    set_param_value('gamma',gamma);
    set_param_value('eta',  eta);
    set_param_value('alpha',alpha);
    set_param_value('omega_W',omega_W);
    set_param_value('C_f_param',C_f_param);
    set_param_value('C_mean_param',C_mean_param);
    set_param_value('i_star',i_star)
    set_param_value('iferg',iferg);
    set_param_value('transferCoef',transferCoef);
    set_param_value('Z_volatility',Z_volatility);
    set_param_value('rho',rho);
    k_W =(1-omega_W)*(1-beta*omega_W)/omega_W;
    tau=(epsilon_W-1)/epsilon_W;

%-----------------------------------
%               MODEL
%-----------------------------------
model;

% Relative demand
C_h=(P_h-P)*(-eta)+C;              

% Demand for foreign goods
C_f=(P_f-P)*(-eta)+C;         

% Aggregated foreign goods       
V_Cf=exp(C_f);      

% Consumer Price Index                
exp(P)^(1-eta)=((1-alpha)*exp(P_h)^(1-eta)+alpha*exp(P_f)^(1-eta)); 

% Production Function
Y=Z+N;                 

% Demand for Home Output    
exp(Y)=(1-alpha)*exp(C_h)+alpha*(exp(P_h)/exp(P_f))^(-gamma)*(iferg*C_f_param+(1-iferg)*exp(steady_state(C_f)));

% National Budget Constraint in Transfer Union
@#if specification==1
    exp(C+P)=exp(Y+P_h)+transferCoef*Z;     
@#endif

% National Budget Constraint in Autarky
@#if specification==2
    exp(C+P)=exp(Y+P_h);     
@#endif

% National Budget Constraint in Complete Markets
@#if specification==3
    exp(C)=((1-iferg)*exp(steady_state(Y+P_h-(1-1/sigma)*P))+iferg*C_mean_param)/exp(P)^(1/sigma);     
    V_YPh=exp(Y+P_h);
    V_P=exp((1-1/sigma)*P);
@#else 
    V_YPh=exp(steady_state(Y)+steady_state(P_h));
    V_P=exp((1-1/sigma)*steady_state(P));
@#endif

% National Budget Constraint in Incomplete Markets (Bond Economy)
@#if specification==4
    exp(Y)*exp(P_h)+(1+i)*B(-1)=exp(C)*exp(P)+B;
    i=i_star+p*B;
    exp(C)^(-sigma)=exp(C(+1))^(-sigma)*beta*(1+i)*exp(P)/exp(P(+1));
@#else
    i=i_star;
    B=0;
@#endif

% Wage inflation
Pi_hat_W=W-W(-1);                                 
                  
%Calvo Wages
exp(V_W)=exp(N)^(1+varphi)+beta*omega_W*exp(V_W(+1));
exp(V_W_tilde)=exp(C)^(-sigma)*exp(N)/exp(P)+beta*omega_W*exp(V_W_tilde(+1));
exp(W_tilde)=chi*epsilon_W/(epsilon_W-1)*tau*exp(V_W)/exp(V_W_tilde);
exp(W)^(1-epsilon_W)=(1-omega_W)*exp(W_tilde)^(1-epsilon_W)+omega_W*exp(W(-1))^(1-epsilon_W);
exp(DW)=(1-omega_W)*(exp(W_tilde)/exp(W))^(-epsilon_W*(1+varphi))+omega_W*(exp(W(-1))/exp(W))^(-epsilon_W*(1+varphi))*exp(DW(-1));

% Firm FOC w.r.t. to output
P_h=W-Z;      

% Productivity Shock Process
Z=rho*Z(-1)+Z_volatility*eps_z;

% Household Utility
U=exp(C)^(1-sigma)/(1-sigma)-chi*exp(N)^(1+varphi)/(1+varphi)*exp(DW);
end;

%check;
initval;
Z=0;
V_W=log(1/(1-beta*omega_W));
V_W_tilde=log(1/(1-beta*omega_W));
V_YPh=1;
V_P=1;

end;
steady;

shocks;
var eps_z; stderr 1;
end;

@#if option_hp==0
    stoch_simul(order=2, nograph,irf=0,nocorr, ar=0,nofunctions, noprint) C N U V_Cf i V_YPh V_P Y ;
@#else
    stoch_simul(order=2, nograph,irf=0,nocorr, ar=0,nofunctions, noprint,hp_filter=1600) C N U V_Cf i V_YPh V_P Y ;
@#endif


