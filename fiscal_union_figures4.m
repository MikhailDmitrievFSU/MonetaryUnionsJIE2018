
clear all
close all
tic

%-------------------------------------------------------
%   DEFINE THE MIN AND MAX VALUES OF ALPHA AND GAMMA
%-------------------------------------------------------
options=2;     %1 - transfer union 2 - autarky 3 - complete 4 - incomplete
N_eta=25;
N_omega_W=25;
omega_W_min=0.5;
omega_W_max=0.95;
eta_min=0.6;
eta_max=5.1;
eta=1.05;
gamma=3.75;
omega_W=0.87;
alpha=0.35;
chi=1;
varphi=3;
sigma=2;
Z_volatility=0.01;
OMEGA_W=[omega_W_min:(omega_W_max-omega_W_min)/(N_eta-1):omega_W_max]; 
ETA=[eta_min:(eta_max-eta_min)/(N_omega_W-1):eta_max];

%-------------------------------------------------------
%        CREATE VECTORS OF ZEROS FOR REPLACEMENT
%-------------------------------------------------------

C_loss=zeros(N_eta,N_omega_W);
transferOptimal=zeros(N_eta,N_omega_W);
LABOR_SS=zeros(N_eta,N_omega_W);
C_SS=zeros(N_eta,N_omega_W);


%-------------------------------------------------------
%               LOOP OVER ALPHA AND ETA
%-------------------------------------------------------


for options=1:2
    iferg=0;
    C_f_param=0;
    C_mean_param=0;
    eta=2;
    alpha=0.35;
    transferCoef=0;
    omega_W=0.87;
    save parameterfile alpha eta gamma chi varphi sigma iferg C_f_param C_mean_param transferCoef omega_W

    if options==1
        dynare fu_dynare -Dspecification=1 -Doption_hp=0 noclearall 
    elseif options==2
        dynare fu_dynare -Dspecification=2 -Doption_hp=0 noclearall 
    elseif options==3
        dynare fu_dynare -Dspecification=3 -Doption_hp=0 noclearall 
    elseif options==4
        dynare fu_dynare -Dspecification=4 -Doption_hp=0 noclearall 
    end 
    C_SS=oo_.steady_state(1,1);
    LABOR_SS=oo_.steady_state(2,1);
        
    V_Cf_0=oo_.steady_state(4,1); 
    i_0=oo_.steady_state(5,1); 
    C_mean_param_0=oo_.steady_state(6,1)/oo_.steady_state(7,1);  
    eps=1;

    for j=1:N_omega_W
        for i=1:N_eta
            close all
            set_param_value('eta',  ETA(1,i));
            set_param_value('omega_W',OMEGA_W(1,j));
            set_param_value('iferg',1)
    
            if options==1
                transferCoef0=(1-gamma)+(1-ETA(1,i))*(1-alpha);
                optionsOpt = optimoptions('fminunc','Algorithm','quasi-newton','StepTolerance',1e-6,'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',200,'Display','final-detailed','FiniteDifferenceStepSize',1e-5);
                objective=@(x)autarky(x,V_Cf_0,i_0,C_mean_param_0, LABOR_SS, C_SS,sigma, chi, varphi);
                [transferOptimal(i,j),C_loss(i,j)]=fminunc(objective, transferCoef0,optionsOpt);
                [ETA(1,i) OMEGA_W(1,j)]
                C_loss(i,j)
            elseif options>1
                C_loss(i,j) =autarky(0,V_Cf_0,i_0,C_mean_param_0, LABOR_SS, C_SS,sigma, chi, varphi);
            end
        end
    end
    C_loss_cell{options}=C_loss;
end

figure( 'Name', 'Gains from access to fiscal transfers under autarky');
surf(OMEGA_W,ETA, 100*(C_loss_cell{2}-C_loss_cell{1}),'FaceAlpha',0.75);
legend('off');
%colormap('jet');
xlabel('Calvo Wage Rigidity (\theta_W)');
ylabel('Trade Elasticity (\eta)');
zlabel('Gains In Permanent Consumption (%)');
title('')
hgexport(gcf, 'graph_eta_thetaW_axes.eps', hgexport('factorystyle'), 'Format', 'eps');



toc
delete *.log *_dynamic.m *_static.m  *_set_auxiliary_variables.m *_results.mat *_results2.mat
delete  fu_dynare.m 
rmdir('fu_dynare', 's')

function f=autarky(transferCoef_loc,V_Cf_loc,i_loc,C_mean_param_loc, LABOR_SS, C_SS, sigma, chi, varphi )
        global oo_ M_ options_
        V_Cf_old=V_Cf_loc;
        i_old=i_loc;
        C_mean_param_old=C_mean_param_loc;
        set_param_value('transferCoef',transferCoef_loc);
        iter=0;
        set_param_value('iferg',1);
        eps=1;
        while eps>1e-8% iter<1
            set_param_value('C_f_param',V_Cf_old);
            set_param_value('i_star',i_old);
            set_param_value('C_mean_param',C_mean_param_old);
            steady;
            var_list_2 = char('C','N','U','V_Cf','i','V_YPh','V_P');
            stoch_simul(var_list_2);
            eps_vec=[abs(oo_.mean(4,1)-V_Cf_old),abs(oo_.mean(5,1)-i_old),abs(oo_.mean(6,1)/oo_.mean(7,1)-C_mean_param_old)];
            eps=max(eps_vec);
            V_Cf_old=oo_.mean(4,1);
            i_old=oo_.mean(5,1); 
            C_mean_param_old=oo_.mean(6,1)/oo_.mean(7,1);
            iter=iter+1;
        end
        iter;
        eps;
        U=oo_.mean(3,1);
        f =1- ((1-sigma).*( U + chi/(1+varphi).*exp(LABOR_SS).^(1+varphi))).^(1./(1-sigma))./exp(C_SS);
        set_param_value('C_f_param',V_Cf_loc)
        steady;
       
end






