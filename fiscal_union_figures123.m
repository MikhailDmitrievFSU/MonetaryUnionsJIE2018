
clear all
close all
tic

%-------------------------------------------------------
%   DEFINE THE MIN AND MAX VALUES OF ALPHA AND GAMMA
%-------------------------------------------------------
options=2;     %1 - transfer union 2 - autarky 3 - complete 4 - incomplete
N_eta=25;
N_alpha=25;
eta_min=0.6;
eta_max=5.1;
alpha_min=0.20;
alpha_max=0.89;
gamma=3.75;
omega_W=0.87;
eta=1.05;
chi=1;
varphi=3;
sigma=2;
Z_volatility=0.01;
ETA=[eta_min:(eta_max-eta_min)/(N_eta-1):eta_max]; 
ALPHA=[alpha_min:(alpha_max-alpha_min)/(N_alpha-1):alpha_max];

%-------------------------------------------------------
%        CREATE VECTORS OF ZEROS FOR REPLACEMENT
%-------------------------------------------------------

C_loss=zeros(N_eta,N_alpha);
transferOptimal=zeros(N_eta,N_alpha);
LABOR_SS=zeros(N_eta,N_alpha);
C_SS=zeros(N_eta,N_alpha);


%-------------------------------------------------------
%               LOOP OVER ALPHA AND ETA
%-------------------------------------------------------


for options=1:4
    iferg=0;
    C_f_param=0;
    C_mean_param=0;
    eta=2;
    alpha=0.35;
    transferCoef=0;
    save parameterfile alpha eta gamma chi varphi sigma iferg C_f_param C_mean_param transferCoef

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

    for j=1:N_alpha
        for i=1:N_eta
            close all
            set_param_value('eta',  ETA(1,i));
            set_param_value('alpha',ALPHA(1,j));
            set_param_value('iferg',1)
    
            if options==1
                transferCoef0=(1-gamma)+(1-ETA(1,i))*(1-ALPHA(1,j));
                optionsOpt = optimoptions('fminunc','Algorithm','quasi-newton','StepTolerance',1e-6,'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',200,'Display','final-detailed','FiniteDifferenceStepSize',1e-5);
                objective=@(x)autarky(x,V_Cf_0,i_0,C_mean_param_0, LABOR_SS, C_SS,sigma, chi, varphi);
                [transferOptimal(i,j),C_loss(i,j)]=fminunc(objective, transferCoef0,optionsOpt);
                [ETA(1,i) ALPHA(1,j)]
                C_loss(i,j)
            elseif options>1
                C_loss(i,j) =autarky(0,V_Cf_0,i_0,C_mean_param_0, LABOR_SS, C_SS,sigma, chi, varphi);
            end
        end
    end
    C_loss_cell{options}=C_loss;
end

figure( 'Name', 'autarky');
surf(ALPHA,ETA, 100*(C_loss_cell{2}-C_loss_cell{1}),'FaceAlpha',0.75);
legend('off');
%colormap('jet');
xlabel('Trade Openness (\alpha)');
xlim([0.2 0.9]);
ylabel('Trade Elasticity (\eta)');
zlabel('Gains In Permanent Consumption (%)');
title('');
hgexport(gcf, 'autarky_gains.eps', hgexport('factorystyle'), 'Format', 'eps');



figure( 'Name', 'complete markets');
surf(ALPHA,ETA, 100*(C_loss_cell{3}-C_loss_cell{1}),'FaceAlpha',0.75);
legend('off');
%colormap('jet');
xlabel('Trade Openness (\alpha)');
xlim([0.2 0.9]);
ylabel('Trade Elasticity (\eta)');
zlabel('Gains In Permanent Consumption (%)');
title('');
hgexport(gcf, 'complete_gains.eps', hgexport('factorystyle'), 'Format', 'eps');


figure( 'Name', 'incomplete markets');
surf(ALPHA,ETA, 100*(C_loss_cell{4}-C_loss_cell{1}),'FaceAlpha',0.75);
legend('off');
%colormap('jet');
xlabel('Trade Openness (\alpha)');
xlim([0.2 0.9]);
ylabel('Trade Elasticity (\eta)');
zlabel('Gains In Permanent Consumption (%)');
title('');
hgexport(gcf, 'incomplete_gains.eps', hgexport('factorystyle'), 'Format', 'eps');



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
        while iter<1%eps>1e-8%eps>1e-7 %iter<5%eps>1e-7 
            set_param_value('C_f_param',V_Cf_old);
            set_param_value('i_star',i_old);
            set_param_value('C_mean_param',C_mean_param_old);
            steady;
            var_list_2 = char('C','N','U','V_Cf','i','V_YPh','V_P','Y');
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






