
clear all
close all
tic

%For table 2 set empirical_shocks=0, Data_specification=-1,
% and for rigidity omega_W=0.001/omega_W=0.75/omega_W=0.87

%For table 3 set empirical_shocks=0, Data_specification=1, 
% and for rigidity omega_W=0.001/omega_W=0.75/omega_W=0.87

%For table 4 set empirical_shocks=1, omega_W=0.75, and 
%Data_specification should change from -1 to 2 for all 4 columns

%options=3;     %1 - transfer union 2 - autarky 3 - complete 4 - incomplete
empirical_shocks=1;
Data_specification=2;  % -1 if Corbo without taking average 0 if Corbo with average adjustment 1 Imbs without adjustment 2 Imbs with adjustment
omega_W=0.75;

if empirical_shocks==1
    rho=0.99;
elseif empirical_shocks==0
    rho=0.95;
end

Names={
'Austria';
'Czech';
'Denmark';
'Finland';
'France';
'Germany';
'Greece';
'Hungary';
'Italy';
'Netherlands';
'Portugal';
'Slovakia';
'Spain';
'Sweden';
'UK';};

GDP_volatility_data_Corbo=[
0.0049
0.0095
0.0072
0.0100
0.0049
0.0058
0.0243
0.0102
0.0039
0.0050
0.0113
0.0143
0.0073
0.0089
0.0076
0.0000];

GDP_volatility_data_Corbo2=[
0.0134
0.0192
0.0154
0.0213
0.0099
0.0172
0.0245
0.0172
0.0136
0.0133
0.0131
0.0229
0.0125
0.0186
0.0135
0.0130];

GDP_volatility_data_Imbs=[
    0.0134020
    0.0213315
    0.0099074
    0.0172433
    0.0244679
    0.0171525
    0.0135894
    0.0130618
    0.0228667
    0.0124966
    0.0185846
    0.0135477
    0.0130268];

GDP_volatility_data_Imbs2=[
0.0049
0.0100
0.0049
0.0058
0.0243
0.0102
0.0039
0.0113
0.0143
0.0073
0.0089
0.0076
0.0000];

Elasticities_Corbo=[
0.55	4.5	3.8	; %Austria
0.70	3.4	3.8	; %Czech
0.40	3.3	3.4	; %Denmark
0.32	3.5	3.4	; %Finland
0.31	3.7	3.8	; %France
0.34	3.7	4.3	; %Germany
0.35	2.9	4.1	; %Greece
0.49	3.3	4.2	; %Hungary
0.22	3.2	3.2	; %Italy
0.62	3.5	3.5	; %Netherlands
0.41	3.3	3.9	; %Portugal
0.30	3.7	3.9	; %Slovakia
0.30	3.4	3.2	; %Spain
0.42	4.2	4.5	; %Sweden
0.37	2.9	3.0	; %UK
0.38	3.5	3.73	; %Euro Avg
]; 

Elasticities_Imbs=[
0.55	1.915	3.906   ; %Austria
     %Czech
	 %Denmark
0.32	3.475	3.393   ; %Finland
0.31	3.133	3.493   ; %France
0.34	3.525	3.803   ; %Germany
0.35	4.83	4.15    ; %Greece
0.49	2.418	3.072	; %Hungary
0.22	3.889	3.753   ; %Italy
	 %Netherlands
0.41	3.552	4.872	; %Portugal
0.30	3.246	2.738   ; %Slovakia
0.30	3.485	4.387	; %Spain
0.42	3.244	3.998	; %Sweden
0.37	3.175	3.561	; %UK
0.38	3.32     3.78     ; %Euro Avg
]; 

if Data_specification==0
    MatrixParameters=Elasticities_Corbo;
    GDP_volatility_data=GDP_volatility_data_Corbo;
end
if Data_specification==1
    MatrixParameters=Elasticities_Imbs;
    GDP_volatility_data=GDP_volatility_data_Imbs;
end
if Data_specification==-1
    MatrixParameters=Elasticities_Corbo;
    GDP_volatility_data=GDP_volatility_data_Corbo2;
end
if Data_specification==2
    MatrixParameters=Elasticities_Imbs;
    GDP_volatility_data=GDP_volatility_data_Imbs2;
end    
volatility_technology_shocks=0.01*ones(length(MatrixParameters),1);

gamma=3.75;
eta=1.05;
varphi=3;
sigma=2;
Z_volatility=0.01;
C_f_param=0;
C_mean_param=0;
transferCoef=0;
alpha=0.5;
iferg=0;

save parameterfile transferCoef alpha eta gamma omega_W Z_volatility C_f_param C_mean_param iferg

C_loss=zeros(length(MatrixParameters),1);
transferOptimal=zeros(length(MatrixParameters),1);
LABOR_SS=zeros(length(MatrixParameters),1);
C_SS=zeros(length(MatrixParameters),1);

if empirical_shocks==1

    save parameterfile transferCoef alpha eta gamma Z_volatility iferg omega_W rho
    dynare fu_dynare -Dspecification=4 -Doption_hp=1 noclearall 
    V_Cf_0=oo_.steady_state(4,1); 
    i_0=oo_.steady_state(5,1); 
    C_SS=oo_.steady_state(1,1);
    LABOR_SS=oo_.steady_state(2,1);    
    V_Cf_0=oo_.steady_state(4,1); 
    i_0=oo_.steady_state(5,1); 
    C_mean_param_0=oo_.steady_state(6,1)/oo_.steady_state(7,1);  
 
    for i=1:length(MatrixParameters)
        set_param_value('gamma',  MatrixParameters(i,3));
        set_param_value('eta',  MatrixParameters(i,2));
        set_param_value('alpha',MatrixParameters(i,1));
        set_param_value('iferg',1)
        compute_variance=1;
        volatility_simulated_GDP(i,1)=autarky(0,V_Cf_0,i_0,C_mean_param_0, LABOR_SS, C_SS,sigma, chi, varphi,compute_variance);
    end
    volatility_technology_shocks=0.01*GDP_volatility_data./volatility_simulated_GDP;
end


for options=1:4
    chi=1;
    iferg=0;
    C_f_param=0;
    C_mean_param=0;
    eta=2;
    alpha=0.5;
    transferCoef=0;
    Z_volatility=0.01;
    compute_variance=0;

    save parameterfile alpha eta gamma chi varphi sigma iferg C_f_param C_mean_param transferCoef omega_W rho

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

        for i=1:length(MatrixParameters)
            close all
            set_param_value('gamma',  MatrixParameters(i,3));
            set_param_value('eta',  MatrixParameters(i,2));
            set_param_value('alpha',MatrixParameters(i,1));
            set_param_value('Z_volatility',volatility_technology_shocks(i,1));
            set_param_value('iferg',1)
    
            if options==1
                transferCoef0=(1-gamma)+(1-eta)*(1-alpha);
                optionsOpt = optimoptions('fminunc','Algorithm','quasi-newton','StepTolerance',1e-6,'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',200,'Display','final-detailed','FiniteDifferenceStepSize',1e-5);
                objective=@(x)autarky(x,V_Cf_0,i_0,C_mean_param_0, LABOR_SS, C_SS,sigma, chi, varphi, compute_variance);
                [transferOptimal(i,1),C_loss(i,1)]=fminunc(objective, transferCoef0,optionsOpt);
                [gamma eta alpha]
                C_loss(i,1)
            elseif options>1
                C_loss(i,1) =autarky(0,V_Cf_0,i_0,C_mean_param_0, LABOR_SS, C_SS,sigma, chi, varphi, 0);
            end
        end
    C_loss_cell{options}=C_loss;
end

gains_autarky=100*(C_loss_cell{2}-C_loss_cell{1});
gains_complete=100*(C_loss_cell{3}-C_loss_cell{1});
gains_incomplete=100*(C_loss_cell{4}-C_loss_cell{1});






toc
delete *.log *_dynamic.m *_static.m  *_set_auxiliary_variables.m *_results.mat *_results2.mat
delete  fu_dynare.m 
rmdir('fu_dynare', 's')

function f=autarky(transferCoef_loc,V_Cf_loc,i_loc,C_mean_param_loc, LABOR_SS, C_SS, sigma, chi, varphi,compute_variance )
        global oo_ M_ options_
        V_Cf_old=V_Cf_loc;
        i_old=i_loc;
        C_mean_param_old=C_mean_param_loc;
        set_param_value('transferCoef',transferCoef_loc);
        iter=0;
        set_param_value('iferg',1);
        eps=1;
        while eps>1e-7%eps>1e-8%eps>1e-7 %iter<5%eps>1e-7 
            set_param_value('C_f_param',V_Cf_old);
            set_param_value('i_star',i_old);
            set_param_value('C_mean_param',C_mean_param_old);
            %i_old
            %C_mean_param_old
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
        if compute_variance==0
            f =1- ((1-sigma).*( U + chi/(1+varphi).*exp(LABOR_SS).^(1+varphi))).^(1./(1-sigma))./exp(C_SS);
        elseif compute_variance==1
            f =sqrt(oo_.var(8,8));
        end
        set_param_value('C_f_param',V_Cf_loc)

       
end






