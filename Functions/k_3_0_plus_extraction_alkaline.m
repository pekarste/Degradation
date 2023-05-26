function [k_3_0_plus_cherevko,k_3_0_plus_damjanovic, k_3_0_plus_damjanovic_log, k_3_0_plus_schalenbach] = k_3_0_plus_extraction_alkaline(k_4_0_plus)
% k_3_0_plus_extraction takes in a value for k_4_0_plus, and calculates a
% value for k_3_0_plus, based on this value. It gives back a value based on
% scohy and damjanovic. It is based on the "Degradation_model" script


%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Data used for fitting r2_alkaline
Cherevko_alkaline = readmatrix("Data\Alkaline\Cherevko_alkaline_polarisation_data.xlsx");% Potential/current density data from Cherevko
Damjanovic_alkaline = readmatrix("Data\Alkaline\Damjanovic_alkaline_polarisation.xlsx"); % Current density/potential data from Damjanovic
%--------------------------------------------------------------------------
Schalenbach_polarisation_alkaline = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_polarisation_curve_alkaline.xlsx");% Polarisation curve from alkaline dissolution
%--------------------------------------------------------------------------
%% Extracted data from the Excel files 

% Alkaline - Cherevko
Cherevko_E_alkaline = Cherevko_alkaline(1:end,1);                           % [V vs RHE] - Potential
Cherevko_i_alkaline = Cherevko_alkaline(1:end,2)*10^(-3+4);                 % [A/m^2] - Current density, originally in mA/cm^2
Cherevko_T_alkaline = 25 + 273;                                             % [K] - Temperature
Cherevko_OH_alkaline = 0.05*2;                                              % [-] - Activity of OH-


% Alkaline - Damjanovic
Damjanovic_E_alkaline = Damjanovic_alkaline(1:end,2);                       % [V vs RHE]          
Damjanovic_i_alkaline = Damjanovic_alkaline(1:end,1)*10^4;                  % [A/m^2] - Originally A/cm^2 (The article contains log(i), but I converted thm
Damjanovic_T_alkaline = 25 + 273;                                           % [K] - Temperature
Damjanovic_OH_alkaline = 1;                                                 % [-] - Activity of OH-
%_-------------------------------------------------------------------------
% Alkaline - Schalenbach
Schalenbach_E_alkaline = Schalenbach_polarisation_alkaline(1:end,1);           % [V vs RHE] - Potential
Schalenbach_i_alkaline = Schalenbach_polarisation_alkaline(1:end,2)*10^(-3+4); % [A/m^2] - Current density, originally in mA/cm^2
Schalenbach_T_alkaline = 25 + 273;                                             % [K] - Temperature
Schalenbach_a_OH_alkaline = 0.05*1;                                              % [-] - Activity of OH-
Schalenbach_sweep_rate = 2*10^(-3);                                            % [V/s] -Schalenbach Sweep rate 
%--------------------------------------------------------------------------

%% Fitting
% Cherevko
[Cherevko_curve_alkaline, Cherevko_gof_alkaline] = ...                                        % This is the expression with rds
    r_2_fit_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear");

% Damjanovic
[Damjanovic_curve_alkaline, Damjanovic_gof_alkaline] = ...                                    % This is the expression with rds
    r_2_fit_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear");

% Damjanovic (log)
[Damjanovic_log_curve_alkaline, Damjanovic_log_gof_alkaline] = ...                            % This is the expression with rds
    r_2_fit_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic");

%--------------------------------------------------------------------------
% Schalenbach - Alkaline 
[Schalenbach_curve_alkaline, Schalenbach_gof_alkaline] = ...                            % This is the expression with rds
    r_2_fit_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_a_OH_alkaline, Schalenbach_T_alkaline, "Linear");
%--------------------------------------------------------------------------

%% %%%%%%%%%%% The data from the Schalenbach article %%%%%%%%%%%%%%%%%%%%%%
% These data is based on the highest anodic peak

Schalenbach_dissolution_CV_linear_data = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_dissolution_linear_alkaline.xlsx");
% Schalenbach dissolution vs time - [ng/cm^2*s]
    
Schalenbach_dissolution_CV_semilog_data = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_semilog_dissolution_alkaline.xlsx");
% Schalenbach dissolution vs time - [ng/cm^2*s]

Schalenbach_dissolution_v_potental_data = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_dissolution_vs_potential_polarisation_alkaline.xlsx");
% Schalenbach dissolution vs potential - [ng/cm^2*s] --- check if constant sweep rate


Schalenbach_dissolution_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,2);    % Schalenbach dissolution data - [ng/cm^2*s]
Schalenbach_time_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,1);           % Schalenbach time data - [s]

Schalenbach_dissolution_mole = Schalenbach_dissolution_CV_linear*10^(-9)*10^(4)/Mm_Ir;  % Changes the units from ng/cm^2*s --> mole/m^2*s

%% %%%%%%%%%%%%%%% Calling the diff equation solver %%%%%%%%%%%%%%%%%%%%%%

[t_cherevko_alkaline, theta_cherevko_alkaline] = diff_equation_solver_alkaline(Schalenbach_time_CV_linear, "value", Cherevko_curve_alkaline, Schalenbach_a_OH_alkaline, Schalenbach_T_alkaline, k_4_0_plus, theta_2_0);
[t_damj_alkaline, theta_damj_alkaline] = diff_equation_solver_alkaline(Schalenbach_time_CV_linear, "value", Damjanovic_curve_alkaline, Schalenbach_a_OH_alkaline, Schalenbach_T_alkaline, k_4_0_plus, theta_2_0);
[t_damj_log_alkaline, theta_damj_log_alkaline] = diff_equation_solver_alkaline(Schalenbach_time_CV_linear, "value", Damjanovic_log_curve_alkaline, Schalenbach_a_OH_alkaline, Schalenbach_T_alkaline, k_4_0_plus, theta_2_0);
%--------------------------------------------------------------------------
[t_schalenbach_alkaline, theta_schalenbach_alkaline] = diff_equation_solver_alkaline(Schalenbach_time_CV_linear, "value", Schalenbach_curve_alkaline, Schalenbach_a_OH_alkaline, Schalenbach_T_alkaline, k_4_0_plus, theta_2_0);
%--------------------------------------------------------------------------

%% Extracting the k_3_0_plus value  
% Normalisation just wouldn't do

fun_cherevko = @(k_3_0_plus, x) gamma*k_3_0_plus*Schalenbach_a_OH_alkaline*interp1(t_cherevko_alkaline,theta_cherevko_alkaline,x);
fun_damjanovic = @(k_3_0_plus, x) gamma*k_3_0_plus*Schalenbach_a_OH_alkaline*interp1(t_damj_alkaline,theta_damj_alkaline,x);
fun_damjanovic_log = @(k_3_0_plus, x) gamma*k_3_0_plus*Schalenbach_a_OH_alkaline*interp1(t_damj_log_alkaline,theta_damj_log_alkaline,x);
fun_schalenbach = @(k_3_0_plus, x) gamma*k_3_0_plus*Schalenbach_a_OH_alkaline*interp1(t_schalenbach_alkaline,theta_schalenbach_alkaline,x);


FT_cherevko = fittype(fun_cherevko, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});
FT_damjanovic = fittype(fun_damjanovic, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});
FT_damjanovic_log = fittype(fun_damjanovic_log, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});
FT_schalenbach = fittype(fun_schalenbach, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', k_4_0_plus, ...
           'StartPoint', k_4_0_plus*10^(-2),...
           'TolFun', 1e-22);                                                % k_3_0_plus
          

[curve_cherevko, gof, output,warnstr,errstr,convmsg] = fit(Schalenbach_time_CV_linear,Schalenbach_dissolution_mole,FT_cherevko,FO);
[curve_damjanovic, gof, output,warnstr,errstr,convmsg] = fit(Schalenbach_time_CV_linear,Schalenbach_dissolution_mole,FT_damjanovic,FO);
[curve_damjanovic_log, gof, output,warnstr,errstr,convmsg] = fit(Schalenbach_time_CV_linear,Schalenbach_dissolution_mole,FT_damjanovic_log,FO);
[curve_schalenbach, gof, output,warnstr,errstr,convmsg] = fit(Schalenbach_time_CV_linear,Schalenbach_dissolution_mole,FT_schalenbach,FO);

k_3_0_plus_cherevko = curve_cherevko.k_3_0_plus;
k_3_0_plus_damjanovic = curve_damjanovic.k_3_0_plus;
k_3_0_plus_damjanovic_log = curve_damjanovic_log.k_3_0_plus;
k_3_0_plus_schalenbach = curve_schalenbach.k_3_0_plus;


end