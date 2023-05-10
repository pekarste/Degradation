function [k_3_0_plus_scohy, k_3_0_plus_damjanovic, k_3_0_plus_damjanovic_log] = k_3_0_plus_extraction(k_4_0_plus)
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

% Data used for fitting of r2
Scohy_Ir_data = readmatrix("Scohy_activated_Ir_LSV.xlsx");                  % Potential/current density data from Scohy
Damjanovic_Ir_data = readmatrix("Damjanovic_Ir_E_vs_log_i.xlsx");           % Current censity/potential data from Damjanovic

% Extracted data from the Excel files 

% Scohy
Scohy_potential = Scohy_Ir_data(1:end,1);                                   % [V vs RHE] - Potential
Scohy_current_density = Scohy_Ir_data(1:end,2)*10^(-3+4);                   % [A/m^2] - Current density, originally in mA/cm^2
Scohy_T = 25 + 273;                                                         % [K] - Temperature
Scohy_a_H_plus = 0.5*2;                                                     % [-] - Activity of H+


%Damjanovic
Damjanovic_potential = Damjanovic_Ir_data(1:end,2);                         % [V vs RHE]          
Damjanovic_current_density = Damjanovic_Ir_data(1:end,1)*10^4;              % [A/m^2] - Originally A/cm^2 (The article contains log(i), but I converted thm
Damjanovic_T = 25 + 273;                                                    % [K] - Temperature
Damjanovic_a_H_plus = 1;                                                    % [-] - Activity of H+


%% %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting the expression of the current based on r_2 to the data
% The r_2_fit returns the curve (fitting results) and the gof.
% The coefficients are contained in the curve

% Scohy
[Scohy_curve, Scohy_gof] = ...                                              % This is the expression with rds
    r_2_fit(Scohy_potential, Scohy_current_density, Scohy_a_H_plus, Scohy_T, "Linear");

%[Scohy_curve_fit, Scohy_gof_fit] = ...                                      % This is the expression with no rds, but k_4 << 1
%    r_fit(Scohy_potential, Scohy_current_density, Scohy_a_H_plus, Scohy_T, "Linear");

% Damjanovic
[Damjanovic_curve, Damjanovic_gof] = ...                                    % This is the expression with rds
    r_2_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Linear");


%[Damjanovic_curve_fit, Damjanovic_gof_fit] = ...                            % This is the expression with no rds, but k_4 << 1
%    r_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Linear");

% Damjanovic (log)
[Damjanovic_log_curve, Damjanovic_log_gof] = ...                            % This is the expression with rds
    r_2_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Logarithmic");

%[Damjanovic_log_curve_fit, Damjanovic_log_gof_fit] = ...                    % This is the expression with no rds, but k_4 << 1
%    r_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Logarithmic");

%% %%%%%%%%%%% The data from the Mayrhofer article %%%%%%%%%%%%%%%%%%%%%%
% These data is based on the highest anodic peak

Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_2.xlsx");    % Mayrhofer dissolution vs time data - [ng/cm^2s]

Mayrhofer_dissolution = Mayrhofer_dissolution_data(5:end,2);                % Mayrhofer dissolution data - [ng/cm^2*s] -- Starting from 5 to remove the tail
Mayrhofer_time = Mayrhofer_dissolution_data(5:end,1);                       % Mayrhofer time data [s] -- Starting from 5 to remove the tail to be consistent

Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s

Mayrhofer_a_H_plus = 0.1*2;                                                 % Concentration of H+ (0.1 M H2SO4)
Mayrhofer_T = 25 + 273.13;                                                  % mayrhofer states room temperature

%% %%%%%%%%%%%%%%% Calling the diff equation solver %%%%%%%%%%%%%%%%%%%%%%

[t_scohy, theta_scohy] = diff_equation_solver(Mayrhofer_time, "value", Scohy_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj, theta_damj] = diff_equation_solver(Mayrhofer_time, "value", Damjanovic_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj_log, theta_damj_log] = diff_equation_solver(Mayrhofer_time, "value", Damjanovic_log_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);


%% Extracting the k_3_0_plus value  
% Normalisation just wouldn't do

fun_scohy = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_scohy,theta_scohy,x);
fun_damjanovic = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_damj,theta_damj,x);
fun_damjanovic_log = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_damj_log,theta_damj_log,x);


FT_scohy = fittype(fun_scohy, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});
FT_damjanovic = fittype(fun_damjanovic, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});
FT_damjanovic_log = fittype(fun_damjanovic_log, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', k_4_0_plus, ...
           'StartPoint', 1e-4,...
           'TolFun', 1e-22);                                                % k_3_0_plus
          

[curve_scohy, gof, output,warnstr,errstr,convmsg] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT_scohy,FO);
[curve_damjanovic, gof, output,warnstr,errstr,convmsg] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT_damjanovic,FO);
[curve_damjanovic_log, gof, output,warnstr,errstr,convmsg] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT_damjanovic_log,FO);

k_3_0_plus_scohy = curve_scohy.k_3_0_plus;
k_3_0_plus_damjanovic = curve_damjanovic.k_3_0_plus;
k_3_0_plus_damjanovic_log = curve_damjanovic_log.k_3_0_plus;

end