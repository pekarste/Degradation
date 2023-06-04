function [t_ode15s_acidic, gamma_theta_ode15s_acidic,potential_ode15s_1 ,gamma_theta_interpol] = time_theta_potential_ode15s_acidic(E_data, i_data, a_H, T, data_type,k_4_0_plus)
%time_theta_potential_ode15s_acidic will take in some a bunch of things
%and giove back arrays of t, theta, and interpolated values
%   Detailed explanation goes here

%% %%%%%%%%%%%%%%%%%%%%%% Acidic model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will couple all the different functions together and work
% like a masterscript. 

%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_OER_SHE = 1.229;                                                          % Standard reduction potential for OER vs SHE - acidic
E_REF_RHE = 0.0;                                                            % Standard redcution potential for HER vs SHE - acidic
E_n = E_OER_SHE - E_REF_RHE;                                                % Standard reduction potential for OER vs RHE
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
%gamma = 8.16*10^(-6);                                                      % mol/m^2 [concentration of active sites]
Marhofer_a_H_acidic = 0.05*1;                                               % [-] - Activity of OH-
theta_2_0 = eps;        
%% %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting the expression of the current based on r_2 to the data
% The r_2_fit returns the curve (fitting results) and the gof.
% The coefficients are contained in the curve

% Cherevko
[curve_acidic, gof_acidic] = ...                                        % This is the expression with rds
    r_2_fit_acidic(E_data, i_data, a_H, T, data_type);
%% %%%%%%%%%%%%%%%%% Data for degradation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_2.xlsx");    % Mayrhofer dissolution vs time data - [ng/cm^2s]

%Mayrhofer_dissolution = Mayrhofer_dissolution_data(5:end,2);                % Mayrhofer dissolution data - [ng/cm^2*s] -- Starting from 5 to remove the tail
Mayrhofer_time = Mayrhofer_dissolution_data(5:end,1);                       % Mayrhofer time data [s] -- Starting from 5 to remove the tail to be consistent

%Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s

%% %%%%%%%%%%%%%%%%%%%%%Solving differental equation%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
[t_ode15s_acidic, gamma_theta_ode15s_acidic] = diff_equation_solver_acidic(Mayrhofer_time, "value", curve_acidic, Marhofer_a_H_acidic, T, k_4_0_plus, theta_2_0);

%% %%%%%%%%%%Transforming time to potential for the ode15s solution
potential_ode15s_1 = CV_potential_acidic(t_ode15s_acidic, "array");

%% Interpolating the solution from ode15s to find values corresponding to
% the measured values since ode15s gives more points
gamma_theta_interpol = interp1(t_ode15s_acidic,gamma_theta_ode15s_acidic,Mayrhofer_time); 

end