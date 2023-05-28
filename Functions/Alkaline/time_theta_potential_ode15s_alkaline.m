function [t_ode15s_alkaline, theta_ode15s_alkaline,potential_ode15s_1 ,theta_interpol] = time_theta_potential_ode15s_alkaline(E_data, i_data, a_OH, T, data_type,k_4_0_plus)
%time_theta_potential_ode15s_alkaline will take in some a bunch of things
%and giove back arrays of t, theta, and interpolated values
%   Detailed explanation goes here

%% %%%%%%%%%%%%%%%%%%%%%% ALKALINE model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will couple all the different functions together and work
% like a masterscript. 


% Will try to use acidic equations but alkaline environemnt
%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_OER_SHE = 0.40;                                                           % Standard reduction potential for OER vs SHE - alkaline
E_REF_SHE = -0.829;                                                         % Standard redcution potential for HER vs SHE - alkaline
E_n = E_OER_SHE - E_REF_SHE;                                                % Standard reduction potential for OER vs RHE - alkaline
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
theta_2_0 = eps;


%% %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting the expression of the current based on r_2 to the data
% The r_2_fit returns the curve (fitting results) and the gof.
% The coefficients are contained in the curve

% Cherevko
[curve_alkaline, gof_alkaline] = ...                                        % This is the expression with rds
    r_2_fit_alkaline(E_data, i_data, a_OH, T, data_type);
%% %%%%%%%%%%%%%%%%% Data for degradation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Schalenbach_dissolution_CV_linear_data = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_dissolution_linear_alkaline.xlsx");
% Schalenbach dissolution vs time - [ng/cm^2*s]

Schalenbach_dissolution_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,2);    % Schalenbach dissolution data - [ng/cm^2*s]
Schalenbach_time_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,1);           % Schalenbach time data - [s]

Schalenbach_dissolution_mole = Schalenbach_dissolution_CV_linear*10^(-9)*10^(4)/Mm_Ir;  % Changes the units from ng/cm^2*s --> mole/m^2*s

%% %%%%%%%%%%%%%%%%%%%%%Solving differental equation%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
[t_ode15s_alkaline, theta_ode15s_alkaline] = diff_equation_solver_alkaline(Schalenbach_time_CV_linear, "value", curve_alkaline, a_OH, T, k_4_0_plus(1), theta_2_0);

%% %%%%%%%%%%Transforming time to potential for the ode15s solution
potential_ode15s_1 = CV_potential_alkaline(t_ode15s_alkaline, "array");

%% Interpolating the solution from ode15s to find values corresponding to
% the measured values since ode15s gives more points
theta_interpol = interp1(t_ode15s_alkaline,theta_ode15s_alkaline,Schalenbach_time_CV_linear); 

end