function d_gamma_theta_d_t = diff_equation_acidic(time,data_type,gamma_theta,curve,a_H,T, k_4_0_plus)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Define Physical Constants
%gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
%a_H2O = 1;
%% Calculating the potential from the time
E = CV_potential_acidic(time, data_type);                                   % Using the CV_potential to make the potential a function of time
% Based on Mayrhofer part 2
                                                                            % data_type makes it possible to use CV_potential with both ode15s and with an array of time
%% Extracting the coefficients from the separate fit

r_2_expression = r_2_acidic(curve, E, a_H, T);                              % Calling the function for the expression for the rate of the reaction
r_4_expression = k_4_0_plus*gamma_theta;                                    % Expressing r_4
d_gamma_theta_d_t = r_2_expression - r_4_expression;                        % Expressing the differential equation
end