function d_theta_d_t = diff_equation_alkaline(time,data_type,theta,curve,a_OH,T, k_4_0_plus)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Define Physical Constants
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

%% Calculating the potential from the time
E = CV_potential_alkaline(time, data_type);                                 % Using the CV_potential to make the potential a function of time
% Based on Schalenbach_dissolution_CV_linear
                                                                            % data_type makes it possible to use CV_potential with both ode15s and with an array of time
%% Extracting the coefficients from the separate fit

r_2_expression = r_2_alkaline(curve, E, a_OH, T);                           % Calling the function for the expression for the rate of the reaction

d_theta_d_t = r_2_expression/gamma - k_4_0_plus*theta;                      % Divinding by gamma since r_2_expressions contains gamma from before, eliminating gamma from the equation

end