function [curve] = chi_square_acidic(time_solver, gamma_theta_solver, k_4_0_plus, time_data, degradation_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%gamma = 8.16*10^(-6);                                                       % [mol/m^^2] 
a_H2O = 1;

fun = @(k_3_0_plus, x) k_3_0_plus.*a_H2O.*interp1(time_solver,gamma_theta_solver, x);
FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', 10^4, ...
           'StartPoint', k_4_0_plus,...
           'TolFun', 1e-20);                                                % k_3_0_plus

[curve, gof, output,warnstr,errstr,convmsg]...
    = fit(time_data,degradation_data,FT,FO);    

end