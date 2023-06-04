function [curve] = chi_square_alkaline(time_solver, gamma_theta_solver, a_OH, k_4_0_plus, time_data, degradation_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%gamma = 8.16*10^(-6);                                                       % [mol/m^^2] 

fun = @(k_3_0_plus, x) k_3_0_plus.*a_OH.^(2).*interp1(time_solver,gamma_theta_solver, x);
FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', 10^4, ...
           'StartPoint', k_4_0_plus,...
           'TolFun', 1e-20);                                                % k_3_0_plus

[curve, gof, output,warnstr,errstr,convmsg]...
    = fit(time_data,degradation_data,FT,FO);    

end