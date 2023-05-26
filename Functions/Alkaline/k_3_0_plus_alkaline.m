function [curve, gof] = k_3_0_plus_alkaline(a_OH,t_ode15s, theta_ode15s, time, dissolution, k_4_0_plus)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Define Physical Constants
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % Number of electrons transferred
E_n = 1.229;                                                                 % V - If it is set to 0.4 the log curve won't fit...?
%E_n = 0.399; %It is this for the reaction, but I don't know why this gives
% bad graphs
a_H2O = 1;                                                                  % [-] - activity of water
a_O2 = 0.21;                                                                % [-] - activity of oxygen
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

%% minimisation
fun = @(k_3_0_plus, x) gamma*k_3_0_plus*a_OH^(2)*interp1(t_ode15s,theta_ode15s,x);
FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', k_4_0_plus, ...
           'StartPoint', 1e-4,...
           'TolFun', 1e-22);                                                % k_3_0_plus
          

[curve, gof, output,warnstr,errstr,convmsg] = fit(time,dissolution,FT,FO);
end