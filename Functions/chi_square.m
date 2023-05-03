function [curve, gof] = chi_square(theta_interpol, degradation_rate) % Should this function also take in the constants?
%chi_square wants to minimize the difference between the measured and calculated values of the dissolution of iridium
%   chi_square takes in the calculated rate of dissolution of irridium, and
%   compares it to the measured values from the Mayrhofer data. I don't
%   know yet if we want to run a minimalization of this based on one of the
%   coefficients or 

%% Define Physical Constants
a_H2O = 1;                                                                  % [-] - activity of water
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
Mm_Ir = 192.2;                                                              % g/mol [SI]
%% Degradation data
Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_2.xlsx");      % The dissolution data aquired from plot digitizer
Mayrhofer_dissolution = Mayrhofer_dissolution_data(1:end,2);                % Mayrhofer dissolution data - [ng/cm^2*s] 
Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s
%d_Ir_dt_meas = degradation_rate;                                  % Rate of dissolution from Mayrhofer

%% Fitting routine - minimize and give k_3_0_plus
fun = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*x;
FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', 10^(2), ...
           'StartPoint', 10^-6);                                              % k_3_0_plus
           


[curve, gof] = fit(theta_interpol,degradation_rate,FT,FO);                               % Curve contains the coefficients and gof some statistical data 

end