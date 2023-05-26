function [curve, gof] = r_2_fit_alkaline(E_data,i_data, a_OH, T, data_type)
%r_2_fit_alkaline Fiting the expressionn for r_2_alkaline to numerical data
%   r_2_fit_alkaline takes in potential and current density and uses it to produce
%   coefficient for the expression of current density derived from kinetics
%   where the second step is rds

%% Define Physical Constants
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % Number of electrons transferred
E_OER_SHE = 0.40;                                                           % Standard reduction potential for OER vs SHE - alkaline
E_REF_RHE = -0.829;                                                         % Standard redcution potential for HER vs SHE - alkaline
E_n = E_OER_SHE - E_REF_RHE;                                                % Standard reduction potential for OER vs RHE - alkaline
a_H2O = 1;                                                                  % [-] - activity of water
a_O2 = 0.21;                                                                % [-] - activity of oxygen
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

E_rev = E_n - ((R*T)/(n*F))*log((a_OH^2)/(a_H2O*sqrt(a_O2)));               % Reversible potential - with activity of O2


% Expression that is to be fitted is:
% i = 2*F*r_2 = 2*F*(gamma*k_20p*exp((1-alpha)*F*(E-E_n)/(R*T)))*a_OH/(1 + (k_1_0/a_OH)*exp((-1)*F*(E-E_n)/(R*T)));
% log(i) = log(2*F*gamma_k_2_0_plus*a_OH) + (1-alpha)*F*(E-E_n)/(R*T)*log(e) - log(1 +
% (k_1_0/a_OH)*exp((-1)*F*(E-E_n)/(R*T)))

%% If the data is linear then do this:

if data_type == "Linear"
    % Data
    E = E_data;                                                             % [V] - Potential vs SHE/RHE
    i = i_data;                                                             % [A/m^2] - Current density
    
    fun = @(k_2_0_plus, alpha, k_1_0, x) n*F*gamma*k_2_0_plus*exp((1-alpha)*F*(x - E_rev)/(R*T))*a_OH./(1 + (k_1_0/a_OH)*exp(-F*(x - E_rev)/(R*T)));
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_2_0_plus', 'alpha', 'k_1_0'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps],...                                  % k_2_0_plus alpha k_1_0
               'Upper',[10^(8), 2, 10^(8)],...                              % k_2_0_plus alpha k_1_0
               'StartPoint',[1*10^(-2), 0.5, 1*10^(4)]);         % Starting point for the coefficients


    [curve, gof] = fit(E,i,FT,FO);                                          % Curve contains the coefficients and gof some statistical data 
    
%% If data_type is logaritmic, then do this:

elseif data_type == "Logarithmic"
    % Data
    E = E_data;                                                             % [V] - Potential vs SHE/RHE
    log_i = log10(i_data);                                                  % [log(A/m^2)] - Current density

    fun = @(k_2_0_plus,alpha, k_1_0,x) log10(n*F*gamma*k_2_0_plus*a_OH) + (1-alpha)*F*(x - E_rev)*log10(exp(1))/(R*T) - log10((1 + (k_1_0/a_OH)*exp(-F*(x - E_rev)/(R*T))));
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_2_0_plus', 'alpha', 'k_1_0'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps],...                                  % k_2_0_plus alpha k_1_0
               'Upper',[10^2, 2, 10^(8)],...                                % k_2_0_plus alpha k_1_0
               'StartPoint',[1*10^(-2), 0.5, 1*10^(4)]);         % Starting point for the coefficients


    [curve, gof] = fit(E,log_i,FT,FO);                                      % Curve contains the coefficients and gof some statistical data 
    
end