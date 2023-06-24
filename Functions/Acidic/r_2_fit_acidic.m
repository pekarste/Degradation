function [curve, gof] = r_2_fit_acidic(E_data,i_data, a_H_plus, T, data_type)
%r_2_fit Fiting the expressionn for r_2 to numerical data
%   r_2_fit takes in potential and curretn density and uses it to produce
%   coefficient for the expression of current density derived from kinetics
%   where the second step is rds

%% Define Physical Constants
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % Number of electrons transferred
E_OER_SHE = 1.229;                                                          % Standard reduction potential for OER vs SHE - acidic
E_REF_RHE = 0.0;                                                            % Standard redcution potential for HER vs SHE - acidic
E_n = E_OER_SHE - E_REF_RHE;                                                % Standard reduction potential for OER vs RHE
a_H2O = 1;                                                                  % [-] - activity of water
a_O2 = 0.21;                                                                % [-] - activity of oxygen
%gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

E_rev = E_n - ((R.*T)./(n.*F)).*log(a_H2O./sqrt(a_O2)) + ((R.*T)./F).*log(a_H_plus);% Reversible potential - with activity of O2


% Expression that is to be fitted is:
% i = 2*F*r_2 = 2*F*(gamma*k_20p*exp((1-alpha)*F*(E-E_n)/(R*T)))/(1 + (a_Hplus/a_H2O)*k_1_0*exp((-1)*F*(E-E_n)/(R*T)));
% log(i) = log(2*F*gamma_k_2_0_plus) + (1-alpha)*F*(E-E_n)/(R*T) - log(1 +
% Hplus/a_H2O)*k_1_0*exp((-1)*F*(E-E_n)/(R*T)))

%% If the data is linear then do this:

if data_type == "Linear"
    % Data
    E = E_data;                                                             % [V] - Potential vs SHE/RHE
    i = i_data;                                                             % [A/m^2] - Current density
    
    fun = @(gamma_k_2_0_plus, alpha, k_1_0, x) n.*F.*gamma_k_2_0_plus.*exp((1-alpha).*F.*(x - E_rev)/(R.*T))./(1 + (a_H_plus./a_H2O).*k_1_0.*exp(-F.*(x - E_rev)./(R.*T)));
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'gamma_k_2_0_plus', 'alpha', 'k_1_0'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps],...                                  % gamma_k_2_0_plus alpha k_1_0
               'Upper',[10^(1), 2, 10^(8)],...                              % gamma_k_2_0_plus alpha k_1_0
               'StartPoint',[5.182*10^(-5), 0.5206, 2.734*10^(4)],...       % Starting point for the coefficients
               'TolFun', 1e-20);                                            % Tolerance         

    [curve, gof] = fit(E,i,FT,FO);                                          % Curve contains the coefficients and gof some statistical data 
    
%% If data_type is logaritmic, then do this:

elseif data_type == "Logarithmic"
    % Data
    E = E_data;                                                             % [V] - Potential vs SHE/RHE
    log_i = log10(i_data);                                                  % [log(A/m^2)] - Current density

    fun = @(gamma_k_2_0_plus,alpha, k_1_0,x) log10(n.*F.*gamma_k_2_0_plus) + (1-alpha).*F.*(x - E_rev).*log10(exp(1))./(R.*T) - log10((1 + (a_H_plus./a_H2O).*k_1_0.*exp(-F.*(x - E_rev)/(R.*T))));
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'gamma_k_2_0_plus', 'alpha', 'k_1_0'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps],...                                  % gamma_k_2_0_plus alpha k_1_0
               'Upper',[10^(1), 2, 10^(8)],...                              % gamma_k_2_0_plus alpha k_1_0
               'StartPoint',[5.182*10^(-5), 0.5206, 2.734*10^(4)],...       % Starting point for the coefficients
               'TolFun', 1e-20);                                            % Tolerance   


    [curve, gof] = fit(E,log_i,FT,FO);                                      % Curve contains the coefficients and gof some statistical data 
    
end