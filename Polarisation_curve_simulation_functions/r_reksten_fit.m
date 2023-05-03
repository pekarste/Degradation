function [curve, gof] = r_reksten_fit(separate_curve, E_data,i_data, a_H_plus, T, data_type)
%r_reksten_fit will fit an expression based on the work of Reksten for the
%total current without assuming any rds steps and using steady state.

%% %%%%%%%%%%%%%%%%%% Description of input arguments %%%%%%%%%%%%%%%%%%%%%
%   curve: This will either be the prevoius fit from Scohy or Damjanovic
%          based ont the first fitting. Can base it on the second fitting
%          (which is an approximation of this fit, but k_1_0_plus has so
%          much uncertainty)
%          Will contain: k_2_0_plus, alpha, k_1_0
%
%   E_data: The potential data from Scohy or Damjanovic
%  
%   i_data: The current density data from Scohy or Damjanovic
%
%   a_H_plus: The activity of H+, we assume it to be equal to the
%          concentratiopn
%  
%   T: Temperature in K
%   data_type: This is just so we can use Damajnovic logarithmic data since
%   that will give a better fit with the Damjanovic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define Physical Constants
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % Number of electrons transferred
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-] - activity of water
a_O2 = 0.21;                                                                % [-] - activity of oxygen
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

E_rev = E_n - ((R*T)/(n*F))*log(a_H2O/sqrt(a_O2)) + ((R*T)/F)*log(a_H_plus);% Reversible potential

K = (a_H_plus^(2)*sqrt(a_O2))/a_H2O;                                        % Thermodynamic constraint

%% %%%%%%%%%%%%%%%%% Importing from separate fit %%%%%%%%%%%%%%%%%%%% %%

k_2_0_plus = separate_curve.k_2_0_plus;                                     % k_2_0_plus from separate fit
k_1_0 = separate_curve.k_1_0;                                               % k_1_0 from separate fit
alpha = separate_curve.alpha;                                               % alpha from separate fit


%% If the data is linear then do this:

if data_type == "Linear"
    % Data
    E = E_data;                                                             % [V] - Potential vs SHE/RHE
    i = i_data;                                                             % [A/m^2] - Current density
    
    % Anonymous function
%     fun = @(k_2_0, k_1_0_plus, k_4_0_plus, x) n*F*gamma*(1 - (((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2))).*(1 + (k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + (k_1_0.*k_2_0.*exp(-2*F*(x - E_rev)/(R*T))).*((a_H_plus)^2)/a_H2O))./...
%     (1./((k_1_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T)))*a_H2O) + 1./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))) + ((k_1_0.*exp(-F*(x - E_rev)/(R*T)))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))))*(a_H_plus/a_H2O) - ... % Reksten has a minus here, but I have a plus...Now I try with a plus and see what happens
%     (((1/k_4_0_plus) - ((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))))./...
%     (gamma*(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2))))*gamma.*(1 + (k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + (k_1_0.*k_2_0.*exp(-2*F*(x - E_rev)/(R*T))).*((a_H_plus)^2)/a_H2O));

    fun = @(k_2_0, k_1_0_plus, k_4_0_plus, x) n*F*gamma.* ...
        (1 - (((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2))).*...
        (1 + k_2_0.*exp(-F*(x - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(x - E_rev)/(R*T)).*k_2_0.*exp(-F*(x - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O))./ ...
        (1./((k_1_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T)))*a_H2O) + 1./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))) + (k_1_0.*exp(-F*(x - E_rev)/(R*T))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))))*(a_H_plus/a_H2O) + ... 
        ((1/k_4_0_plus) - ((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T)))./...
        (gamma*(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2)))).*...
        gamma.*(1 + k_2_0.*exp(-F*(x - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(x - E_rev)/(R*T)).*k_2_0.*exp(-F*(x - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O));
    %% Doing the fitting
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_2_0', 'k_1_0_plus', 'k_4_0_plus'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps],...                                  % k_2_0, k_1_0_plus, k_4_0_plus
               'Upper',[10^(10), 10^(10), 10^(10)],...                      % k_2_0, k_1_0_plus, k_4_0_plus
               'StartPoint',[10^2, 10^3, 10^3]);                            % Starting point for the coefficients


    [curve, gof] = fit(E,i,FT,FO);                                          % Curve contains the coefficients and gof some statistical data 

%% If data_type is logaritmic, then do this

elseif data_type == "Logarithmic"
    % Data
    E = E_data;                                                             % [V] - Potential vs SHE/RHE
    log_i = log10(i_data);                                                  % [log(A/m^2)] - Current density

    fun = @(k_2_0, k_1_0_plus, k_4_0_plus, x) log10(n*F*gamma) + ...
        log10((1 - (((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2))).*...
        (1 + k_2_0.*exp(-F*(x - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(x - E_rev)/(R*T)).*k_2_0.*exp(-F*(x - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O))) - ...
        log10((1./((k_1_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T)))*a_H2O) + 1./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))) + (k_1_0.*exp(-F*(x - E_rev)/(R*T))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))))*(a_H_plus/a_H2O) + ... 
        ((1/k_4_0_plus) - ((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T)))./...
        (gamma*(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2)))).*...
        gamma.*(1 + k_2_0.*exp(-F*(x - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(x - E_rev)/(R*T)).*k_2_0.*exp(-F*(x - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O)));
    
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_2_0', 'k_1_0_plus', 'k_4_0_plus'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps],...                                  % k_4_0_plus, k_2_0, k_1_0_plus
               'Upper',[10^(10), 10^(10), 10^(10)],...                      % k_4_0_plus, k_2_0, k_1_0_plus
               'StartPoint',[10^3, 10^1, 10^3]);                            % Starting point for the coefficients

    [curve, gof] = fit(E,log_i,FT,FO);                                      % Curve contains the coefficients and gof some statistical data 

end

% @(k_2_0,k_1_0_plus, k_4_0_plus,x) log10(n*F*gamma) + ...
%         log10((1 - ((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(1 + (1/(k_1_0*k_2_0*K))*((k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2)).*...
%         (1 + k_2_0.*exp(-F*(x - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(x - E_rev)/(R*T)).*k_2_0.*exp(-F*(x - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O))) - ...
%         log10((1./((k_1_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T)))*a_H2O) + 1./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))) + (k_1_0.*exp(-F*(x - E_rev)/(R*T))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))))*(a_H_plus/a_H2O) - ... 
%         ((1/k_4_0_plus) - ((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T)))./...
%         (gamma*(1 + (1/(k_1_0*k_2_0*K))*((k_2_0_plus.*exp((1 - alpha)*F*(x - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2)))).*...
%         gamma.*(1 + k_2_0.*exp(-F*(x - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(x - E_rev)/(R*T)).*k_2_0.*exp(-F*(x - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O))),..