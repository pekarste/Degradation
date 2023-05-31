function [curve, gof] = r_rekstek_fit_no_curve(E_data,i_data, a_OH, Temperature, data_type)
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
E_OER_SHE = 0.40;                                                           % Standard reduction potential for OER vs SHE - alkaline
E_REF_RHE = -0.829;                                                         % Standard redcution potential for HER vs SHE - alkaline
E_n = E_OER_SHE - E_REF_RHE;                                                % Standard reduction potential for OER vs RHE - alkaline
a_H2O = 1;                                                                  % [-] - activity of water
a_O2 = 0.21;                                                                % [-] - activity of oxygen
%gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
T = Temperature;
E_rev = E_n - ((R*T)/(n*F))*log((a_OH^2)/(a_H2O*sqrt(a_O2)));               % Reversible potential - with activity of O2

K = (sqrt(a_O2).*a_H2O)./a_OH^(2);% =1/(k_1_0*k_2_0*k_4_0)                  % Thermodynamic constraint on rate constants - perhaps this must be a separate function


%% If data_type i linear
if data_type == "Linear"
    % Data
    %E = E_data;                                                             % [V] - Potential vs SHE/RHE
    %i = i_data;                                                             % [A/m^2] - Current density

    % Anonymous function
    fun = @(gamma, k_1_0, k_2_0, k_1_0_plus, alpha_1, k_2_0_plus, alpha_2, k_4_0_plus, x) n.*F*gamma.*a_OH.*(1 - B_reksten_no_curve(k_1_0, k_2_0, a_OH, T, x).*( 1 + k_2(k_2_0, a_OH, T, x).*(a_H2O./a_OH) + k_1(k_1_0, a_OH, T, x).*k_2(k_2_0, a_OH, T, x).*(a_H2O./(a_OH^(2))) ))./...
        (1./k_1_plus(k_1_0_plus, alpha_1, a_OH, T, x) + 1./k_2_plus(k_2_0_plus, alpha_2, a_OH, T, x) + k_1(k_1_0, a_OH, T, x)./(k_2_plus(k_2_0_plus, alpha_2, a_OH, T, x).*a_OH) + gamma.*A_reksten_no_curve(k_1_0, k_2_0, k_2_0_plus, k_4_0_plus, alpha_2, a_OH, T, x, gamma).*( 1 + k_2(k_2_0, a_OH, T, x).*(a_H2O./a_OH) + k_1(k_1_0, a_OH, T, x).*k_2(k_2_0, a_OH, T, x).*(a_H2O./(a_OH^(2))) )  );
    % Doing the fitting
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'gamma', 'k_1_0', 'k_2_0', 'k_1_0_plus','alpha_1','k_2_0_plus','alpha_2', 'k_4_0_plus'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps, eps, eps, eps, eps, eps],...                     % gamma, k_1_0, k_2_0, k_1_0_plus, alpha_1, k_2_0_plus, alpha_2, k_4_0_plus,
               'Upper',[10^(-2), 10^(10), 10^(10), 10^(10), 2, 10^(10), 2, 10^(10)],... % gamma, k_1_0, k_2_0, k_1_0_plus, alpha_1, k_2_0_plus, alpha_2, k_4_0_plus,
               'StartPoint',[10^(-5), 10^3, 10^1, 10^3, 0.5, 10^0, 0.5, 10^3],...       % Starting point for the coefficients
               'TolFun', 10^(-22));                                         % Tolerance    


    [curve, gof] = fit(E_data,i_data,FT,FO);                                          % Curve contains the coefficients and gof some statistical data 
%% If data_type is logarithmic
elseif data_type == "Logarithmic"
    %E = E_data;                                                             % [V] - Potential vs SHE/RHE
    %log_i = log10(i_data);                                                             % [A/m^2] - Current density

    % Anonymous function
    fun = @(gamma, k_1_0, k_2_0, k_1_0_plus, alpha_1, k_2_0_plus, alpha_2, k_4_0_plus, x) log10(n.*F*gamma.*a_OH) + log10((1 - B_reksten_no_curve(k_1_0, k_2_0, a_OH, T, x).*( 1 + k_2(k_2_0, a_OH, T, x).*(a_H2O./a_OH) + k_1(k_1_0, a_OH, T, x).*k_2(k_2_0, a_OH, T, x).*(a_H2O./(a_OH^(2))) ))) - ...
        log10((1./k_1_plus(k_1_0_plus, alpha_1, a_OH, T, x) + 1./k_2_plus(k_2_0_plus, alpha_2, a_OH, T, x) + k_1(k_1_0, a_OH, T, x)./(k_2_plus(k_2_0_plus, alpha_2, a_OH, T, x).*a_OH) + gamma.*A_reksten_no_curve(k_1_0, k_2_0, k_2_0_plus, k_4_0_plus, alpha_2, a_OH, T, x, gamma).*( 1 + k_2(k_2_0, a_OH, T, x).*(a_H2O./a_OH) + k_1(k_1_0, a_OH, T, x).*k_2(k_2_0, a_OH, T, x).*(a_H2O./(a_OH^(2))) )  ));
    % Doing the fitting
    FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'gamma', 'k_1_0', 'k_2_0', 'k_1_0_plus','alpha_1','k_2_0_plus','alpha_2', 'k_4_0_plus'});

    FO = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[eps, eps, eps, eps, eps, eps, eps, eps],...                     % gamma, k_1_0, k_2_0, k_1_0_plus, alpha_1, k_2_0_plus, alpha_2, k_4_0_plus,
               'Upper',[10^(-2), 10^(10), 10^(10), 10^(10), 2, 10^(10), 2, 10^(10)],... % gamma, k_1_0, k_2_0, k_1_0_plus, alpha_1, k_2_0_plus, alpha_2, k_4_0_plus,
               'StartPoint',[10^(-5), 10^5, 10^3, 10^2, 0.5, 10^0, 0.5, 10^2],...       % Starting point for the coefficients
               'TolFun', 10^(-22));                                         % Tolerance    


    [curve, gof] = fit(E_data,log10(i_data),FT,FO);                                          % Curve contains the coefficients and gof some statistical data 
end

