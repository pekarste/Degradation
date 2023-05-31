function [A,B] = A_B_reksten(separate_curve, k_2_0, k_4_0_plus, a_OH, T, E)
%A_B_reksten will be a helping function to "current_density_reksten" by
%calculating the values for A and B from the expression from Reksten

%% %%%%%%%%%%%%%%%%%% Description of input arguments %%%%%%%%%%%%%%%%%%%%%
%   k_2_0: the ration between the chemical backward and forward rate 
%          constant for step 2
%   k_4_0: the ration between the chemical backward and forward rate
%          constant for step 4 (which is the last step where oxygen is
%          produced (Given by thermodynamic constraint outside of function)
%   k_4_0_plus: The chemical forward rate constant for step 4
%   a_H_plus: The activity of H+, we assume it to be equal to the
%          concentration
%   T: Temperature in K
%   E: The potential
%
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
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
E_rev = E_n - ((R*T)/(n*F))*log((a_OH^2)/(a_H2O*sqrt(a_O2)));               % Reversible potential - with activity of O2
K = (sqrt(a_O2).*a_H2O)./a_OH^(2);% =1/(k_1_0*k_2_0*k_4_0)                  % Thermodynamic constraint on rate constants - perhaps this must be a separate function
%% %%%%%%%%%%%%%%%%% Importing from separate fit %%%%%%%%%%%%%%%%%%%% %%

k_2_0_plus = separate_curve.k_2_0_plus;                                     % k_2_0_plus from separate fit
k_1_0 = separate_curve.k_1_0;                                               % k_1_0 from separate fit
alpha = separate_curve.alpha;                                               % alpha from separate fit

%% Defining necessary functions


k_2_plus = k_2_0_plus.*exp((1-alpha).*F.*(E - E_rev)./(R*T));
k_2 = k_2_0.*exp(-F.*(E - E_rev)./(R*T));
k_4_0 = 1/(k_1_0*k_2_0*K);                                                  % Using the thermodynamic constraint on k_4_0

%% Defining A and B

A = ((1./k_4_0_plus) - (k_4_0.*sqrt(a_O2))./k_2_plus.*a_OH)./...
    (gamma.*(1 + k_4_0.*(k_2.*(a_H2O./a_OH) + 1).*sqrt(a_O2)));

B = (k_4_0.*sqrt(a_O2))./(1 + k_4_0.*(k_2.*(a_H2O./a_OH) + 1).*sqrt(a_O2));
end