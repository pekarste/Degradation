function [k_1, k_1_plus, k_2, k_2_plus] = rate_constants_func_of_E(k_1_0, k_1_0_plus, k_2_0, k_2_0_plus, a_H_plus, alpha, T, E)
%rate_constants_func_of_E takes in the chemical forward rates and the
%ratio between the backward and forward for both electrochemical steps and
%gives out the rate constants with their potential dependence

%% %%%%%%%%%%%%%%%%%% Description of input arguments %%%%%%%%%%%%%%%%%%%%%
%   k_1_0: the ration between the chemical backward and forward rate 
%          constant for step 1 (Obtained from a separate fit)
%   k_2_0: the ration between the chemical backward and forward rate 
%          constant for step 2
%  
%   k_1_0_plus: The chemical forward rate constant for step 1
%   k_2_0_plus: The chemical forward rate constant for step 2 (obtained
%               from a separate fit)
%
%   a_H_plus: The activity of H+, we assume it to be equal to the
%          concentratiopn
%   alpha = the symmetry factor (obtained from a separate fit)
%   T: Temperature in K
%   E: The potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%% Defining constants %%%%%%%%%%%%%%%%%%%%%%%%% %%
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % [-] - number of electrons transferred in total
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-]
a_O2 = 0.21;                                                                % [-] - activity of oxygen

%% %%%%%%%%%%% Calculating the reversible potential %%%%%%%%%%%%%%%% %%
E_rev = E_n - ((R*T)/(n*F))*log(a_H2O/sqrt(a_O2)) + (R*T/F)*log(a_H_plus);  % [V] - Reversibel potential as a function of pH

%% %%%%%%%%%%% Defining the rate constants %%%%%%%%%%%%%%%%%% %%
k_1 = k_1_0.*exp(-F.*(E - E_rev)./(R*T));                                     % [-] - k_1_minus/k_1_plus
k_1_plus = k_1_0_plus.*exp((1 - alpha).*F.*(E - E_rev)./(R*T));                % [-] - forward rate constant for step 1

k_2 = k_2_0.*exp(-F.*(E - E_rev)./(R*T));                                     % [-] - k_2_minus/k_2_plus
k_2_plus = k_2_0_plus.*exp((1 - alpha).*F.*(E - E_rev)./(R*T));                % [-] - forward rate constant for step 2

end