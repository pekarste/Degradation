function current_density = current_density_reksten(separate_curve, k_1_0_plus, k_2_0, k_4_0_plus, a_H_plus, T, E)
%curent_density_reksten will calculate the current density based on the
%total expression from Reksten without any rds steps

%% %%%%%%%%%%%%%%%%%% Description of input arguments %%%%%%%%%%%%%%%%%%%%%
%   k_1_0: the ration between the chemical backward and forward rate 
%          constant for step 1 (Obtained from a separate fit)
%   k_2_0: the ration between the chemical backward and forward rate 
%          constant for step 2
%   k_4_0: the ration between the chemical backward and forward rate
%          constant for step 4 (which is the last step where oxygen is
%          produced)
%
%   k_1_0_plus: The chemical forward rate constant for step 1
%   k_2_0_plus: The chemical forward rate constant for step 2 (obtained
%               from a separate fit)
%   k_4_0_plus: The chemical forward rate constant for step 4
%
%   a_H_plus: The activity of H+, we assume it to be equal to the
%          concentratiopn
%   alpha = the symmetry factor (obtained from a separate fit)
%   T: Temperature in K
%   E: The potential
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%% Defining constants %%%%%%%%%%%%%%%%%%%%%%%%% %%
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % [-] - number of electrons transferred in total
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-]
a_O2 = 0.21;                                                                % [-] - activity of oxygen
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]


%% %%%%%%%%%%%%%%%%% Importing from separate fit %%%%%%%%%%%%%%%%%%%% %%

k_2_0_plus = separate_curve.k_2_0_plus;                                     % k_2_0_plus from separate fit
k_1_0 = separate_curve.k_1_0;                                               % k_1_0 from separate fit
alpha = separate_curve.alpha;                                               % alpha from separate fit


%% %%%% Defining the rate constants with their potential dependence %%%% %%

[k_1, k_1_plus, k_2, k_2_plus] = ...
    rate_constants_func_of_E(k_1_0, k_1_0_plus, k_2_0, k_2_0_plus, a_H_plus, alpha, T, E); 

%% %%%%%%%%%%%%%%%% Defining A and B %%%%%%%%%%%%%%%%%% %%

[A, B] = A_B_reksten(separate_curve, k_2_0, k_4_0_plus, a_H_plus, T, E);    % Calling A_B_reksten to give A and B
E_rev = E_n - ((R*T)/(n*F))*log(a_H2O/sqrt(a_O2)) + (R*T/F)*log(a_H_plus);  % [V] - Reversibel potential as a function of pH
K = (a_H_plus^(2)*sqrt(a_O2))/a_H2O;                                        % Thermodynamic constraint on rate constants - perhaps this must be a separate function
%% Calculating the current density

% current_density = n*F*gamma*(1 - B.*(1 + k_2.*a_H_plus + k_1.*k_2.*((a_H_plus)^2)./a_H2O))./...
%     (1./(k_1_plus.*a_H2O) + 1./(k_2_plus) + (k_1./k_2_plus).*(a_H_plus/a_H2O) + ... % Reksten has a minus here, but I have a plus...Now I try with a plus and see what happens
%     A.*gamma.*(1 + k_2.*a_H_plus + k_1.*k_2.*((a_H_plus)^2)/a_H2O));
current_density = n*F*gamma.* ...
        (1 - (((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(E - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2))).*...
        (1 + k_2_0.*exp(-F*(E - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(E - E_rev)/(R*T)).*k_2_0.*exp(-F*(E - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O))./ ...
        (1./((k_1_0_plus.*exp((1 - alpha)*F*(E - E_rev)/(R*T)))*a_H2O) + 1./(k_2_0_plus.*exp((1 - alpha)*F*(E - E_rev)/(R*T))) + (k_1_0.*exp(-F*(E - E_rev)/(R*T))./(k_2_0_plus.*exp((1 - alpha)*F*(E - E_rev)/(R*T))))*(a_H_plus/a_H2O) + ... 
        ((1/k_4_0_plus) - ((1/(k_1_0*k_2_0*K))*sqrt(a_O2))./(k_2_0_plus.*exp((1 - alpha)*F*(E - E_rev)/(R*T)))./...
        (gamma*(1 + (1/(k_1_0*k_2_0*K))*((k_2_0.*exp(-F*(E - E_rev)/(R*T))).*a_H_plus + 1)*sqrt(a_O2)))).*...
        gamma.*(1 + k_2_0.*exp(-F*(E - E_rev)/(R*T)).*a_H_plus + k_1_0.*exp(-F*(E - E_rev)/(R*T)).*k_2_0.*exp(-F*(E - E_rev)/(R*T)).*((a_H_plus)^2)/a_H2O));
end