function r_2_expression = r_2(curve,E, a_H_plus, T)
% r_2 takes in the coefficient and potential and returns the expression for
% the rate of the OER reaction 

%% Define Physical Constants
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % [-]
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-]
a_O2 = 0.21;                                                                % [-]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

E_rev = E_n - ((R*T)/(n*F))*log(a_H2O/sqrt(a_O2)) + ((R*T)/F)*log(a_H_plus);
%% Expressing the rate of step 2
k_2_0_plus = curve.k_2_0_plus;                                              % Extracting k_2_0_plus from the separate fit
k_1_0 = curve.k_1_0;                                                        % Extracting k_1_0 from the separate fit
alpha = curve.alpha;                                                        % Extracting talpha from the separate fit

r_2_expression = gamma*k_2_0_plus*exp((1-alpha)*F*(E - E_rev)/(R*T))./...
    (1 + (a_H_plus/a_H2O)*k_1_0*exp(-F*(E - E_rev)/(R*T)));
end