function k_1 = k_1(k_1_0, a_OH, T, E)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% %%%%%%%%%%%%%%%% Defining constants %%%%%%%%%%%%%%%%%%%%%%%%% %%
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
n = 2;                                                                      % [-] - number of electrons transferred in total
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-]
a_O2 = 0.21;                                                                % [-] - activity of oxygen

%% %%%%%%%%%%% Calculating the reversible potential %%%%%%%%%%%%%%%% %%
E_rev = E_n - ((R*T)/(n*F))*log((a_OH^2)/(a_H2O*sqrt(a_O2)));               % Reversible potential - with activity of O2

k_1 = k_1_0.*exp(-F.*(E - E_rev)./(R.*T));

end