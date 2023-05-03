% This script will be used to try and simulate the polarisation curves
% based on the expression from Reksten and perhaps the one from Marshall as
% well. This is becuase we wont do another fitting, but rather find some
% values for k_4 which are within reasonable values


%% Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-] - activity of water
a_O2 = 0.21;                                                                % [-] - activity of oxygen
Mm_Ir = 192.2;                                                              % g/mol [SI]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
n = 2;                                                                      % [-] - number of electrons being transferred


%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Data used for fitting of r2
Scohy_Ir_data = readmatrix("Scohy_activated_Ir_LSV.xlsx");                  % Potential/current density data from Scohy
Damjanovic_Ir_data = readmatrix("Damjanovic_Ir_E_vs_log_i.xlsx");           % Current censity/potential data from Damjanovic

% Extracted data from the Excel files 

% Scohy
Scohy_potential = Scohy_Ir_data(1:end,1);                                   % [V vs RHE] - Potential
Scohy_current_density = Scohy_Ir_data(1:end,2)*10^(-3+4);                   % [A/m^2] - Current density, originally in mA/cm^2
Scohy_T = 25 + 273;                                                         % [K] - Temperature
Scohy_a_H_plus = 0.5*2;                                                     % [-] - Activity of H+


%Damjanovic
Damjanovic_potential = Damjanovic_Ir_data(1:end,2);                         % [V vs RHE]          
Damjanovic_current_density = Damjanovic_Ir_data(1:end,1)*10^4;              % [A/m^2] - Originally A/cm^2 (The article contains log(i), but I converted thm
Damjanovic_T = 25 + 273;                                                    % [K] - Temperature
Damjanovic_a_H_plus = 1;                                                    % [-] - Activity of H+


%% %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting the expression of the current based on r_2 to the data
% The r_2_fit returns the curve (fitting results) and the gof.
% The coefficients are contained in the curve
% The r_fit is based on the approximation of the large equation from
% Reksten

% Scohy
[Scohy_curve, Scohy_gof] = ...                                              % This is the expression with rds
    r_2_fit(Scohy_potential, Scohy_current_density, Scohy_a_H_plus, Scohy_T, "Linear");

[Scohy_curve_fit, Scohy_gof_fit] = ...                                      % This is the expression with no rds, but k_4 << 1
    r_fit(Scohy_potential, Scohy_current_density, Scohy_a_H_plus, Scohy_T, "Linear");

% Damjanovic
[Damjanovic_curve, Damjanovic_gof] = ...                                    % This is the expression with rds
    r_2_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Linear");


[Damjanovic_curve_fit, Damjanovic_gof_fit] = ...                            % This is the expression with no rds, but k_4 << 1
    r_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Linear");

% Damjanovic (log)
[Damjanovic_log_curve, Damjanovic_log_gof] = ...                            % This is the expression with rds
    r_2_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Logarithmic");

[Damjanovic_log_curve_fit, Damjanovic_log_gof_fit] = ...                    % This is the expression with no rds, but k_4 << 1
    r_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Logarithmic");

%% %%%%%%%%%%%%%%%%%%%%%%%% REKSTEN %%%%%%%%%%%%%%%%%%%%% %%

% Introducing the fitting from the large equation by using the separate
% fits from earlier

[Scohy_reksten_curve, Scohy_reksten_gof] = r_reksten_fit(Scohy_curve, Scohy_potential, Scohy_current_density, Scohy_a_H_plus, Scohy_T, "Linear");
[Damjanovic_reksten_curve, Damjanovic_reksten_gof] = r_reksten_fit(Damjanovic_curve, Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Linear");
[Damjanovic_log_reksten_curve, Damjanovic_log_reksten_gof] = r_reksten_fit(Damjanovic_curve, Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Logarithmic");



%current_density_scohy = current_density_reksten(Scohy_curve, Scohy_reksten_curve.k_1_0_plus, Scohy_reksten_curve.k_2_0, Scohy_reksten_curve.k_4_0_plus, Scohy_a_H_plus, Scohy_T, Scohy_potential);
%current = current_density_reksten(curve, k_1_0_plus, k_2_0, k_4_0_plus, a_H_plus, T, E)


% Comments while looking at the rate constants
%   k_1_0_plus: The chemical rate constant for the first step
%               This one has the highest uncertainty. Does not affect curve
%               at all with increasing values (10^1 --> 10^10), while the others are based
%               on Scohy or 10^1.
%               When it goes lower than 10 (10^1 --> 10^-3) it goes below
%               the Scohy curve while the other values are 10^1 or based on
%               Scohy
%
%   k_2_0: The ration between the chemical rate constants for the second
%          step - k_2_0_minus/k_2_0_plus
%          Very vulnerable to an increase in the value above 10^1. Reasonable since the Scohy predicts a k_2_0_plus = 0.116
%          and it is reasonable to assume that the forward rate is much
%          faster than the backward rate --> k_2_0 should be smaller
%          Values below 10^1 (10^-1 and 10^0.5 gives very low curve...)?
%           This rate constant also decides how k_4_0 changes
%
%   k_4_0_plus: The chemical rate constant for the forward rate of the last
%               step.
%               values above 10^1 gives the same graph if the other
%               constants are based on Scohy or 10^1.
%               Lower than 10^1 while the others are Scohy or 10^1 gives
%               absurd curve
%
%
% ---- Å endre alle de tre konstantene til 10^2 ga en overraskende bra
% kurve... jeg tipper det må være muligå kjøre en fitting på denne kurven
% her. Kanskje det gir meg mye usikkerhet, men det må være smartere enn å
% sjekke manuelt 10^0, 10^2, 10^3 -- Best



%% Scohy for reference
k_2_0_scohy = Scohy_reksten_curve.k_2_0;
k_1_0_plus_scohy = Scohy_reksten_curve.k_1_0_plus;
k_4_0_plus_scohy = Scohy_reksten_curve.k_4_0_plus;

k_4_0_plus_array = [10^(-8), 10^(-6), 10^(-4), 10^(-2), 10^0]*k_4_0_plus_scohy;
k_1_0_plus_array = [10^-5, 10^-4, 10^-3, 10^-2, 10^0]*k_1_0_plus_scohy;
k_2_0_array = [10^(-8), 10^(-5), 10^0, 10^5, 10^8]*k_2_0_scohy;

figure()
plot(Scohy_potential, Scohy_current_density)
hold on
for i = 1:length(k_2_0_array)
   current_density = current_density_reksten(Scohy_curve, k_1_0_plus_scohy, k_2_0_array(i), k_4_0_plus_scohy, Scohy_a_H_plus, Scohy_T, Scohy_potential);
   plot(Scohy_potential, current_density)
end
hold off
legend("Data", "small", "next small", "middle", "next big", "big")
xlabel("potential")
ylabel("current density")
title("Scohy - k_2_0")

figure()
plot(Scohy_potential, Scohy_current_density)
hold on
for i = 1:length(k_4_0_plus_array)
   current_density = current_density_reksten(Scohy_curve, k_1_0_plus_scohy, k_2_0_scohy, k_4_0_plus_array(i), Scohy_a_H_plus, Scohy_T, Scohy_potential);
   plot(Scohy_potential, current_density)
end
%plot(Scohy_reksten_curve)
hold off
legend("Data", "small", "next small", "middle", "next big", "big")
xlabel("potential")
ylabel("current density")
title("Scohy - k_4_0_plus")

figure()
plot(Scohy_potential, Scohy_current_density)
hold on
for i = 1:length(k_1_0_plus_array)
   current_density = current_density_reksten(Scohy_curve, k_1_0_plus_array(i), k_2_0_scohy, k_4_0_plus_scohy, Scohy_a_H_plus, Scohy_T, Scohy_potential);
   plot(Scohy_potential, current_density)
end
hold off
legend("Data", "small", "next small", "middle", "next big", "big")
xlabel("potential")
ylabel("current density")
title("Scohy - k_1_0_plus")

figure()
plot(Damjanovic_potential, Damjanovic_current_density)
hold on
plot(Damjanovic_reksten_curve)
hold off
xlabel("potential")
ylabel("current density")

figure()
plot(Damjanovic_potential, log10(Damjanovic_current_density))
hold on
plot(Damjanovic_log_reksten_curve)
hold off
xlabel("potential")
ylabel("current density")

%% %%%%%%%%%% Marshall %%%%%%%%%%%% %%