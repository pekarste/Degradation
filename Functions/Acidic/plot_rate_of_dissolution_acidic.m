% This will be a plotting script where I will plot r_3 as a function of
% time with dissolution rate

%% %%%%%%%%%%%%%%%%%%%%%% Acidic model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will couple all the different functions together and work
% like a masterscript. 

%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_OER_SHE = 1.229;                                                          % Standard reduction potential for OER vs SHE - acidic
E_REF_RHE = 0.0;                                                            % Standard redcution potential for HER vs SHE - acidic
E_n = E_OER_SHE - E_REF_RHE;                                                % Standard reduction potential for OER vs RHE
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
theta_2_0 = eps;

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Guess for k_4_0_plus
k_4_0_plus = [10^(-1) 10^(-2) 10^(-3)];

% Data used for fitting r2_acidic
Scohy_acidic = readmatrix("Data\Acidic\Scohy_activated_Ir_LSV.xlsx");% Potential/current density data from Scohy
Damjanovic_acidic = readmatrix("Data\Acidic\Damjanovic_Ir_E_vs_log_i.xlsx"); % Current density/potential data from Damjanovic
Cherevko_acidic = readmatrix("Data\Acidic\Cherevko_polarisation.xlsx");% Polarisation curve from Cherevko acidic
%--------------------------------------------------------------------------
%% Extracted data from the Excel files 

% Acidic - Scohy
Scohy_E_acidic = Scohy_acidic(1:end,1);                                     % [V vs RHE] - Potential
Scohy_i_acidic = Scohy_acidic(1:end,2)*10^(-3+4);                           % [A/m^2] - Current density, originally in mA/cm^2
Scohy_T_acidic = 25 + 273;                                                  % [K] - Temperature
Scohy_a_H_plus = 0.5*2;                                                     % [-] - Activity of H+

% Acidic - Damjanovic
Damjanovic_E_acidic = Damjanovic_acidic(1:end,2);                           % [V vs RHE]          
Damjanovic_i_acidic = Damjanovic_acidic(1:end,1)*10^4;                      % [A/m^2] - Originally A/cm^2 (The article contains log(i), but I converted thm
Damjanovic_T_acidic = 25 + 273;                                             % [K] - Temperature
Damjanovic_a_H_plus = 1;                                                    % [-] - Activity of H+

% Acidic - Cherevko
Cherevko_E_acidic = Cherevko_acidic(1:end,1);                               % [V vs RHE] - Potential
Cherevko_i_acidic = Cherevko_acidic(1:end,2)*10^(-3+4);                     % [A/m^2] - Current density, originally in mA/cm^2
Cherevko_T_acidic = 25 + 273;                                               % [K] - Temperature
Cherevko_a_H_plus = 0.1*2;                                                  % [-] - Activity of H+

%--------------------------------------------------------------------------

Mayrhofer_potential_data = readmatrix("Mayrhofer_potential.xlsx");          % Mayrhofer potential vs time data - Don't really use this... computes it based on scan rate and initial potential
Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_2.xlsx");    % Mayrhofer dissolution vs time data - [ng/cm^2s]

Mayrhofer_dissolution = Mayrhofer_dissolution_data(5:end,2);                % Mayrhofer dissolution data - [ng/cm^2*s] -- Starting from 5 to remove the tail
Mayrhofer_time = Mayrhofer_dissolution_data(5:end,1);                       % Mayrhofer time data [s] -- Starting from 5 to remove the tail to be consistent

Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s

sweep_rate = 10*10^(-3);                                                    % Mayrhofer Sweep rate [V/s]
Mayrhofer_a_H_plus = 0.1*2;                                                 % Concentration of H+ (0.1 M H2SO4)
Mayrhofer_T = 25 + 273.13;                                                  % mayrhofer states room temperature%% ################## Theta vs Time #############################

%-------------------------------------------------------------------------
% Transforming time to potential fot the interpolation
potential_interpol_acidic = CV_potential_acidic(Mayrhofer_time, "array");

% Creating a string element for the legends
string_array_1 = sprintf('$r_{3}$($k^{0}_{4+}$ = %.2f)', round(k_4_0_plus(1), 5));
string_array_2 = sprintf('$r_{3}$($k^{0}_{4+}$ = %.3f)', round(k_4_0_plus(2), 5));
string_array_3 = sprintf('$r_{3}$($k^{0}_{4+}$ = %.4f)', round(k_4_0_plus(3), 5));
string_array_4 = '$\frac{d Ir}{d t}$';

%% ####################     Scohy     ################################

% Finding necessary values
[t_scohy_1, gamma_theta_scohy_1, potential_scohy_1, theta_scohy_interpol_1] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(1));
[t_scohy_2, gamma_theta_scohy_2, potential_scohy_2, theta_scohy_interpol_2] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(2));
[t_scohy_3, gamma_theta_scohy_3, potential_scohy_3, theta_scohy_interpol_3] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(3));

% Fitting to the degradation data
scohy_1_curve = chi_square_acidic(t_scohy_1, gamma_theta_scohy_1, k_4_0_plus(1), Mayrhofer_time, Mayrhofer_dissolution_mole);
scohy_2_curve = chi_square_acidic(t_scohy_2, gamma_theta_scohy_2, k_4_0_plus(2), Mayrhofer_time, Mayrhofer_dissolution_mole);
scohy_3_curve = chi_square_acidic(t_scohy_3, gamma_theta_scohy_3, k_4_0_plus(3), Mayrhofer_time, Mayrhofer_dissolution_mole);

scohy_1_string = sprintf('$k^{0}_{3+}$ = %.3g', scohy_1_curve.k_3_0_plus);
scohy_2_string = sprintf('$k^{0}_{3+}$ = %.3g', scohy_2_curve.k_3_0_plus);
scohy_3_string = sprintf('$k^{0}_{3+}$ = %.3g', scohy_3_curve.k_3_0_plus);

% Plotting
figure('Name','Scohy: r_3')
fig_scohy_1 = plot(scohy_1_curve, 'red');
hold on
fig_scohy_2 = plot(scohy_2_curve, "green");                           % Creating a fig to stor the plot of the curve fit (cfit element)
fig_scohy_3 = plot(scohy_3_curve, "blue");                            % Creating a fig to stor the plot of the curve fit (cfit element)
plot(Mayrhofer_time, Mayrhofer_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')

set(fig_scohy_1,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_scohy_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_scohy_3,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfithold off

ax_scohy_acidic = gca; % current axes                                  % Creating an ax with gca such that the fontsize can be changed
ax_scohy_acidic.XAxis.FontSize = 15;                                   % Changing the tick size on the x-axis
ax_scohy_acidic.YAxis.FontSize = 15;                                   % Changing the tick size on the y-axis

annotation('textbox', [.15 .80 .1 .1], 'String',["Scohy -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)

xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('$r_{3}$ - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([0 3.5*10^(-8)])

annotation('textbox', [.46 .52 .1 .1], 'String',scohy_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', -10);
annotation('textbox', [.62 .54 .1 .1], 'String',scohy_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -3);
annotation('textbox', [.65 .38 .1 .1], 'String',scohy_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -20);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% #####################     Damjanovic       #############################

% Finding necessary values
[t_damj_1, gamma_theta_damj_1, potential_damj_1, theta_damj_interpol_1] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(1));
[t_damj_2, gamma_theta_damj_2, potential_damj_2, theta_damj_interpol_2] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(2));
[t_damj_3, gamma_theta_damj_3, potential_damj_3, theta_damj_interpol_3] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(3));

% Fitting to the degradation data
damj_1_curve = chi_square_acidic(t_damj_1, gamma_theta_damj_1, k_4_0_plus(1), Mayrhofer_time, Mayrhofer_dissolution_mole);
damj_2_curve = chi_square_acidic(t_damj_2, gamma_theta_damj_2, k_4_0_plus(2), Mayrhofer_time, Mayrhofer_dissolution_mole);
damj_3_curve = chi_square_acidic(t_damj_3, gamma_theta_damj_3, k_4_0_plus(3), Mayrhofer_time, Mayrhofer_dissolution_mole);

damj_1_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_1_curve.k_3_0_plus);
damj_2_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_2_curve.k_3_0_plus);
damj_3_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_3_curve.k_3_0_plus);

% Plotting
figure('Name', 'Damjanovic: r_3')                                           % Creating figure
fig_damj_1 = plot(damj_1_curve, 'red');
hold on
fig_damj_2 = plot(damj_2_curve, "green");                                   % Creating a fig to stor the plot of the curve fit (cfit element)
fig_damj_3 = plot(damj_3_curve, "blue");                                    % Creating a fig to stor the plot of the curve fit (cfit element)
plot(Mayrhofer_time, Mayrhofer_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')

set(fig_damj_1,'lineWidth',1);                                              % Changing the linewidth of the curve of the cfit
set(fig_damj_2,'lineWidth',1);                                              % Changing the linewidth of the curve of the cfit
set(fig_damj_3,'lineWidth',1);                                              % Changing the linewidth of the curve of the cfithold off

ax_damj_acidic = gca; % current axes                                      % Creating an ax with gca such that the fontsize can be changed
ax_damj_acidic.XAxis.FontSize = 15;                                       % Changing the tick size on the x-axis
ax_damj_acidic.YAxis.FontSize = 15;                                       % Changing the tick size on the y-axis

annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)

xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('$r_{3}$ - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([0 3.5*10^(-8)])
annotation('textbox', [.46 .15 .1 .1], 'String',damj_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', -10);
annotation('textbox', [.62 .46 .1 .1], 'String',damj_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -3);
annotation('textbox', [.62 .30 .1 .1], 'String',damj_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -20);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% ###################    Damjanovic log     ##############################

% Finding necessary values
[t_damj_log_1, gamma_theta_damj_log_1, potential_damj_log_1, theta_damj_log_interpol_1] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(1));
[t_damj_log_2, gamma_theta_damj_log_2, potential_damj_log_2, theta_damj_log_interpol_2] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(2));
[t_damj_log_3, gamma_theta_damj_log_3, potential_damj_log_3, theta_damj_log_interpol_3] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(3));

% Fitting to the degradation data
damj_log_1_curve = chi_square_acidic(t_damj_log_1, gamma_theta_damj_log_1, k_4_0_plus(1), Mayrhofer_time, Mayrhofer_dissolution_mole);
damj_log_2_curve = chi_square_acidic(t_damj_log_2, gamma_theta_damj_log_2, k_4_0_plus(2), Mayrhofer_time, Mayrhofer_dissolution_mole);
damj_log_3_curve = chi_square_acidic(t_damj_log_3, gamma_theta_damj_log_3, k_4_0_plus(3), Mayrhofer_time, Mayrhofer_dissolution_mole);

damj_log_1_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_log_1_curve.k_3_0_plus);
damj_log_2_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_log_2_curve.k_3_0_plus);
damj_log_3_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_log_3_curve.k_3_0_plus);

% Plotting
figure('Name', 'Damjanovic log: r_3')                                       % Creating figure
fig_damj_log_1 = plot(damj_log_1_curve, 'red');
hold on
fig_damj_log_2 = plot(damj_log_2_curve, "green");                           % Creating a fig to stor the plot of the curve fit (cfit element)
fig_damj_log_3 = plot(damj_log_3_curve, "blue");                            % Creating a fig to stor the plot of the curve fit (cfit element)
plot(Mayrhofer_time, Mayrhofer_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')

set(fig_damj_log_1,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_damj_log_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_damj_log_3,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfithold off

ax_damj_log_acidic = gca; % current axes                                  % Creating an ax with gca such that the fontsize can be changed
ax_damj_log_acidic.XAxis.FontSize = 15;                                   % Changing the tick size on the x-axis
ax_damj_log_acidic.YAxis.FontSize = 15;                                   % Changing the tick size on the y-axis

annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic log -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)

xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('$r_{3}$ - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([0 3.5*10^(-8)])

annotation('textbox', [.46 .15 .1 .1], 'String',damj_log_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', -10);
annotation('textbox', [.62 .46 .1 .1], 'String',damj_log_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -3);
annotation('textbox', [.62 .30 .1 .1], 'String',damj_log_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -20);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% #########################     Cherevko     ##########################

% Finding necessary values
[t_cherevko_1, gamma_theta_cherevko_1, potential_cherevko_1, theta_cherevko_interpol_1] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(1));
[t_cherevko_2, gamma_theta_cherevko_2, potential_cherevko_2, theta_cherevko_interpol_2] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(2));
[t_cherevko_3, gamma_theta_cherevko_3, potential_cherevko_3, theta_cherevko_interpol_3] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(3));

% Fitting to the degradation data
cherevko_1_curve = chi_square_acidic(t_cherevko_1, gamma_theta_cherevko_1, k_4_0_plus(1), Mayrhofer_time, Mayrhofer_dissolution_mole);
cherevko_2_curve = chi_square_acidic(t_cherevko_2, gamma_theta_cherevko_2, k_4_0_plus(2), Mayrhofer_time, Mayrhofer_dissolution_mole);
cherevko_3_curve = chi_square_acidic(t_cherevko_3, gamma_theta_cherevko_3, k_4_0_plus(3), Mayrhofer_time, Mayrhofer_dissolution_mole);

cherevko_1_string = sprintf('$k^{0}_{3+}$ = %.3g', cherevko_1_curve.k_3_0_plus);
cherevko_2_string = sprintf('$k^{0}_{3+}$ = %.3g', cherevko_2_curve.k_3_0_plus);
cherevko_3_string = sprintf('$k^{0}_{3+}$ = %.3g', cherevko_3_curve.k_3_0_plus);

% Plotting
figure('Name', 'Cherevko: r_3')                                       % Creating figure
fig_cherevko_1_2 = plot(cherevko_1_curve, 'red');
hold on
fig_cherevko_2_2 = plot(cherevko_2_curve, "green");                        % Creating a fig to stor the plot of the curve fit (cfit element)
fig_cherevko_3_2 = plot(cherevko_3_curve, "blue");                         % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_cherevko_1_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_cherevko_2_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_cherevko_3_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfithold off
plot(Mayrhofer_time, Mayrhofer_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')
ax_cherevko_acidic = gca; % current axes                                  % Creating an ax with gca such that the fontsize can be changed
ax_cherevko_acidic.XAxis.FontSize = 12;                                   % Changing the tick size on the x-axis
ax_cherevko_acidic.YAxis.FontSize = 12;                                   % Changing the tick size on the y-axis
annotation('textbox', [.15 .80 .1 .1], 'String',["Cherevko -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)
xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('Dissolution - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([0 3.5*10^(-8)])

annotation('textbox', [.46 .15 .1 .1], 'String',cherevko_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', -5);
annotation('textbox', [.62 .46 .1 .1], 'String',cherevko_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -3);
annotation('textbox', [.62 .30 .1 .1], 'String',cherevko_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -20);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------