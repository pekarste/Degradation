% This will be a plotting script where I will plot r_3 as a function of
% time with dissolution rate

%% %%%%%%%%%%%%%%%%%%%%%% ALKALINE model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will couple all the different functions together and work
% like a masterscript. 


% Will try to use acidic equations but alkaline environemnt
%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_OER_SHE = 0.40;                                                           % Standard reduction potential for OER vs SHE - alkaline
E_REF_SHE = -0.829;                                                         % Standard redcution potential for HER vs SHE - alkaline
E_n = E_OER_SHE - E_REF_SHE;                                                % Standard reduction potential for OER vs RHE - alkaline
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
%gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
theta_2_0 = eps;

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Guess for k_4_0_plus
k_4_0_plus = [10^(-1) 10^(-2) 10^(-3)];

% Data used for fitting r2_alkaline
Cherevko_alkaline = readmatrix("Data\Alkaline\Cherevko_alkaline_polarisation_data.xlsx");% Potential/current density data from Cherevko
Damjanovic_alkaline = readmatrix("Data\Alkaline\Damjanovic_alkaline_polarisation.xlsx"); % Current density/potential data from Damjanovic
%--------------------------------------------------------------------------
Schalenbach_polarisation_alkaline = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_polarisation_curve_alkaline.xlsx");% Polarisation curve from alkaline dissolution
%--------------------------------------------------------------------------
%% Extracted data from the Excel files 

% Alkaline - Cherevko
Cherevko_E_alkaline = Cherevko_alkaline(1:end,1);                           % [V vs RHE] - Potential
Cherevko_i_alkaline = Cherevko_alkaline(1:end,2)*10^(-3+4);                 % [A/m^2] - Current density, originally in mA/cm^2
Cherevko_T_alkaline = 25 + 273;                                             % [K] - Temperature
Cherevko_OH_alkaline = 0.05*2;                                              % [-] - Activity of OH-


% Alkaline - Damjanovic
Damjanovic_E_alkaline = Damjanovic_alkaline(1:end,2);                       % [V vs RHE]          
Damjanovic_i_alkaline = Damjanovic_alkaline(1:end,1)*10^4;                  % [A/m^2] - Originally A/cm^2 (The article contains log(i), but I converted thm
Damjanovic_T_alkaline = 25 + 273;                                           % [K] - Temperature
Damjanovic_OH_alkaline = 1;                                                 % [-] - Activity of OH-
%_-------------------------------------------------------------------------
% Alkaline - Schalenbach
Schalenbach_E_alkaline = Schalenbach_polarisation_alkaline(1:end,1);           % [V vs RHE] - Potential
Schalenbach_i_alkaline = Schalenbach_polarisation_alkaline(1:end,2)*10^(-3+4); % [A/m^2] - Current density, originally in mA/cm^2
Schalenbach_T_alkaline = 25 + 273;                                             % [K] - Temperature
Schalenbach_OH_alkaline = 0.05*1;                                              % [-] - Activity of OH-
Schalenbach_sweep_rate = 2*10^(-3);                                            % [V/s] -Schalenbach Sweep rate 
%--------------------------------------------------------------------------

Schalenbach_dissolution_CV_linear_data = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_dissolution_linear_alkaline.xlsx");
% Schalenbach dissolution vs time - [ng/cm^2*s]

Schalenbach_dissolution_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,2);    % Schalenbach dissolution data - [ng/cm^2*s]
Schalenbach_time_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,1);           % Schalenbach time data - [s]
Schalenbach_dissolution_mole = Schalenbach_dissolution_CV_linear*10^(-9)*10^(4)/Mm_Ir;  % Changes the units from ng/cm^2*s --> mole/m^2*s

%% ################## Theta vs Time #############################

% Transforming time to potential fot the interpolation
potential_interpol = CV_potential_alkaline(Schalenbach_time_CV_linear, "array");

% Creating a string element for the legends
string_array_1 = sprintf('$r_{3}$($k^{0}_{4+}$ = %.2f)', round(k_4_0_plus(1), 5));
string_array_2 = sprintf('$r_{3}$($k^{0}_{4+}$ = %.3f)', round(k_4_0_plus(2), 5));
string_array_3 = sprintf('$r_{3}$($k^{0}_{4+}$ = %.4f)', round(k_4_0_plus(3), 5));
string_array_4 = '$\frac{d Ir}{d t}$';

%% ####################     Cherevko     ################################

% Finding necessary values
[t_cherevko_1, gamma_theta_cherevko_1, potential_cherevko_1, theta_cherevko_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(1));
[t_cherevko_2, gamma_theta_cherevko_2, potential_cherevko_2, theta_cherevko_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(2));
[t_cherevko_3, gamma_theta_cherevko_3, potential_cherevko_3, theta_cherevko_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(3));

% Fitting to the degradation data
cherevko_1_curve = chi_square_alkaline(t_cherevko_1, gamma_theta_cherevko_1, Schalenbach_OH_alkaline, k_4_0_plus(1), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
cherevko_2_curve = chi_square_alkaline(t_cherevko_2, gamma_theta_cherevko_2, Schalenbach_OH_alkaline, k_4_0_plus(2), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
cherevko_3_curve = chi_square_alkaline(t_cherevko_3, gamma_theta_cherevko_3, Schalenbach_OH_alkaline, k_4_0_plus(3), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);

cherevko_1_string = sprintf('$k^{0}_{3+}$ = %.3g', cherevko_1_curve.k_3_0_plus);
cherevko_2_string = sprintf('$k^{0}_{3+}$ = %.3g', cherevko_2_curve.k_3_0_plus);
cherevko_3_string = sprintf('$k^{0}_{3+}$ = %.3g', cherevko_3_curve.k_3_0_plus);

% Plotting
figure('Name','Cherevko: r_3')
fig_cherevko_1 = plot(cherevko_1_curve, 'red');
hold on
fig_cherevko_2 = plot(cherevko_2_curve, "green");                           % Creating a fig to stor the plot of the curve fit (cfit element)
fig_cherevko_3 = plot(cherevko_3_curve, "blue");                            % Creating a fig to stor the plot of the curve fit (cfit element)
plot(Schalenbach_time_CV_linear, Schalenbach_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')

set(fig_cherevko_1,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_cherevko_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_cherevko_3,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfithold off

ax_cherevko_alkaline = gca; % current axes                                  % Creating an ax with gca such that the fontsize can be changed
ax_cherevko_alkaline.XAxis.FontSize = 15;                                   % Changing the tick size on the x-axis
ax_cherevko_alkaline.YAxis.FontSize = 15;                                   % Changing the tick size on the y-axis

annotation('textbox', [.15 .80 .1 .1], 'String',["Cherevko -", "Alkaline"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)


xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('$r_{3}$ - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([0 17*10^(-9)])

annotation('textbox', [.46 .52 .1 .1], 'String',cherevko_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', 0);
annotation('textbox', [.62 .54 .1 .1], 'String',cherevko_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -3);
annotation('textbox', [.65 .38 .1 .1], 'String',cherevko_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -25);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% #####################     Damjanovic       #############################

% Finding necessary values
[t_damj_1, gamma_theta_damj_1, potential_damj_1, theta_damj_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(1));
[t_damj_2, gamma_theta_damj_2, potential_damj_2, theta_damj_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(2));
[t_damj_3, gamma_theta_damj_3, potential_damj_3, theta_damj_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(3));

% Fitting to the degradation data
damj_1_curve = chi_square_alkaline(t_damj_1, gamma_theta_damj_1, Schalenbach_OH_alkaline, k_4_0_plus(1), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
damj_2_curve = chi_square_alkaline(t_damj_2, gamma_theta_damj_2, Schalenbach_OH_alkaline, k_4_0_plus(2), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
damj_3_curve = chi_square_alkaline(t_damj_3, gamma_theta_damj_3, Schalenbach_OH_alkaline, k_4_0_plus(3), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);

damj_1_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_1_curve.k_3_0_plus);
damj_2_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_2_curve.k_3_0_plus);
damj_3_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_3_curve.k_3_0_plus);

% Plotting
figure('Name', 'Damjanovic: r_3')                                           % Creating figure
fig_damj_1 = plot(damj_1_curve, 'red');
hold on
fig_damj_2 = plot(damj_2_curve, "green");                                   % Creating a fig to stor the plot of the curve fit (cfit element)
fig_damj_3 = plot(damj_3_curve, "blue");                                    % Creating a fig to stor the plot of the curve fit (cfit element)
plot(Schalenbach_time_CV_linear, Schalenbach_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')

set(fig_damj_1,'lineWidth',1);                                              % Changing the linewidth of the curve of the cfit
set(fig_damj_2,'lineWidth',1);                                              % Changing the linewidth of the curve of the cfit
set(fig_damj_3,'lineWidth',1);                                              % Changing the linewidth of the curve of the cfithold off

ax_damj_alkaline = gca; % current axes                                      % Creating an ax with gca such that the fontsize can be changed
ax_damj_alkaline.XAxis.FontSize = 15;                                       % Changing the tick size on the x-axis
ax_damj_alkaline.YAxis.FontSize = 15;                                       % Changing the tick size on the y-axis

annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic -", "Alkaline"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)

xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('$r_{3}$ - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([0 17*10^(-9)])

annotation('textbox', [.46 .15 .1 .1], 'String',damj_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', 0);
annotation('textbox', [.62 .46 .1 .1], 'String',damj_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -5);
annotation('textbox', [.62 .30 .1 .1], 'String',damj_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -25);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% ###################    Damjanovic log     ##############################

% Finding necessary values
[t_damj_log_1, gamma_theta_damj_log_1, potential_damj_log_1, theta_damj_log_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(1));
[t_damj_log_2, gamma_theta_damj_log_2, potential_damj_log_2, theta_damj_log_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(2));
[t_damj_log_3, gamma_theta_damj_log_3, potential_damj_log_3, theta_damj_log_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(3));

% Fitting to the degradation data
damj_log_1_curve = chi_square_alkaline(t_damj_log_1, gamma_theta_damj_log_1, Schalenbach_OH_alkaline, k_4_0_plus(1), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
damj_log_2_curve = chi_square_alkaline(t_damj_log_2, gamma_theta_damj_log_2, Schalenbach_OH_alkaline, k_4_0_plus(2), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
damj_log_3_curve = chi_square_alkaline(t_damj_log_3, gamma_theta_damj_log_3, Schalenbach_OH_alkaline, k_4_0_plus(3), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);

damj_log_1_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_log_1_curve.k_3_0_plus);
damj_log_2_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_log_2_curve.k_3_0_plus);
damj_log_3_string = sprintf('$k^{0}_{3+}$ = %.3g', damj_log_3_curve.k_3_0_plus);

% Plotting
figure('Name', 'Damjanovic log: r_3')                                       % Creating figure
fig_damj_log_1 = plot(damj_log_1_curve, 'red');
hold on
fig_damj_log_2 = plot(damj_log_2_curve, "green");                           % Creating a fig to stor the plot of the curve fit (cfit element)
fig_damj_log_3 = plot(damj_log_3_curve, "blue");                            % Creating a fig to stor the plot of the curve fit (cfit element)
plot(Schalenbach_time_CV_linear, Schalenbach_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')

set(fig_damj_log_1,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_damj_log_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_damj_log_3,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfithold off

ax_damj_log_alkaline = gca; % current axes                                  % Creating an ax with gca such that the fontsize can be changed
ax_damj_log_alkaline.XAxis.FontSize = 15;                                   % Changing the tick size on the x-axis
ax_damj_log_alkaline.YAxis.FontSize = 15;                                   % Changing the tick size on the y-axis

annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic log -", "Alkaline"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)

xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('$r_{3}$ - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([0 17*10^(-9)])

annotation('textbox', [.46 .15 .1 .1], 'String',damj_log_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', 0);
annotation('textbox', [.62 .46 .1 .1], 'String',damj_log_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -20);
annotation('textbox', [.62 .30 .1 .1], 'String',damj_log_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -30);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% #########################     Schalenbach     ##########################

% Finding necessary values
[t_schalenbach_1, gamma_theta_schalenbach_1, potential_schalenbach_1, theta_schalenbach_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(1));
[t_schalenbach_2, gamma_theta_schalenbach_2, potential_schalenbach_2, theta_schalenbach_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(2));
[t_schalenbach_3, gamma_theta_schalenbach_3, potential_schalenbach_3, theta_schalenbach_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(3));

% Fitting to the degradation data
schalenbach_1_curve = chi_square_alkaline(t_schalenbach_1, gamma_theta_schalenbach_1, Schalenbach_OH_alkaline, k_4_0_plus(1), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
schalenbach_2_curve = chi_square_alkaline(t_schalenbach_2, gamma_theta_schalenbach_2, Schalenbach_OH_alkaline, k_4_0_plus(2), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);
schalenbach_3_curve = chi_square_alkaline(t_schalenbach_3, gamma_theta_schalenbach_3, Schalenbach_OH_alkaline, k_4_0_plus(3), Schalenbach_time_CV_linear, Schalenbach_dissolution_mole);

schalenbach_1_string = sprintf('$k^{0}_{3+}$ = %.3g', schalenbach_1_curve.k_3_0_plus);
schalenbach_2_string = sprintf('$k^{0}_{3+}$ = %.3g', schalenbach_2_curve.k_3_0_plus);
schalenbach_3_string = sprintf('$k^{0}_{3+}$ = %.3g', schalenbach_3_curve.k_3_0_plus);

% Plotting
figure('Name', 'Schalenbach: r_3')                                       % Creating figure
fig_schalenbach_1_2 = plot(schalenbach_1_curve, 'red');
hold on
fig_schalenbach_2_2 = plot(schalenbach_2_curve, "green");                        % Creating a fig to stor the plot of the curve fit (cfit element)
fig_schalenbach_3_2 = plot(schalenbach_3_curve, "blue");                         % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_schalenbach_1_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_schalenbach_2_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
set(fig_schalenbach_3_2,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfithold off
plot(Schalenbach_time_CV_linear, Schalenbach_dissolution_mole,...           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')
ax_schalenbach_alkaline = gca; % current axes                                  % Creating an ax with gca such that the fontsize can be changed
ax_schalenbach_alkaline.XAxis.FontSize = 15;                                   % Changing the tick size on the x-axis
ax_schalenbach_alkaline.YAxis.FontSize = 15;                                   % Changing the tick size on the y-axis
annotation('textbox', [.15 .80 .1 .1], 'String',["Schalenbach -", "Alkaline"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({string_array_1, string_array_2, string_array_3, ''},...                                                             % Creating a legend for the graphs
    'Position', [.67 .77 .1 .1], 'Interpreter','latex', 'FontSize',15)
xlabel('Time t - [s]','Interpreter','latex', 'FontSize', 15)                % Creating x-label
ylabel('Dissolution - [$\frac{mol}{m^{2}s}$]',...                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([0 17*10^(-9)])

annotation('textbox', [.46 .15 .1 .1], 'String',schalenbach_1_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','red', 'Rotation', 0);
annotation('textbox', [.62 .46 .1 .1], 'String',schalenbach_2_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','blue', 'Rotation', -5);
annotation('textbox', [.62 .30 .1 .1], 'String',schalenbach_3_string,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color','green', 'Rotation', -30);
annotation('textbox', [.37 .65 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------