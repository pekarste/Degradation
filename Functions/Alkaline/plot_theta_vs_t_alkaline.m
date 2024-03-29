% This will be a plotting script where I will plot theta_2 as a function of
% time

%% %%%%%%%%%%%%%%%%%%%%%% ALKALINE model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will couple all the different functions together and work
% like a masterscript. 

%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_OER_SHE = 0.40;                                                           % Standard reduction potential for OER vs SHE - alkaline
E_REF_SHE = -0.829;                                                         % Standard redcution potential for HER vs SHE - alkaline
E_n = E_OER_SHE - E_REF_SHE;                                                % Standard reduction potential for OER vs RHE - alkaline
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
theta_2_0 = eps;

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Guess for k_4_0_plus
k_4_0_plus = [10^(-1) 10^(-2) 10^(-3)];

% Data used for fitting r2_alkaline
Cherevko_alkaline = readmatrix("Data\Alkaline\Cherevko_alkaline_polarisation_data.xlsx");% Potential/current density data from Cherevko
Damjanovic_alkaline = readmatrix("Data\Alkaline\Damjanovic_alkaline_polarisation_3.xlsx"); % Current density/potential data from Damjanovic
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

Schalenbach_dissolution_CV_linear_data = readmatrix("Data\Alkaline\Schalenbach\Schalenbach_dissolution_alkaline_peak_1.xlsx");
% Schalenbach dissolution vs time - [ng/cm^2*s]

Schalenbach_dissolution_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,2);    % Schalenbach dissolution data - [ng/cm^2*s]
Schalenbach_time_CV_linear = Schalenbach_dissolution_CV_linear_data(1:end,1);           % Schalenbach time data - [s]
Schalenbach_dissolution_mole = Schalenbach_dissolution_CV_linear*10^(-9)*10^(4)/Mm_Ir;  % Changes the units from ng/cm^2*s --> mole/m^2*s

% Creating string elements for x and y labels
x_label_string = '$t$ - [s]';
y_label_string = '$\Gamma\theta_{2}$(t) - [mol m$^{-2}$]';
y_label_string_2 = '$E$ - [V vs RHE]';

% Colour blind pallette
Orange          = [.90 .60 .0];                                        % Orange                                        
Reddish_purple  = [.80 .60 .70];                                       % Reddish purple
Sky_blue        = [.35 .70 .90];                                       % Sky blue
%% ################## Theta vs Time #############################

% Transforming time to potential fot the interpolation
potential_interpol_alkaline = CV_potential_alkaline(Schalenbach_time_CV_linear, "array");

% Creating a string element for the legends
string_array_1 = sprintf('$k^{0}_{4+}$ = %.1f $s^{-1}$', round(k_4_0_plus(1), 5));
string_array_2 = sprintf('$k^{0}_{4+}$ = %.2f $s^{-1}$', round(k_4_0_plus(2), 5));
string_array_3 = sprintf('$k^{0}_{4+}$ = %.3f $s^{-1}$', round(k_4_0_plus(3), 5));
string_array_4 = 'E(t)';

%% Cherevko
[t_cherevko_1, gamma_theta_cherevko_1, potential_cherevko_1, theta_cherevko_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(1));
[t_cherevko_2, gamma_theta_cherevko_2, potential_cherevko_2, theta_cherevko_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(2));
[t_cherevko_3, gamma_theta_cherevko_3, potential_cherevko_3, theta_cherevko_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(3));


figure('Name', 'Cherevko: theta_2 vs time')                                     % Creating figure
%yyaxis left
plot(t_cherevko_1, gamma_theta_cherevko_1, "Color", Orange)                            % Plots the line for 1
hold on
scatter(Schalenbach_time_CV_linear, theta_cherevko_interpol_1,...               % Scatter interpolated values for 1
    45,Orange, 'o','filled')                                                                      
plot(t_cherevko_2, gamma_theta_cherevko_2, "Color", Reddish_purple)                           % Plots the line for 2
scatter(Schalenbach_time_CV_linear, theta_cherevko_interpol_2,...               % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_cherevko_3, gamma_theta_cherevko_3, "Color", Sky_blue)                          % Plots the line for 3
scatter(Schalenbach_time_CV_linear, theta_cherevko_interpol_3,...               % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  

ax_cherevko_alkaline = gca; % current axes                                      % Creating an ax with gca such that the fontsize can be changed
ax_cherevko_alkaline.TickDir = "out";
box off
%hold off

ax_cherevko_alkaline.XAxis.FontSize = 15;                                       % Changing the tick size on the x-axis
ax_cherevko_alkaline(1).YAxis.FontSize = 15;                                       % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([min(gamma_theta_cherevko_3)*0 max(gamma_theta_cherevko_3)])
yline(max(gamma_theta_cherevko_3))

yyaxis right
ax_cherevko_alkaline.YAxis(2).FontSize = 15;
ax_cherevko_alkaline.YAxis(2).Color = 'black';
plot(Schalenbach_time_CV_linear, potential_interpol_alkaline,...                         % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')

ylabel(y_label_string_2,'Interpreter','latex')                                % Label for second y_axis
%ylim([min(potential_interpol_alkaline)*0 max(potential_interpol_alkaline)])

% annotation('textbox', [.15 .80 .1 .1], 'String',["Cherevko -", "Alkaline"],...  % Creating an annotation, textbox, with the rsquare value from the cfit
%     'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
% legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...  % Creating a legend for the graphs
%     'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.15 .75 .1 .1], 'String',["Cherevko -", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.60 .14 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Orange);
annotation('textbox', [.55 .50 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Reddish_purple, 'Rotation',-45);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Sky_blue, 'Rotation', -20);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% Damjanovic
[t_damj_1, gamma_theta_damj_1, potential_damj_1, theta_damj_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(1));
[t_damj_2, gamma_theta_damj_2, potential_damj_2, theta_damj_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(2));
[t_damj_3, gamma_theta_damj_3, potential_damj_3, theta_damj_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(3));


figure('Name', 'Damjanovic: theta_2 vs time')                                   % Creating figure
%yyaxis left
plot(t_damj_1, gamma_theta_damj_1, "Color", Orange)                                    % Plots the line for 1
hold on
scatter(Schalenbach_time_CV_linear, theta_damj_interpol_1,...                   % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_damj_2, gamma_theta_damj_2, "Color", Reddish_purple)                             % Plots the line for 2
scatter(Schalenbach_time_CV_linear, theta_damj_interpol_2,...                   % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_damj_3, gamma_theta_damj_3, "Color", Sky_blue)                            % Plots the line for 3
scatter(Schalenbach_time_CV_linear, theta_damj_interpol_3,...                   % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_damj_alkaline = gca; % current axes                                          % Creating an ax with gca such that the fontsize can be changed
ax_damj_alkaline.TickDir = 'out';
box off
ax_damj_alkaline.XAxis.FontSize = 15;                                           % Changing the tick size on the x-axis
ax_damj_alkaline.YAxis.FontSize = 15;                                           % Changing the tick size on the y-axis

xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([0, max(gamma_theta_damj_3)])
yline(max(gamma_theta_damj_3))

yyaxis right
ax_damj_alkaline.YAxis(2).FontSize = 15;
ax_damj_alkaline.YAxis(2).Color = 'black';
plot(Schalenbach_time_CV_linear, potential_interpol_alkaline,...                         % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
%ylim([min(gamma_theta_damj_3)*0 max(gamma_theta_damj_3)])
ylabel(y_label_string_2,'Interpreter','latex')                                % Label for second y_axis
% annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic -", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
%     'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
% legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...  % Creating a legend for the graphs
%     'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic -", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.60 .14 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Orange);
annotation('textbox', [.55 .50 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Reddish_purple, 'Rotation',-45);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Sky_blue, 'Rotation', -20);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% Damjanovic log
[t_damj_log_1, gamma_theta_damj_log_1, potential_damj_log_1, theta_damj_log_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(1));
[t_damj_log_2, gamma_theta_damj_log_2, potential_damj_log_2, theta_damj_log_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(2));
[t_damj_log_3, gamma_theta_damj_log_3, potential_damj_log_3, theta_damj_log_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(3));


figure('Name', 'Damjanovic log: theta_2 vs time')                                   % Creating figure
%yyaxis left
plot(t_damj_log_1, gamma_theta_damj_log_1, "Color", Orange)                                    % Plots the line for 1
hold on
scatter(Schalenbach_time_CV_linear, theta_damj_log_interpol_1,...                   % Scatter interpolated values for 1
    45, Orange, 'o', 'filled')                                                                      
plot(t_damj_log_2, gamma_theta_damj_log_2, "Color", Reddish_purple)                                   % Plots the line for 2
scatter(Schalenbach_time_CV_linear, theta_damj_log_interpol_2,...                   % Scatter interpolated values for 2
    45, Reddish_purple, 'square', 'filled')     
plot(t_damj_log_3, gamma_theta_damj_log_3, "Color", Sky_blue)                                  % Plots the line for 3
scatter(Schalenbach_time_CV_linear, theta_damj_log_interpol_3,...                   % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_damj_log_alkaline = gca; % current axes                                          % Creating an ax with gca such that the fontsize can be changed
ax_damj_log_alkaline.TickDir = 'out';
box off
ax_damj_log_alkaline.XAxis.FontSize = 15;                                           % Changing the tick size on the x-axis
ax_damj_log_alkaline.YAxis(1).FontSize = 15;                                           % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([min(gamma_theta_damj_log_3)*0 max(gamma_theta_damj_log_3)])
yline(max(gamma_theta_damj_log_3))

yyaxis right
ax_damj_log_alkaline.YAxis(2).FontSize = 15;
ax_damj_log_alkaline.YAxis(2).Color = 'black';
plot(Schalenbach_time_CV_linear, potential_interpol_alkaline,...                         % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel(y_label_string_2,'Interpreter','latex')                                % Label for second y_axis
% annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic log-", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
%     'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
% legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...  % Creating a legend for the graphs
%     'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic log-", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.60 .17 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Orange);
annotation('textbox', [.57 .52 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Reddish_purple, 'Rotation',-45);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Sky_blue, 'Rotation', -20);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------
%
%% Schalenbach
[t_schalenbach_1, gamma_theta_schalenbach_1, potential_schalenbach_1, theta_schalenbach_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(1));
[t_schalenbach_2, gamma_theta_schalenbach_2, potential_schalenbach_2, theta_schalenbach_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(2));
[t_schalenbach_3, gamma_theta_schalenbach_3, potential_schalenbach_3, theta_schalenbach_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(3));


figure('Name', 'Schalenbach: theta_2 vs time')                                   % Creating figure
%yyaxis left
plot(t_schalenbach_1, gamma_theta_schalenbach_1, "Color", Orange)                       % Plots the line for 1
hold on
scatter(Schalenbach_time_CV_linear, theta_schalenbach_interpol_1,...             % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_schalenbach_2, gamma_theta_schalenbach_2, "Color", Reddish_purple)                      % Plots the line for 2
scatter(Schalenbach_time_CV_linear, theta_schalenbach_interpol_2,...             % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_schalenbach_3, gamma_theta_schalenbach_3, "Color", Sky_blue)                     % Plots the line for 3
scatter(Schalenbach_time_CV_linear, theta_schalenbach_interpol_3,...             % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_schalenbach_alkaline = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_schalenbach_alkaline.TickDir = 'out';
box off
ax_schalenbach_alkaline.XAxis.FontSize = 15;                                     % Changing the tick size on the x-axis
ax_schalenbach_alkaline.YAxis(1).FontSize = 15;                                     % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([min(gamma_theta_schalenbach_3)*0 max(gamma_theta_schalenbach_3)])
yline(max(gamma_theta_schalenbach_3))

yyaxis right
ax_schalenbach_alkaline.YAxis(2).FontSize = 15;                                     % Changing the tick size on the y-axis
ax_schalenbach_alkaline.YAxis(2).Color = 'black';                                     % Changing the tick size on the y-axis
plot(Schalenbach_time_CV_linear, potential_interpol_alkaline,...                         % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel(y_label_string_2,'Interpreter','latex')                                % Label for second y_axis
% annotation('textbox', [.15 .80 .1 .1], 'String',["Schalenbach-", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
%     'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
% legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...  % Creating a legend for the graphs
%     'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.15 .80 .1 .1], 'String',["Schalenbach-", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.60 .15 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Orange);
annotation('textbox', [.55 .5 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Reddish_purple, 'Rotation',-50);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',Sky_blue, 'Rotation', -20);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',16, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% Schalenbach

figure('Name', 'Schalenbach: theta_2 vs time')                                   % Creating figure
%yyaxis left
plot(t_schalenbach_1, gamma_theta_schalenbach_1, "Color", "red")                       % Plots the line for 1
hold on
scatter(Schalenbach_time_CV_linear, theta_schalenbach_interpol_1,...             % Scatter interpolated values for 1
    45,"red", 'o', 'filled')                                                                      
plot(t_schalenbach_2, gamma_theta_schalenbach_2, "Color", "blue")                      % Plots the line for 2
scatter(Schalenbach_time_CV_linear, theta_schalenbach_interpol_2,...             % Scatter interpolated values for 2
    45,"blue", 'square', 'filled')     
plot(t_schalenbach_3, gamma_theta_schalenbach_3, "Color", "green")                     % Plots the line for 3
scatter(Schalenbach_time_CV_linear, theta_schalenbach_interpol_3,...             % Scatter the interpolated values for 3
    45, "green", 'diamond', 'filled')  
hold off
ax_schalenbach_alkaline = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_schalenbach_alkaline.XAxis.FontSize = 12;                                     % Changing the tick size on the x-axis
ax_schalenbach_alkaline.YAxis.FontSize = 12;                                     % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Schalenbach_time_CV_linear(1) Schalenbach_time_CV_linear(end)])
ylim([min(gamma_theta_schalenbach_3)*0 max(gamma_theta_schalenbach_3)])

yyaxis right
%ax_2_damj_alkaline = gca;
%ax_2_damj_alkaline.YAxis.FontSize = 12;
plot(Schalenbach_time_CV_linear, Schalenbach_dissolution_mole,...                         % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex')                                % Label for second y_axis
annotation('textbox', [.15 .80 .1 .1], 'String',["Schalenbach-", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
legend({'',string_array_1, '', string_array_2, '', string_array_3, "$r_{diss}$"},...  % Creating a legend for the graphs
    'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------