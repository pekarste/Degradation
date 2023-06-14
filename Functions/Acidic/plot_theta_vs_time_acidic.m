% This will be a plotting script where I will plot theta_2 as a function of
% time

%% %%%%%%%%%%%%%%%%%%%%%% ACIDIC model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
theta_2_0 = eps;

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Guess for k_4_0_plus
k_4_0_plus = [10^(-1) 10^(-2) 10^(-3)];

% Data used for fitting r2_alkaline
Scohy_acidic = readmatrix("Data\Acidic\Scohy_activated_Ir_LSV.xlsx");% Potential/current density data from Scohy
Damjanovic_acidic = readmatrix("Data\Acidic\Damjanovic_Ir_E_vs_log_i_acidic.xlsx"); % Current density/potential data from Damjanovic
Cherevko_acidic = readmatrix("Data\Acidic\Cherevko_polarisation.xlsx");% Polarisation curve from Cherevko acidic
%--------------------------------------------------------------------------
Mayrhofer_acidic = readmatrix("Data\Acidic\Mayrhofer_polarisation_cycle_20.xlsx");% Polarisation curve for the iridium cycled 20 times and the same used in the degradation with slow sweep rate
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

% Acidc - Mayrhofer
Mayrhofer_E_acidic = Mayrhofer_acidic(1:end,1);                             % [V vs RHE] - Potential
Mayrhofer_i_acidic = Mayrhofer_acidic(1:end,2)*10^(-3+4);                   % [A/m^2] - Current density, originally in mA/cm^2
Mayrhofer_T_acidic = 25+273.13;                                             % [K] - Temperature
Mayrhofer_a_H_plus = 0.1*2;                                                 % Concentration of H+ (0.1 M H2SO4)
%--------------------------------------------------------------------------

%Mayrhofer_potential_data = readmatrix("Mayrhofer_potential.xlsx");          % Mayrhofer potential vs time data - Don't really use this... computes it based on scan rate and initial potential
Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_slow_second.xlsx");    % Mayrhofer dissolution vs time data - [ng/cm^2s]

Mayrhofer_dissolution = Mayrhofer_dissolution_data(4:end-4,2);                % Mayrhofer dissolution data - [ng/cm^2*s] -- Starting from 5 to remove the tail
Mayrhofer_time = Mayrhofer_dissolution_data(4:end-4,1);                       % Mayrhofer time data [s] -- Starting from 5 to remove the tail to be consistent

Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s

sweep_rate = 2*10^(-3);                                                    % Mayrhofer Sweep rate [V/s]
%Mayrhofer_a_H_plus = 0.1*2;                                                 % Concentration of H+ (0.1 M H2SO4)
Mayrhofer_T = 25 + 273.13;                                                  % mayrhofer states room temperature%% ################## Theta vs Time #############################

%-------------------------------------------------------------------------
% Transforming time to potential fot the interpolation
potential_interpol_acidic = CV_potential_acidic(Mayrhofer_time, "array");

% Creating a string element for the legends
string_array_1 = sprintf('$k^{0}_{4+}$ = %.1f', round(k_4_0_plus(1), 5));
string_array_2 = sprintf('$k^{0}_{4+}$ = %.2f', round(k_4_0_plus(2), 5));
string_array_3 = sprintf('$k^{0}_{4+}$ = %.3f', round(k_4_0_plus(3), 5));
string_array_4 = 'E(t)';

% Creating string elements for x and y labels
x_label_string = '$t$ - [s]';
y_label_string = '$\Gamma\theta_{2}$(t) - [mol m$^{-2}$]';
y_label_string_2 = '$E$ - [V vs RHE]';

% Colour blind pallette
Orange          = [.90 .60 .0];                                        % Orange                                        
Reddish_purple  = [.80 .60 .70];                                       % Reddish purple
Sky_blue        = [.35 .70 .90];                                       % Sky blue
%% Scohy
[t_scohy_1, gamma_theta_scohy_1, potential_scohy_1, gamma_theta_scohy_interpol_1] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(1));
[t_scohy_2, gamma_theta_scohy_2, potential_scohy_2, gamma_theta_scohy_interpol_2] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(2));
[t_scohy_3, gamma_theta_scohy_3, potential_scohy_3, gamma_theta_scohy_interpol_3] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(3));


figure('Name', 'Scohy: theta_2 vs time')                                    % Creating figure
%yyaxis left
plot(t_scohy_1, gamma_theta_scohy_1, "Color", Orange)                        % Plots the line for 1
hold on
scatter(Mayrhofer_time, gamma_theta_scohy_interpol_1,...                    % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_scohy_2, gamma_theta_scohy_2, "Color", Reddish_purple)                       % Plots the line for 2
scatter(Mayrhofer_time, gamma_theta_scohy_interpol_2,...                    % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_scohy_3, gamma_theta_scohy_3, "Color", Sky_blue)                      % Plots the line for 3
scatter(Mayrhofer_time, gamma_theta_scohy_interpol_3,...                    % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_scohy_acidic = gca; % current axes                                       % Creating an ax with gca such that the fontsize can be changed
ax_scohy_acidic.XAxis.FontSize = 15;                                        % Changing the tick size on the x-axis
ax_scohy_acidic.YAxis(1).FontSize = 15;                                        % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([min(gamma_theta_scohy_3)*0 max(gamma_theta_scohy_3)])

yyaxis right
plot(Mayrhofer_time, potential_interpol_acidic,...                          % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ax_scohy_acidic.YAxis(2).FontSize = 15;                                        % Changing the tick size on the y-axis
ax_scohy_acidic.YAxis(2).Color = 'black';                                        % Changing the tick size on the y-axis
ylabel(y_label_string_2,'Interpreter','latex')                              % Label for second y_axis
%legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...  % Creating a legend for the graphs
%    'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)
ylim([min(potential_interpol_acidic) max(potential_interpol_acidic)])

annotation('textbox', [.15 .55 .1 .1], 'String',["Scohy -", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.60 .14 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange);
annotation('textbox', [.35 .40 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% Damjanovic
[t_damj_1, gamma_theta_damj_1, potential_damj_1, gamma_theta_damj_interpol_1] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(1));
[t_damj_2, gamma_theta_damj_2, potential_damj_2, gamma_theta_damj_interpol_2] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(2));
[t_damj_3, gamma_theta_damj_3, potential_damj_3, gamma_theta_damj_interpol_3] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(3));


figure('Name', 'Damjanovic: theta_2 vs time')                               % Creating figure
%yyaxis left
plot(t_damj_1, gamma_theta_damj_1, "Color", Orange)                          % Plots the line for 1
hold on
scatter(Mayrhofer_time, gamma_theta_damj_interpol_1,...                     % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_damj_2, gamma_theta_damj_2, "Color", Reddish_purple)                         % Plots the line for 2
scatter(Mayrhofer_time, gamma_theta_damj_interpol_2,...                     % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_damj_3, gamma_theta_damj_3, "Color", Sky_blue)                        % Plots the line for 3
scatter(Mayrhofer_time, gamma_theta_damj_interpol_3,...                     % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_damj_acidic = gca; % current axes                                        % Creating an ax with gca such that the fontsize can be changed
ax_damj_acidic.XAxis.FontSize = 15;                                         % Changing the tick size on the x-axis
ax_damj_acidic.YAxis(1).FontSize = 15;                                         % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([min(gamma_theta_damj_3)*0 max(gamma_theta_damj_3)])

yyaxis right
plot(Mayrhofer_time, potential_interpol_acidic,...                          % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ax_damj_acidic.YAxis(2).FontSize = 15;                                         % Changing the tick size on the y-axis
ax_damj_acidic.YAxis(2).Color = 'black';                                         % Changing the tick size on the y-axis
ylabel(y_label_string_2,'Interpreter','latex')            % Label for second y_axis
ylim([min(potential_interpol_acidic) max(potential_interpol_acidic)])
%annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic -", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
%    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
%legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...% Creating a legend for the graphs
%    'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.52 .32 .1 .1], 'String',["Damjanovic -", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.55 .17 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange);
annotation('textbox', [.40 .50 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% Damjanovic log
[t_damj_log_1, gamma_theta_damj_log_1, potential_damj_log_1, gamma_theta_damj_log_interpol_1] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(1));
[t_damj_log_2, gamma_theta_damj_log_2, potential_damj_log_2, gamma_theta_damj_log_interpol_2] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(2));
[t_damj_log_3, gamma_theta_damj_log_3, potential_damj_log_3, gamma_theta_damj_log_interpol_3] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(3));


figure('Name', 'Damjanovic log: theta_2 vs time')                           % Creating figure
%yyaxis left
plot(t_damj_log_1, gamma_theta_damj_log_1, "Color", Orange)                  % Plots the line for 1
hold on
scatter(Mayrhofer_time, gamma_theta_damj_log_interpol_1,...                 % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_damj_log_2, gamma_theta_damj_log_2, "Color", Reddish_purple)                 % Plots the line for 2
scatter(Mayrhofer_time, gamma_theta_damj_log_interpol_2,...                 % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_damj_log_3, gamma_theta_damj_log_3, "Color", Sky_blue)                % Plots the line for 3
scatter(Mayrhofer_time, gamma_theta_damj_log_interpol_3,...                 % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_damj_log_acidic = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_damj_log_acidic.XAxis.FontSize = 15;                                     % Changing the tick size on the x-axis
ax_damj_log_acidic.YAxis(1).FontSize = 15;                                     % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([min(gamma_theta_damj_log_3)*0 max(gamma_theta_damj_log_3)])

yyaxis right
plot(Mayrhofer_time, potential_interpol_acidic,...                          % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ax_damj_log_acidic.YAxis(2).FontSize = 15;                                     % Changing the tick size on the y-axis
ax_damj_log_acidic.YAxis(2).Color = 'black';                                     % Changing the tick size on the y-axis
ylabel(y_label_string_2,'Interpreter','latex')                              % Label for second y_axis
ylim([min(potential_interpol_acidic) max(potential_interpol_acidic)])
% annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic log-", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
%     'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
% legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...  % Creating a legend for the graphs
%     'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.52 .32 .1 .1], 'String',["Damjanovic log-", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.55 .17 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange);
annotation('textbox', [.40 .50 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------


%% Cherevko
[t_cherevko_1, gamma_theta_cherevko_1, potential_cherevko_1, gamma_theta_cherevko_interpol_1] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(1));
[t_cherevko_2, gamma_theta_cherevko_2, potential_cherevko_2, gamma_theta_cherevko_interpol_2] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(2));
[t_cherevko_3, gamma_theta_cherevko_3, potential_cherevko_3, gamma_theta_cherevko_interpol_3] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(3));


figure('Name', 'Cherevko: theta_2 vs time')                                 % Creating figure
%yyaxis left
plot(t_cherevko_1, gamma_theta_cherevko_1, "Color", Orange)                  % Plots the line for 1
hold on
scatter(Mayrhofer_time, gamma_theta_cherevko_interpol_1,...                 % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_cherevko_2, gamma_theta_cherevko_2, "Color", Reddish_purple)                 % Plots the line for 2
scatter(Mayrhofer_time, gamma_theta_cherevko_interpol_2,...                 % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_cherevko_3, gamma_theta_cherevko_3, "Color", Sky_blue)                % Plots the line for 3
scatter(Mayrhofer_time, gamma_theta_cherevko_interpol_3,...                 % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_cherevko_acidic = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_cherevko_acidic.XAxis.FontSize = 15;                                     % Changing the tick size on the x-axis
ax_cherevko_acidic.YAxis(1).FontSize = 15;                                     % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([min(gamma_theta_cherevko_3)*0 max(gamma_theta_cherevko_3)])

yyaxis right
plot(Mayrhofer_time, potential_interpol_acidic,...                          % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ax_cherevko_acidic.YAxis(2).FontSize = 15;                                     % Changing the tick size on the y-axis
ax_cherevko_acidic.YAxis(2).Color = 'black';                                     % Changing the tick size on the y-axis
ylabel(y_label_string_2,'Interpreter','latex')                              % Label for second y_axis
ylim([min(potential_interpol_acidic) max(potential_interpol_acidic)])
% annotation('textbox', [.15 .80 .1 .1], 'String',["Cherevko-", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
%     'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
% legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},... % Creating a legend for the graphs
%     'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.52 .35 .1 .1], 'String',["Cherevko -", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.50 .17 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange);
annotation('textbox', [.37 .55 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% Cherevko

figure('Name', 'Cherevko: theta_2 vs time')                                 % Creating figure
%yyaxis left
plot(t_cherevko_1, gamma_theta_cherevko_1, "Color", Orange)                  % Plots the line for 1
hold on
scatter(Mayrhofer_time, gamma_theta_cherevko_interpol_1,...                 % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_cherevko_2, gamma_theta_cherevko_2, "Color", Reddish_purple)                 % Plots the line for 2
scatter(Mayrhofer_time, gamma_theta_cherevko_interpol_2,...                 % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_cherevko_3, gamma_theta_cherevko_3, "Color", Sky_blue)                % Plots the line for 3
scatter(Mayrhofer_time, gamma_theta_cherevko_interpol_3,...                 % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
hold off
ax_cherevko_acidic = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_cherevko_acidic.XAxis.FontSize = 12;                                     % Changing"green" the tick size on the x-axis
ax_cherevko_acidic.YAxis.FontSize = 12;                                     % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
%ylim([min(gamma_theta_cherevko_3)*0 max(gamma_theta_cherevko_3)])

yyaxis right
%ax_2_damj_alkaline = gca;
%ax_2_damj_alkaline.YAxis.FontSize = 12;
plot(Mayrhofer_time, Mayrhofer_dissolution_mole,...                         % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel(y_label_string_2,'Interpreter','latex')                              % Label for second y_axis
% annotation('textbox', [.15 .80 .1 .1], 'String',["Cherevko-", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
%     'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
% legend({'',string_array_1, '', string_array_2, '', string_array_3, "$r_{diss}$"},...  % Creating a legend for the graphs
%     'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)


%--------------------------------------------------------------------------

%% Mayrhofer
[t_mayrhofer_1, gamma_theta_mayrhofer_1, potential_mayrhofer_1, gamma_theta_mayrhofer_interpol_1] =...
    time_theta_potential_ode15s_acidic(Mayrhofer_E_acidic, Mayrhofer_i_acidic, Mayrhofer_a_H_plus, Mayrhofer_T_acidic, "Linear", k_4_0_plus(1));
[t_mayrhofer_2, gamma_theta_mayrhofer_2, potential_mayrhofer_2, gamma_theta_mayrhofer_interpol_2] =...
    time_theta_potential_ode15s_acidic(Mayrhofer_E_acidic, Mayrhofer_i_acidic, Mayrhofer_a_H_plus, Mayrhofer_T_acidic, "Linear", k_4_0_plus(2));
[t_mayrhofer_3, gamma_theta_mayrhofer_3, potential_mayrhofer_3, gamma_theta_mayrhofer_interpol_3] =...
    time_theta_potential_ode15s_acidic(Mayrhofer_E_acidic, Mayrhofer_i_acidic, Mayrhofer_a_H_plus, Mayrhofer_T_acidic, "Linear", k_4_0_plus(3));

figure('Name', 'Mayrhofer: theta_2 vs time')                                    % Creating figure
%yyaxis left
plot(t_mayrhofer_1, gamma_theta_mayrhofer_1, "Color", Orange)                        % Plots the line for 1
hold on
scatter(Mayrhofer_time, gamma_theta_mayrhofer_interpol_1,...                    % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(t_mayrhofer_2, gamma_theta_mayrhofer_2, "Color", Reddish_purple)                       % Plots the line for 2
scatter(Mayrhofer_time, gamma_theta_mayrhofer_interpol_2,...                    % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(t_mayrhofer_3, gamma_theta_mayrhofer_3, "Color", Sky_blue)                      % Plots the line for 3
scatter(Mayrhofer_time, gamma_theta_mayrhofer_interpol_3,...                    % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_mayrhofer_acidic = gca; % current axes                                       % Creating an ax with gca such that the fontsize can be changed
ax_mayrhofer_acidic.XAxis.FontSize = 15;                                        % Changing the tick size on the x-axis
ax_mayrhofer_acidic.YAxis(1).FontSize = 15;                                        % Changing the tick size on the y-axis
xlabel(x_label_string,'Interpreter','latex')
ylabel(y_label_string,'Interpreter','latex')
xlim([Mayrhofer_time(1) Mayrhofer_time(end)])
ylim([min(gamma_theta_mayrhofer_3)*0 max(gamma_theta_mayrhofer_3)])

yyaxis right
plot(Mayrhofer_time, potential_interpol_acidic,...                          % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ax_mayrhofer_acidic.YAxis(2).FontSize = 15;                                        % Changing the tick size on the y-axis
ax_mayrhofer_acidic.YAxis(2).Color = 'black';                                        % Changing the tick size on the y-axis
ylabel(y_label_string_2,'Interpreter','latex')                              % Label for second y_axis
%legend({'',string_array_1, '', string_array_2, '', string_array_3, "E(t)"},...  % Creating a legend for the graphs
%    'Position', [.2375 .55 .1 .1],'Interpreter','latex', 'FontSize',15)
ylim([min(potential_interpol_acidic) max(potential_interpol_acidic)])

annotation('textbox', [.15 .55 .1 .1], 'String',["Mayrhofer -", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.60 .14 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange);
annotation('textbox', [.35 .40 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple);
annotation('textbox', [.60 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue);
annotation('textbox', [.70 .55 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',[.5 .5 .5]);
