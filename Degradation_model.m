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
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Guess for k_4_0_plus
k_4_0_plus = 10^-2;

% Data used for fitting of r2_acidic
Scohy_acidic = readmatrix("Data\Acidic\Scohy_activated_Ir_LSV.xlsx");       % Potential/current density data from Scohy
Damjanovic_acidic = readmatrix("Data\Acidic\Damjanovic_Ir_E_vs_log_i_acidic.xlsx");% Current density/potential data from Damjanovic
Cherevko_acidic = readmatrix("Data\Acidic\Cherevko_polarisation.xlsx");     % Potential/current density data from Cherevko
%--------------------------------------------------------------------------
Mayrhofer_acidic = readmatrix("Data\Acidic\Mayrhofer_polarisation_cycle_20.xlsx");% Polarisation curve for the iridium cycled 20 times and the same used in the degradation with slow sweep rate
%% Extracted data from the Excel files 

% Acidic - Scohy
Scohy_potential = Scohy_acidic(1:end,1);                                    % [V vs RHE] - Potential
Scohy_current_density = Scohy_acidic(1:end,2)*10^(-3+4);                    % [A/m^2] - Current density, originally in mA/cm^2
Scohy_T = 25 + 273;                                                         % [K] - Temperature
Scohy_a_H_plus = 0.5*2;                                                     % [-] - Activity of H+

% Acidic - Damjanovic
Damjanovic_potential = Damjanovic_acidic(1:end,2);                          % [V vs RHE]          
Damjanovic_current_density = Damjanovic_acidic(1:end,1)*10^4;               % [A/m^2] - Originally A/cm^2 (The article contains log(i), but I converted thm
Damjanovic_T = 25 + 273;                                                    % [K] - Temperature
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
%% %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting the expression of the current based on r_2 to the data
% The r_2_fit returns the curve (fitting results) and the gof.
% The coefficients are contained in the curve

% Scohy
[Scohy_curve, Scohy_gof] = ...                                              % This is the expression with rds
    r_2_fit_acidic(Scohy_potential, Scohy_current_density, Scohy_a_H_plus, Scohy_T, "Linear");

[Scohy_curve_fit, Scohy_gof_fit] = ...                                      % This is the expression with no rds, but k_4 << 1
    r_fit(Scohy_potential, Scohy_current_density, Scohy_a_H_plus, Scohy_T, "Linear");

% Damjanovic
[Damjanovic_curve, Damjanovic_gof] = ...                                    % This is the expression with rds
    r_2_fit_acidic(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Linear");


[Damjanovic_curve_fit, Damjanovic_gof_fit] = ...                            % This is the expression with no rds, but k_4 << 1
    r_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Linear");

% Damjanovic (log)
[Damjanovic_log_curve, Damjanovic_log_gof] = ...                            % This is the expression with rds
    r_2_fit_acidic(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Logarithmic");

[Damjanovic_log_curve_fit, Damjanovic_log_gof_fit] = ...                    % This is the expression with no rds, but k_4 << 1
    r_fit(Damjanovic_potential, Damjanovic_current_density, Damjanovic_a_H_plus, Damjanovic_T, "Logarithmic");

% Cherevko - Acidic
[Cherevko_curve_acidic, Cherevko_gof_acidic] = ...                                              % This is the expression with rds
    r_2_fit_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear");

% Cherevko - Acidic
[Mayrhofer_curve_acidic, Mayrhofer_gof_acidic] = ...                                              % This is the expression with rds
    r_2_fit_acidic(Mayrhofer_E_acidic, Mayrhofer_i_acidic, Mayrhofer_a_H_plus, Mayrhofer_T_acidic, "Linear");
% Plotting the fitting with the data

%-------------------------------------------------------------------------
% Scohy - Acidic
figure("Name", "Scohy Fitting Acidic")                                                                    % Creates figure

scatter(Scohy_potential, Scohy_current_density, 45, "filled", "blue", 'o')  % Scatter plot of the sampled values from Scohy
hold on
fig_scohy_fit = plot(Scohy_curve, "black");                                 % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_scohy_fit,'lineWidth',1);                                           % Changing the linewidth of the curve of the cfit
yline(Scohy_current_density(end))
xline(Scohy_potential(end))
%hold off

legend({'Data', 'Fitting'},...                                              % Creating a legend for the graphs
    'Position', [.2 .65 .1 .1], 'Interpreter','latex', 'FontSize',15)

str_scohy = sprintf("$R^{2}$ = %.5f", round(Scohy_gof.rsquare, 5));         % Creating a string element for the annotation
annotation('textbox', [.15 .8 .1 .1], 'String',str_scohy,...                % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);


ax_scohy = gca; % current axes                                              % Creating an ax with gca such that the fontsize can be changed
ax_scohy.TickDir = 'out';
box off
ax_scohy.XAxis.FontSize = 15;                                               % Changing the tick size on the x-axis
ax_scohy.YAxis.FontSize = 15;                                               % Changing the tick size on the y-axis

xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize', 15)  % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
annotation('textbox', [.15 .80 .1 .1], 'String',["Scohy -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
xlim([Scohy_potential(1) Scohy_potential(end)])
ylim([Scohy_current_density(1)*0 Scohy_current_density(end)])


% Damjanovic - Acidic
figure("Name", "Damjanovic Fitting Acidic")                                                                    % Creates figure

scatter(Damjanovic_potential, Damjanovic_current_density, 45,...            % Scatter plot of the sampled values from Damjanovic
    [0.4940 0.1840 0.5560], "^", "filled")
hold on
fig_damjanovic = plot(Damjanovic_curve, "black");                           % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_damjanovic,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
yline(Damjanovic_current_density(end))
xline(Damjanovic_potential(end))
hold off

legend({'Data', 'Fitting'},...                                              % Creating a legend for the graphs
    'Position', [.2 .65 .1 .1], 'Interpreter','latex', 'FontSize',15)
str_damjanovic = ...                                                        % Creating a string element for the annotation
    sprintf("$R^{2}$ = %.5f", round(Damjanovic_gof.rsquare, 5));
annotation('textbox', [.15 .8 .1 .1], 'String',str_damjanovic,...           % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex','FitBoxToText','on', 'FontSize',15);

ax_damjanovic = gca;
ax_damjanovic.TickDir = 'out';
box off
ax_damjanovic.XAxis.FontSize = 15;                                          % Changing the tick size on the x-axis
ax_damjanovic.YAxis.FontSize = 15;                                          % Changing the tick size on the y-axis

xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize',15)   % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                % Creating y-label
    'Interpreter','latex', 'FontSize',15)
annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
xlim([Damjanovic_potential(1) Damjanovic_potential(end)])
ylim([Damjanovic_current_density(1)*0 Damjanovic_current_density(end)])



% Damjanovic log - Acidic
figure("Name", "Damjanovic log Fitting Acidic")                                                                    % Creating figure

scatter(Damjanovic_potential, log10(Damjanovic_current_density), 45,...     % Scatter plot of the sampled values from Damjanovic
    [0.9290 0.6940 0.1250],"v", "filled")
hold on
fig_damjanovic_log = plot(Damjanovic_log_curve, "black");                   % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_damjanovic_log,'lineWidth',1);                                      % Changing the linewidth of the curve of the cfit
yline(log10(Damjanovic_current_density(end)))
xline(Damjanovic_potential(end))
hold off

legend({'Data', 'Fitting'}, 'Position', [.2 .65 .1 .1],...                  % Creating a legend for the graphs
    'Interpreter','latex', 'FontSize',15)
str_damjanovic_log = ...                                                    % Creating a string element for the annotation
    sprintf("$R^{2}$ = %.5f", round(Damjanovic_log_gof.rsquare, 5));
annotation('textbox', [.15 .8 .1 .1], 'String',str_damjanovic_log,...       % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);

ax_damjanovic_log = gca; % current axes                                         % Creating an ax with gca such that the fontsize can be changed
ax_damjanovic_log.TickDir = 'out';
box off
ax_damjanovic_log.XAxis.FontSize = 15;                                      % Changing the tick size on the x-axis
ax_damjanovic_log.YAxis.FontSize = 15;                                      % Changing the tick size on the y-axis

xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize',15)   % Creating x-label
ylabel('$\log_{10}$ of current density - $\log{i}$/[$Am^{-2}$]',...         % Creating y-label
    'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.15 .80 .1 .1], 'String',["Damjanovic log -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
xlim([Damjanovic_potential(1) Damjanovic_potential(end)])
ylim([-5 log10(Damjanovic_current_density(end))])


% Cherevko - Acidic
figure("Name", "Cherevko Fitting Acidic")                                                                    % Creating figure
scatter(Cherevko_E_acidic, Cherevko_i_acidic, 45,...     % Scatter plot of the sampled values from Damjanovic
    "filled", "red", "square")
hold on
fig_cherevko_acidic = plot(Cherevko_curve_acidic, "black");                   % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_cherevko_acidic,'lineWidth',1);                                      % Changing the linewidth of the curve of the cfit
yline(Cherevko_i_acidic(end))
xline(Cherevko_E_acidic(end))
hold off

legend({'Data', 'Fitting'}, 'Position', [.2 .65 .1 .1],...                  % Creating a legend for the graphs
    'Interpreter','latex', 'FontSize',15)
str_cherevko_acidic = ...                                                    % Creating a string element for the annotation
    sprintf("$R^{2}$ = %.5f", round(Cherevko_gof_acidic.rsquare, 5));
annotation('textbox', [.15 .8 .1 .1], 'String',str_cherevko_acidic,...       % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);

ax_cherevko_acidic = gca;
ax_cherevko_acidic.TickDir = 'out';
box off
ax_cherevko_acidic.XAxis.FontSize = 15;                                      % Changing the tick size on the x-axis
ax_cherevko_acidic.YAxis.FontSize = 15;                                      % Changing the tick size on the y-axis

xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize',15)   % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                % Creating y-label
    'Interpreter','latex', 'FontSize',15)
annotation('textbox', [.15 .80 .1 .1], 'String',["Cherevko -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);

xlim([Cherevko_E_acidic(1) Cherevko_E_acidic(end)])
ylim([Cherevko_i_acidic(1)*0 Cherevko_i_acidic(end)])

% Mayrhofer - Acidic
figure("Name", "Mayrhofer Fitting Acidic")                                                                    % Creating figure
scatter(Mayrhofer_E_acidic, Mayrhofer_i_acidic, 45,...     % Scatter plot of the sampled values from Damjanovic
    "filled", "red", "square")
hold on

fig_mayrhofer_acidic = plot(Mayrhofer_curve_acidic, "black");                   % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_mayrhofer_acidic,'lineWidth',1);                                      % Changing the linewidth of the curve of the cfit
yline(Mayrhofer_i_acidic(end))
xline(Mayrhofer_E_acidic(end))
hold off

legend({'Data', 'Fitting'}, 'Position', [.2 .65 .1 .1],...                  % Creating a legend for the graphs
    'Interpreter','latex', 'FontSize',15)
str_mayrhofer_acidic = ...                                                    % Creating a string element for the annotation
    sprintf("$R^{2}$ = %.5f", round(Mayrhofer_gof_acidic.rsquare, 5));
annotation('textbox', [.15 .8 .1 .1], 'String',str_mayrhofer_acidic,...       % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
ax_mayrhofer_acidic = gca; % current axes                                     % Creating an ax with gca such that the fontsize can be changed
ax_mayrhofer_acidic.TickDir = 'out';
box off
ax_mayrhofer_acidic.XAxis.FontSize = 15;                                      % Changing the tick size on the x-axis
ax_mayrhofer_acidic.YAxis.FontSize = 15;                                      % Changing the tick size on the y-axis


xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize',15)   % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                % Creating y-label
    'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.15 .80 .1 .1], 'String',["Mayrhofer -", "Acidic"],... % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);

xlim([Mayrhofer_E_acidic(1) Mayrhofer_E_acidic(end)])
ylim([Mayrhofer_i_acidic(1)*0 Mayrhofer_i_acidic(end)])



%% %%%%%%%%%%% The data from the Mayrhofer article %%%%%%%%%%%%%%%%%%%%%%
% These data is based on the highest anodic peak

Mayrhofer_potential_data = readmatrix("Mayrhofer_potential.xlsx");          % Mayrhofer potential vs time data - Don't really use this... computes it based on scan rate and initial potential
Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_2.xlsx");    % Mayrhofer dissolution vs time data - [ng/cm^2s]

Mayrhofer_dissolution = Mayrhofer_dissolution_data(5:end,2);                % Mayrhofer dissolution data - [ng/cm^2*s] -- Starting from 5 to remove the tail
Mayrhofer_time = Mayrhofer_dissolution_data(5:end,1);                       % Mayrhofer time data [s] -- Starting from 5 to remove the tail to be consistent

Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s

sweep_rate = 10*10^(-3);                                                    % Mayrhofer Sweep rate [V/s]
Mayrhofer_a_H_plus = 0.1*2;                                                 % Concentration of H+ (0.1 M H2SO4)
Mayrhofer_T = 25 + 273.13;                                                  % mayrhofer states room temperature

%% %%%%%%%%%%%%%%% Calling the diff equation solver %%%%%%%%%%%%%%%%%%%%%%

[t_scohy, theta_scohy] = diff_equation_solver(Mayrhofer_time, "value", Scohy_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj, theta_damj] = diff_equation_solver(Mayrhofer_time, "value", Damjanovic_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj_log, theta_damj_log] = diff_equation_solver(Mayrhofer_time, "value", Damjanovic_log_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_cherevko_acidic, theta_cherevko_acidic] = diff_equation_solver(Mayrhofer_time, "value", Cherevko_curve_acidic, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);


%%  Making plots of the ode15s solution and the interpolation 

% Transforming time to potential for the ode15s solution
potential_scohy_ode15s = CV_potential(t_scohy, "array");
potential_damj_ode15s = CV_potential(t_damj, "array");
potential_damj_log_ode15s = CV_potential(t_damj_log, "array");
potential_cherevko_ode15s_acidic = CV_potential(t_cherevko_acidic, "array");

% Interpolating the solution from ode15s to find values corresponding to
theta_interpol_scohy = interp1(t_scohy,theta_scohy,Mayrhofer_time);         % vq contains the interpolated values for theta from the
theta_interpol_damj = interp1(t_damj,theta_damj,Mayrhofer_time);            % solution of ode15s based on a value for k4
theta_interpol_damj_log = interp1(t_damj_log,theta_damj_log,Mayrhofer_time);% They do all contain equally many values as Mayrhofer_time has entries
theta_interpol_cherevko_acidic = interp1(t_cherevko_acidic,theta_cherevko_acidic,Mayrhofer_time);% They do all contain equally many values as Mayrhofer_time has entries


% Transforming time to potential fot the interpolation
potential_interpol = CV_potential(Mayrhofer_time, "array");


%%
%%%%%%%%%%%%%%%%%%%%%%%% Jobber med å tilpasse plots
% Voltammogram plots
figure()
plot(potential_scohy_ode15s, theta_scohy, "Color", "blue")
hold on
scatter(potential_interpol, theta_interpol_scohy, 20, "blue", 'o')
plot(potential_damj_ode15s, theta_damj, "Color", "green")
scatter(potential_interpol, theta_interpol_damj, 20, "green", '+')
%plot(potential_damj_log_ode15s, theta_damj_log, "Color", "red")
%scatter(potential_interpol, theta_interpol_damj_log, 20, "red", 'x')
plot(potential_cherevko_ode15s_acidic, theta_cherevko_acidic, "Color", "red")
scatter(potential_interpol, theta_interpol_cherevko_acidic, 20, "red", 'x')
legend(["Scohy fit", "Scohy fit interpol", "Damjanovic fit", "Damjanovic fit interpol","Cherevko-fit", "Cherevko fit interpol"], Location = "best")
xlabel('Potential - E/[$V$]','Interpreter','latex')
ylabel('$\theta_{2}(E)$ - [$-$]','Interpreter','latex')
%title("Voltammogram of $\theta_{2}$",'Interpreter','latex')

% Regular plots
figure()
plot(t_scohy, theta_scohy, "Color", "blue")
hold on
scatter(Mayrhofer_time, theta_interpol_scohy, 20,"blue", 'o')
plot(t_damj, theta_damj, "Color", "green")
scatter(Mayrhofer_time, theta_interpol_damj, 20,"green", '+')
%plot(t_damj_log, theta_damj_log, "Color", "red")
%scatter(Mayrhofer_time, theta_interpol_damj_log, 20, "red", 'x')
plot(t_cherevko_acidic, theta_cherevko_acidic, "Color", "red")
scatter(Mayrhofer_time, theta_interpol_cherevko_acidic, 20, "red", 'x')
hold off
xlabel('Time -t [$s$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
ylim([0 max(theta_scohy)*1.1])
%title("$\theta_{2}$ as function of t",'Interpreter','latex')

yyaxis right
plot(Mayrhofer_time, potential_interpol, 'Marker', 'o')
ylabel('E - [$V vs RHE$]','Interpreter','latex')
ylim([0.4 max(potential_interpol)*1.1])
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Cherevko fit ode15s", "Cherevko fit interpol", "potential regime"], Location = "best")


%% Plots of degradation rates and the rate and the solution from the solver

figure()
plot(Mayrhofer_time, Mayrhofer_dissolution_mole, 'Marker', 'o')
ylabel('$\frac{d Ir}{dt}$ - [$\frac{mol}{cm^{2}s}$]','Interpreter','latex')
xlabel('time - [$s$]','Interpreter','latex')
title("Rate",'Interpreter','latex')

yyaxis right
plot(Mayrhofer_time, potential_interpol, 'Marker', 'o')
ylabel('E - [$V vs RHE$]','Interpreter','latex')
legend(["r_{diss}", "potential regime"], Location = "best")
%legend('r_{diss}', Location = 'best')


% Theta and dissolution rate as a function of time
figure()
%cla reset
title("$\theta_{2}(t)$ and $\frac{d Ir}{d t}(t)$",'Interpreter','latex')
plot(t_scohy, theta_scohy, "Color", "blue")
hold on
scatter(Mayrhofer_time, theta_interpol_scohy ,20,"blue" ,'o')
plot(t_damj, theta_damj,"Color","green")
scatter(Mayrhofer_time, theta_interpol_damj, 20,"green" ,'+')
%plot(t_damj_log, theta_damj_log, "Color", "red")
%scatter(Mayrhofer_time, theta_interpol_damj_log, 20,"red" ,'x')
plot(t_cherevko_acidic, theta_cherevko_acidic, "Color", "red")
scatter(Mayrhofer_time, theta_interpol_cherevko_acidic, 20,"red" ,'x')
hold off
xlabel('time - [$s$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')

yyaxis right
plot(Mayrhofer_time, Mayrhofer_dissolution_mole, 'Marker', 'o')
ylabel('$\frac{d Ir}{dt}$ - [$\frac{mol}{cm^{2}s}$]','Interpreter','latex')
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Cherevko ode15s", "Cherevko fit interpol", "r_{diss}"], Location = "best")
%legend('r_{diss}', Location = 'best')

title("Theta and dissolution rate as a function of potential")
figure()
%cla reset
title("$\theta_{2}(E)$ and $\frac{d Ir}{d t}(E)$",'Interpreter','latex')
plot(potential_scohy_ode15s, theta_scohy, "Color", "blue")
hold on
scatter(potential_interpol, theta_interpol_scohy, 20, "blue", 'o')
plot(potential_damj_ode15s, theta_damj, "Color","green")
scatter(potential_interpol, theta_interpol_damj, 20, "green", '+')
%plot(potential_damj_log_ode15s, theta_damj_log, "Color", "red")
%scatter(potential_interpol, theta_interpol_damj_log, 20, "red", 'x')
plot(potential_cherevko_ode15s_acidic, theta_cherevko_acidic, "Color", "red")
scatter(potential_interpol, theta_interpol_cherevko_acidic, 20, "red", 'x')
hold off
xlabel('Potential - E [$V$]','Interpreter','latex')
ylabel('$\theta_{2}(E)$ - [$-$]','Interpreter','latex')
title("Voltammogram of $\theta_{2}$",'Interpreter','latex')

yyaxis right
plot(potential_interpol, Mayrhofer_dissolution_mole, 'Marker', 'o')
ylabel('$\frac{d Ir}{dt}$ - [$\frac{mol}{cm^{2}s}$]','Interpreter','latex')
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Cherevko fit ode15s", "Cherevko fit interpol", "r_{diss}"], Location = "best")

%%   
% Normalisation just wouldn't do

fun_scohy = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_scohy,theta_scohy,x);
FT_scohy = fittype(fun_scohy, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

fun_damj = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_damj,theta_damj,x);
FT_damj = fittype(fun_damj, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

fun_damj_log = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_damj_log,theta_damj_log,x);
FT_damj_log = fittype(fun_damj_log, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

fun_cherevko = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_cherevko_acidic,theta_cherevko_acidic,x);
FT_cherevko = fittype(fun_cherevko, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', k_4_0_plus, ...
           'StartPoint', 1e-4,...
           'TolFun', 1e-22);                                                % k_3_0_plus
          

[curve_scohy_3] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT_scohy,FO);    
[curve_damj_3] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT_damj,FO);    
[curve_damj_log_3] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT_damj_log,FO);    
[curve_cherevko_3] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT_cherevko,FO);    

figure()
plot(curve_cherevko_3)
hold on
plot(Mayrhofer_time, Mayrhofer_dissolution_mole)
xlabel("Time")
ylabel("r_{3}")