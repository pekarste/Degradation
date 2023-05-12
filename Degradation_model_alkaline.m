%% %%%%%%%%%%%%%%%%%%%%%% ALKALINE model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will couple all the different functions together and work
% like a masterscript. 

%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_n = 0.40;                                                                % V
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Guess for k_4_0_plus
k_4_0_plus = 10^-2;

% Data used for fitting r2_alkaline
Cherevko_alkaline = readmatrix("Data\Alkaline\Cherevko_alkaline_polarisation_data.xlsx");% Potential/current density data from Cherevko
Damjanovic_alkaline = readmatrix("Data\Alkaline\Damjanovic_alkaline_polarisation.xlsx"); % Current density/potential data from Damjanovic

%% Extracted data from the Excel files 

% Alkaline - Cherevko
Cherevko_E_alkaline = Cherevko_alkaline(1:end,1);                           % [V vs RHE] - Potential
Cherevko_i_alkaline = Cherevko_alkaline(1:end,2)*10^(-3+4);                 % [A/m^2] - Current density, originally in mA/cm^2
Cherevko_T_alkaline = 25 + 273;                                             % [K] - Temperature
Cherevko_OH_alkaline = 0.05*2;                                              % [-] - Activity of H+


% Alkaline - Damjanovic
Damjanovic_E_alkaline = Damjanovic_alkaline(1:end,2);                       % [V vs RHE]          
Damjanovic_i_alkaline = Damjanovic_alkaline(1:end,1)*10^4;                  % [A/m^2] - Originally A/cm^2 (The article contains log(i), but I converted thm
Damjanovic_T_alkaline = 25 + 273;                                           % [K] - Temperature
Damjanovic_OH_alkaline = 1;                                                 % [-] - Activity of H+

%% %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting the expression of the current based on r_2 to the data
% The r_2_fit returns the curve (fitting results) and the gof.
% The coefficients are contained in the curve

% Scohy
[Cherevko_curve_alkaline, Cherevko_gof_alkaline] = ...                                        % This is the expression with rds
    r_2_fit_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear");

% Damjanovic
[Damjanovic_curve_alkaline, Damjanovic_gof_alkaline] = ...                                    % This is the expression with rds
    r_2_fit_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear");

% Damjanovic (log)
[Damjanovic_log_curve_alkaline, Damjanovic_log_gof_alkaline] = ...                            % This is the expression with rds
    r_2_fit_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic");

%% Plotting the curve fittings

figure("Name","Cherevko Alkaline Fitting")                                                                    % Creates figure
scatter(Cherevko_E_alkaline, Cherevko_i_alkaline, 45, "filled", "blue", 'o')  % Scatter plot of the sampled values from Scohy
hold on
fig_cherevko_fit_alkaline = plot(Cherevko_curve_alkaline, "black");                                 % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_cherevko_fit_alkaline,'lineWidth',1);                                           % Changing the linewidth of the curve of the cfit
hold off
ax_cherevko_alkaline = gca; % current axes                                              % Creating an ax with gca such that the fontsize can be changed
legend({'Data', 'Fitting'},...                                              % Creating a legend for the graphs
    'Position', [.2 .65 .1 .1], 'Interpreter','latex', 'FontSize',15)
str_cherevko_alkaline = sprintf("$R^{2}$ = %.5f", round(Cherevko_gof_alkaline.rsquare, 5));         % Creating a string element for the annotation
annotation('textbox', [.15 .8 .1 .1], 'String',str_cherevko_alkaline,...                % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
ax_cherevko_alkaline.XAxis.FontSize = 12;                                               % Changing the tick size on the x-axis
ax_cherevko_alkaline.YAxis.FontSize = 12;                                               % Changing the tick size on the y-axis
xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize', 15)  % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                % Creating y-label
    'Interpreter','latex', 'FontSize', 15)                                  

figure("Name","Damjanovic Alkaline Fitting")                                                                    % Creates figure
scatter(Damjanovic_E_alkaline, Damjanovic_i_alkaline, 45,...            % Scatter plot of the sampled values from Damjanovic
    "filled", "green", "square" )
hold on
fig_damjanovic_alkaline = plot(Damjanovic_curve_alkaline, "black");                           % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_damjanovic_alkaline,'lineWidth',1);                                          % Changing the linewidth of the curve of the cfit
hold off
ax_damjanovic_alkaline = gca; % current axes                                         % Creating an ax with gca such that the fontsize can be changed
legend({'Data', 'Fitting'},...                                              % Creating a legend for the graphs
    'Position', [.2 .65 .1 .1], 'Interpreter','latex', 'FontSize',15)
str_damjanovic_alkaline = ...                                                        % Creating a string element for the annotation
    sprintf("$R^{2}$ = %.5f", round(Damjanovic_gof_alkaline.rsquare, 5));
annotation('textbox', [.15 .8 .1 .1], 'String',str_damjanovic_alkaline,...           % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex','FitBoxToText','on', 'FontSize',15);
ax_damjanovic_alkaline.XAxis.FontSize = 12;                                          % Changing the tick size on the x-axis
ax_damjanovic_alkaline.YAxis.FontSize = 12;                                          % Changing the tick size on the y-axis
xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize',15)   % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                % Creating y-label
    'Interpreter','latex', 'FontSize',15)

figure("Name", "Damjanovic Alkaline Fitting Log")                                                                    % Creating figure
%plot(Damjanovic_log_curve_fit)
scatter(Damjanovic_E_alkaline, log10(Damjanovic_i_alkaline), 45,...     % Scatter plot of the sampled values from Damjanovic
    "filled", "red", "^")
hold on
fig_damjanovic_log_alkaline = plot(Damjanovic_log_curve_alkaline, "black");                   % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_damjanovic_log_alkaline,'lineWidth',1);                                      % Changing the linewidth of the curve of the cfit
hold off
ax_damjanovic_log_alkaline = gca; % current axes                                     % Creating an ax with gca such that the fontsize can be changed
legend({'Data', 'Fitting'}, 'Position', [.2 .65 .1 .1],...                  % Creating a legend for the graphs
    'Interpreter','latex', 'FontSize',15)
str_damjanovic_log_alkaline = ...                                                    % Creating a string element for the annotation
    sprintf("$R^{2}$ = %.5f", round(Damjanovic_log_gof_alkaline.rsquare, 5));
annotation('textbox', [.15 .8 .1 .1], 'String',str_damjanovic_log_alkaline,...       % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
ax_damjanovic_log_alkaline.XAxis.FontSize = 12;                                      % Changing the tick size on the x-axis
ax_damjanovic_log_alkaline.YAxis.FontSize = 12;                                      % Changing the tick size on the y-axis
xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize',15)   % Creating x-label
ylabel('$\log_{10}$ of current density - $\log{i}$/[$Am^{-2}$]',...         % Creating y-label
    'Interpreter','latex', 'FontSize',15)

%% %%%%%%%%%%% The data from the Mayrhofer article %%%%%%%%%%%%%%%%%%%%%%
% These data is based on the highest anodic peak

Mayrhofer_potential_data = readmatrix("Mayrhofer_potential.xlsx");          % Mayrhofer potential vs time data - Don't really use this... computes it based on scan rate and initial potential
Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_2.xlsx");    % Mayrhofer dissolution vs time data - [ng/cm^2s]

Mayrhofer_dissolution = Mayrhofer_dissolution_data(5:end,2);                % Mayrhofer dissolution data - [ng/cm^2*s] -- Starting from 5 to remove the tail
Mayrhofer_time = Mayrhofer_dissolution_data(5:end,1);                       % Mayrhofer time data [s] -- Starting from 5 to remove the tail to be consistent

Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s

sweep_rate = 10*10^(-3);                                                    % Mayrhofer Sweep rate [V/s]
Mayrhofer_a_H_plus = 10^(-14)/(0.1*2);                                      % Concentration of H+ (0.1 M H2SO4)
Mayrhofer_T = 25 + 273.13;                                                  % mayrhofer states room temperature

%% %%%%%%%%%%%%%%% Calling the diff equation solver %%%%%%%%%%%%%%%%%%%%%%

[t_cherevko_alkaline, theta_cherevko_alkaline] = diff_equation_solver_alkaline(Mayrhofer_time, "value", Cherevko_curve_alkaline, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj_alkaline, theta_damj_alkaline] = diff_equation_solver_alkaline(Mayrhofer_time, "value", Damjanovic_curve_alkaline, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj_log_alkaline, theta_damj_log_alkaline] = diff_equation_solver_alkaline(Mayrhofer_time, "value", Damjanovic_log_curve_alkaline, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);


%%  Making plots of the ode15s solution and the interpolation 

% Transforming time to potential for the ode15s solution
potential_cherevko_alkaline_ode15s = CV_potential(t_cherevko_alkaline, "array");
potential_damj_alkaline_ode15s = CV_potential(t_damj_alkaline, "array");
potential_damj_log_alkaline_ode15s = CV_potential(t_damj_log_alkaline, "array");

% Interpolating the solution from ode15s to find values corresponding to
theta_interpol_cherevko_alkaline = interp1(t_cherevko_alkaline,theta_cherevko_alkaline,Mayrhofer_time);         % vq contains the interpolated values for theta from the
theta_interpol_damj_alkaline = interp1(t_damj_alkaline,theta_damj_alkaline,Mayrhofer_time);            % solution of ode15s based on a value for k4
theta_interpol_damj_log_alkaline = interp1(t_damj_log_alkaline,theta_damj_log_alkaline,Mayrhofer_time);% They do all contain equally many values as Mayrhofer_time has entries

% Transforming time to potential fot the interpolation
potential_interpol = CV_potential(Mayrhofer_time, "array");


%%
%%%%%%%%%%%%%%%%%%%%%%%% Jobber med å tilpasse plots
% Voltammogram plots
figure()
plot(potential_cherevko_alkaline_ode15s, theta_cherevko_alkaline, "Color", "blue")
hold on
scatter(potential_interpol, theta_interpol_cherevko_alkaline, 20, "blue", 'o')
plot(potential_damj_alkaline_ode15s, theta_damj_alkaline, "Color", "green")
scatter(potential_interpol, theta_interpol_damj_alkaline, 20, "green", '+')
plot(potential_damj_log_alkaline_ode15s, theta_damj_log_alkaline, "Color", "red")
scatter(potential_interpol, theta_interpol_damj_log_alkaline, 20, "red", 'x')
legend(["Scohy fit", "Scohy fit interpol", "Damjanovic fit", "Damjanovic fit interpol","Damjanovic log-fit", "Damjanovic log fit interpol"], Location = "best")
xlabel('Potential - E/[$V$]','Interpreter','latex')
ylabel('$\theta_{2}(E)$ - [$-$]','Interpreter','latex')
%title("Voltammogram of $\theta_{2}$",'Interpreter','latex')

% Regular plots
figure()
plot(t_cherevko_alkaline, theta_cherevko_alkaline, "Color", "blue")
hold on
scatter(Mayrhofer_time, theta_interpol_cherevko_alkaline, 20,"blue", 'o')
plot(t_damj_alkaline, theta_damj_alkaline, "Color", "green")
scatter(Mayrhofer_time, theta_interpol_damj_alkaline, 20,"green", '+')
plot(t_damj_log_alkaline, theta_damj_log_alkaline, "Color", "red")
scatter(Mayrhofer_time, theta_interpol_damj_log_alkaline, 20, "red", 'x')
hold off
xlabel('Time -t [$s$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
ylim([0 max(theta_cherevko_alkaline)*1.1])
%title("$\theta_{2}$ as function of t",'Interpreter','latex')

yyaxis right
plot(Mayrhofer_time, potential_interpol, 'Marker', 'o')
ylabel('E - [$V vs RHE$]','Interpreter','latex')
ylim([0.4 max(potential_interpol)*1.1])
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Damjanovic log fit ode15s", "Damjanovic log fit interpol", "potential regime"], Location = "best")


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
%yyaxis left
plot(t_cherevko_alkaline, theta_cherevko_alkaline, "Color", "blue")
hold on
scatter(Mayrhofer_time, theta_interpol_cherevko_alkaline ,20,"blue" ,'o')
plot(t_damj_alkaline, theta_damj_alkaline,"Color","green")
scatter(Mayrhofer_time, theta_interpol_damj_alkaline, 20,"green" ,'+')
plot(t_damj_log_alkaline, theta_damj_log_alkaline, "Color", "red")
scatter(Mayrhofer_time, theta_interpol_damj_log_alkaline, 20,"red" ,'x')
hold off
xlabel('time - [$s$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')

yyaxis right
plot(Mayrhofer_time, Mayrhofer_dissolution_mole, 'Marker', 'o')
ylabel('$\frac{d Ir}{dt}$ - [$\frac{mol}{cm^{2}s}$]','Interpreter','latex')
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Damjanovic log fit ode15s", "Damjanovic log fit interpol", "r_{diss}"], Location = "best")
%legend('r_{diss}', Location = 'best')

title("Theta and dissolution rate as a function of potential")
figure()
%cla reset
title("$\theta_{2}(E)$ and $\frac{d Ir}{d t}(E)$",'Interpreter','latex')
plot(potential_cherevko_alkaline_ode15s, theta_cherevko_alkaline, "Color", "blue")
hold on
scatter(potential_interpol, theta_interpol_cherevko_alkaline, 20, "blue", 'o')
plot(potential_damj_alkaline_ode15s, theta_damj_alkaline, "Color","green")
scatter(potential_interpol, theta_interpol_damj_alkaline, 20, "green", '+')
plot(potential_damj_log_alkaline_ode15s, theta_damj_log_alkaline, "Color", "red")
scatter(potential_interpol, theta_interpol_damj_log_alkaline, 20, "red", 'x')
hold off
xlabel('Potential - E [$V$]','Interpreter','latex')
ylabel('$\theta_{2}(E)$ - [$-$]','Interpreter','latex')
title("Voltammogram of $\theta_{2}$",'Interpreter','latex')

yyaxis right
plot(potential_interpol, Mayrhofer_dissolution_mole, 'Marker', 'o')
ylabel('$\frac{d Ir}{dt}$ - [$\frac{mol}{cm^{2}s}$]','Interpreter','latex')
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Damjanovic log fit ode15s", "Damjanovic log fit interpol", "r_{diss}"], Location = "best")

%%   
% Normalisation just wouldn't do

fun = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_cherevko_alkaline,theta_cherevko_alkaline,x);
FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', k_4_0_plus, ...
           'StartPoint', 1e-4,...
           'TolFun', 1e-22);                                                % k_3_0_plus
          

[curve, gof, output,warnstr,errstr,convmsg] = fit(Mayrhofer_time,Mayrhofer_dissolution_mole,FT,FO);    

figure()
plot(curve)
hold on
plot(Mayrhofer_time, Mayrhofer_dissolution_mole)
xlabel("Time")
ylabel("r_{3}")