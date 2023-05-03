%% %%%%%%%%%%%%%%%%%%%%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will couple all the different functions together and work
% like a masterscript. 

%% Define Physical Constants

R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
E_n = 1.229;                                                                % V
a_H2O = 1;                                                                  % [-]
Mm_Ir = 192.2;                                                              % g/mol [SI]
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]

%% %%%%%%%%%%%%%%%% DATA for the fitting %%%%%%%%%%%%%%%%%%%% 

% Guess for k_4_0_plus
k_4_0_plus = 10^-2;

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

figure()
plot(Scohy_curve_fit)
hold on
plot(Scohy_curve)
plot(Scohy_potential, Scohy_current_density, "Color","blue")
hold off
legend({'Non_rds', 'rds', 'data'})
xlabel('E')
ylabel('i')

figure()
plot(Damjanovic_curve_fit)
hold on
plot(Damjanovic_curve)
plot(Damjanovic_potential, Damjanovic_current_density, "Color","blue")
hold off
legend({'Non_rds', 'rds', 'data'})
xlabel('E')
ylabel('i')

figure()
plot(Damjanovic_log_curve_fit)
hold on
plot(Damjanovic_log_curve)
plot(Damjanovic_potential, log10(Damjanovic_current_density), "Color","blue")
hold off
legend({'Non_rds', 'rds', 'data'})
xlabel('E')
ylabel('i')

%% %%%%%%%%%%% The data from the Mayrhofer article %%%%%%%%%%%%%%%%%%%%%%
% These data is based on the highest anodic peak

Mayrhofer_potential_data = readmatrix("Mayrhofer_potential.xlsx");          % Mayrhofer potential vs time data - Don't really use this... computes it based on scan rate and initial potential
Mayrhofer_dissolution_data = readmatrix("Mayrhofer_dissolution_2.xlsx");    % Mayrhofer dissolution vs time data - [ng/cm^2s]

Mayrhofer_dissolution = Mayrhofer_dissolution_data(1:end,2);                % Mayrhofer dissolution data - [ng/cm^2*s] 
Mayrhofer_time = Mayrhofer_dissolution_data(1:end,1);                       % Mayrhofer time data [s] - 

Mayrhofer_dissolution_mole = Mayrhofer_dissolution*10^(-9)*10^(4)/Mm_Ir;    % Changes the units from ng/cm^2*s --> mole/m^2*s

sweep_rate = 10*10^(-3);                                                    % Mayrhofer Sweep rate [V/s]
Mayrhofer_a_H_plus = 0.1*2;                                                 % Concentration of H+ (0.1 M H2SO4)
Mayrhofer_T = 25 + 273.13;                                                  % mayrhofer states room temperature

%% %%%%%%%%%%%%%%% Calling the diff equation solver %%%%%%%%%%%%%%%%%%%%%%

[t_scohy, theta_scohy] = diff_equation_solver(Mayrhofer_time, "value", Scohy_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj, theta_damj] = diff_equation_solver(Mayrhofer_time, "value", Damjanovic_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);
[t_damj_log, theta_damj_log] = diff_equation_solver(Mayrhofer_time, "value", Damjanovic_log_curve, Mayrhofer_a_H_plus, Mayrhofer_T, k_4_0_plus, eps);


%%  Making plots of the ode15s solution and the interpolation 

% Transforming time to potential for the ode15s solution
potential_scohy_ode15s = CV_potential(t_scohy, "array");
potential_damj_ode15s = CV_potential(t_damj, "array");
potential_damj_log_ode15s = CV_potential(t_damj_log, "array");

% Interpolating the solution from ode15s to find values corresponding to
theta_interpol_scohy = interp1(t_scohy,theta_scohy,Mayrhofer_time);         % vq contains the interpolated values for theta from the
theta_interpol_damj = interp1(t_damj,theta_damj,Mayrhofer_time);            % solution of ode15s based on a value for k4
theta_interpol_damj_log = interp1(t_damj_log,theta_damj_log,Mayrhofer_time);% They do all contain equally many values as Mayrhofer_time has entries

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
plot(potential_damj_log_ode15s, theta_damj_log, "Color", "red")
scatter(potential_interpol, theta_interpol_damj_log, 20, "red", 'x')
legend(["Scohy fit", "Scohy fit interpol", "Damjanovic fit", "Damjanovic fit interpol","Damjanovic log-fit", "Damjanovic log fit interpol"], Location = "best")
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
plot(t_damj_log, theta_damj_log, "Color", "red")
scatter(Mayrhofer_time, theta_interpol_damj_log, 20, "red", 'x')
hold off
xlabel('Time -t [$s$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
ylim([0 max(theta_scohy)*1.1])
%title("$\theta_{2}$ as function of t",'Interpreter','latex')

yyaxis right
plot(Mayrhofer_time, potential_interpol, 'Marker', 'o')
ylabel('E - [$V vs RHE$]','Interpreter','latex')
ylim([0.4 max(potential_interpol)*1.1])
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Damjanovic log fit ode15s", "Damjanovic log fit interpol", "potential regime"], Location = "best")


% % Regular plots
% figure()
% plot(Mayrhofer_time+dt*0, theta_interpol_scohy,'Marker','o')
% hold on
% plot(Mayrhofer_time+dt*0, theta_interpol_damj, 'Marker','+')
% plot(Mayrhofer_time+dt*0, theta_interpol_damj_log, 'Marker','x')
% legend(["Scohy fit", "Damjanovic fit", "Damjanovic log fit"], Location = "best")
% xlabel('Time -t [$s$]','Interpreter','latex')
% ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
% title("$\theta$ as function of t (interpol)",'Interpreter','latex')
% 
% % Voltammogram plots
% figure()
% plot(potential_interpol, theta_interpol_scohy)
% hold on
% plot(potential_interpol, theta_interpol_damj)
% plot(potential_interpol, theta_interpol_damj_log)
% legend(["Scohy fit", "Damjanovic fit", "Damjanovic log fit"], Location = "best")
% xlabel('Potential - E [$V$]','Interpreter','latex')
% ylabel('$\theta_{2}(E)$ - [$-$]','Interpreter','latex')
% title("Voltammogram of $\theta_{2}$ (interpol)",'Interpreter','latex')



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
plot(t_scohy, theta_scohy, "Color", "blue")
hold on
scatter(Mayrhofer_time, theta_interpol_scohy ,20,"blue" ,'o')
plot(t_damj, theta_damj,"Color","green")
scatter(Mayrhofer_time, theta_interpol_damj, 20,"green" ,'+')
plot(t_damj_log, theta_damj_log, "Color", "red")
scatter(Mayrhofer_time, theta_interpol_damj_log, 20,"red" ,'x')
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
plot(potential_scohy_ode15s, theta_scohy, "Color", "blue")
hold on
scatter(potential_interpol, theta_interpol_scohy, 20, "blue", 'o')
plot(potential_damj_ode15s, theta_damj, "Color","green")
scatter(potential_interpol, theta_interpol_damj, 20, "green", '+')
plot(potential_damj_log_ode15s, theta_damj_log, "Color", "red")
scatter(potential_interpol, theta_interpol_damj_log, 20, "red", 'x')
hold off
xlabel('Potential - E [$V$]','Interpreter','latex')
ylabel('$\theta_{2}(E)$ - [$-$]','Interpreter','latex')
title("Voltammogram of $\theta_{2}$",'Interpreter','latex')

yyaxis right
plot(potential_interpol, Mayrhofer_dissolution_mole, 'Marker', 'o')
ylabel('$\frac{d Ir}{dt}$ - [$\frac{mol}{cm^{2}s}$]','Interpreter','latex')
legend(["Scohy fit ode15s", "Scohy fit interpol", "Damjanovic fit ode15s", "Damjanovic fit interpol","Damjanovic log fit ode15s", "Damjanovic log fit interpol", "r_{diss}"], Location = "best")

%%   

fun = @(k_3_0_plus, x) gamma*k_3_0_plus*a_H2O*interp1(t_scohy,theta_scohy,x);
FT = fittype(fun, 'independent',{'x'}, 'coefficients',{'k_3_0_plus'});

FO = fitoptions('Method','NonLinearLeastSquares',...
           'Lower', eps,...                                                 % k_3_0_plus
           'Upper', 10^(2), ...
           'StartPoint', 1e-4,...
           'TolFun', 1e-20);                                              % k_3_0_plus
          


[curve, gof, output,warnstr,errstr,convmsg] = fit(Mayrhofer_time(5:end),Mayrhofer_dissolution_mole(5:end),FT,FO);    

figure()
plot(curve)
hold on
plot(Mayrhofer_time, Mayrhofer_dissolution_mole)
xlabel("Time")
ylabel("r_{3}")
%% Extracting the k_3_0_+ for a given value of k_4_0_+

% Fitting the calculated values to the measured values
% [scohy_degradation_curve, scohy_degradation_gof] = chi_square(theta_interpol_scohy(5:end), Mayrhofer_dissolution_mole(5:end));                     % The fitting from Scohy solution with the measured data
% [damjanovic_degradation_curve, damjanovic_degradation_gof] = chi_square(theta_interpol_damj(5:end), Mayrhofer_dissolution_mole(5:end));            % The fitting from Damjanovic solution with the measured data
% [damjanovic_log_degradation_curve, damjanovic_log_degradation_gof] = chi_square(theta_interpol_damj_log(5:end), Mayrhofer_dissolution_mole(5:end));% The fitting from Damjanovic log solution woth the measured data
% 
% % Extracting k_3_0_plus from the fitting
% scohy_k_3_0_plus = scohy_degradation_curve.k_3_0_plus;                            % k_3_0_plus from scohy solution
% damj_k_3_0_plus = damjanovic_degradation_curve.k_3_0_plus;                        % k_3_0_plus from damj solution
% damj_log_k_3_0_plus = damjanovic_log_degradation_curve.k_3_0_plus;                % k_3_0_plus from damj log solution
% 
% % calculating the degradation rate with the value of k_3_0_plus
% scohy_deg_rate = gamma*theta_interpol_scohy*a_H2O*scohy_k_3_0_plus;                  
% damj_deg_rate = gamma*theta_interpol_damj*a_H2O*damj_k_3_0_plus;
% damj_log_deg_rate = gamma*theta_interpol_damj_log*a_H2O*damj_log_k_3_0_plus;
% 
% figure()
% plot(theta_scohy, scohy_k_3_0_plus*gamma*a_H2O*theta_scohy)
% xlabel("Theta")
% ylabel("r_{3}")
% 
% figure()
% plot(t_scohy, scohy_k_3_0_plus*gamma*a_H2O*theta_scohy)
% xlabel("t")
% ylabel("r_{3}")
% 
% figure()
% plot(t_scohy, theta_scohy)
% xlabel("t")
% ylabel("theta_{2}")
% 
% % Plotting the new degradation rate
% figure()
% title("$\frac{d Ir}{d t}(t)$ vs time",'Interpreter','latex')
% plot(Mayrhofer_time, scohy_deg_rate, '-o', "Color", "blue")
% hold on
% %scatter(Mayrhofer_time, theta_interpol_scohy ,20,"blue" ,'o')
% plot(Mayrhofer_time, damj_deg_rate, '-x' ,"Color","green")
% %scatter(Mayrhofer_time, theta_interpol_damj, 20,"green" ,'+')
% plot(Mayrhofer_time, damj_log_deg_rate, '-+', "Color", "red")
% %scatter(Mayrhofer_time, theta_interpol_damj_log, 20,"red" ,'x')
% plot(Mayrhofer_time, Mayrhofer_dissolution_mole, '--', "Color", "magenta")
% hold off
% xlabel('time - [$s$]','Interpreter','latex')
% ylabel('$\frac{d Ir}{d t}(t)$ - [$-$]','Interpreter','latex')