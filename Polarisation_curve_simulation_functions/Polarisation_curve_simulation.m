% This script will be used to try and simulate the polarisation curves
% based on the expression from Reksten and perhaps the one from Marshall as
% well. This is becuase we wont do another fitting, but rather find some
% values for k_4 which are within reasonable values

% This was originally acidic, but because of the recent descovery of the
% alkaline model, I will turn this alkaline as well


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
k_4_0_plus = 10^(-2);

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
Schalenbach_a_OH_alkaline = 0.05*1;                                              % [-] - Activity of OH-
Schalenbach_sweep_rate = 2*10^(-3);                                            % [V/s] -Schalenbach Sweep rate 
%--------------------------------------------------------------------------

%% %%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting the expression of the current based on r_2 to the data
% The r_2_fit returns the curve (fitting results) and the gof.
% The coefficients are contained in the curve

% Cherevko
[Cherevko_curve_alkaline, Cherevko_gof_alkaline] = ...                                        % This is the expression with rds
    r_2_fit_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear");

% Damjanovic
[Damjanovic_curve_alkaline, Damjanovic_gof_alkaline] = ...                                    % This is the expression with rds
    r_2_fit_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear");

% Damjanovic (log)
[Damjanovic_log_curve_alkaline, Damjanovic_log_gof_alkaline] = ...                            % This is the expression with rds
    r_2_fit_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic");

%--------------------------------------------------------------------------
% Schalenbach - Alkaline 
[Schalenbach_curve_alkaline, Schalenbach_gof_alkaline] = ...                            % This is the expression with rds
    r_2_fit_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_a_OH_alkaline, Schalenbach_T_alkaline, "Linear");
%--------------------------------------------------------------------------

%% %%%%%%%%%%%%%%%%%%%%%%%% REKSTEN %%%%%%%%%%%%%%%%%%%%% %%

% Introducing the fitting from the large equation by using the separate
% fits from earlier

[Cherevko_curve_alkaline_reksten, Cherevko_gof_alkaline_reksten] = r_reksten_fit(Cherevko_curve_alkaline, Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear");
[Damjanovic_curve_alkaline_reksten, Damjanovic_gof_alkaline_reksten] = r_reksten_fit(Damjanovic_curve_alkaline, Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear");
[Damjanovic_log_curve_alkaline_reksten, Damjanovic_log_gof_alkaline_reksten] = r_reksten_fit(Damjanovic_log_curve_alkaline, Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic");
[Schalenbach_curve_alkaline_reksten, Schalenbach_gof_alkaline_reksten] = r_reksten_fit(Schalenbach_curve_alkaline, Schalenbach_E_alkaline,Schalenbach_i_alkaline, Schalenbach_a_OH_alkaline, Schalenbach_T_alkaline, "Linear");


%% Plotting the curves from the fitting

%--------------------------------------------------------------------------
% Cherevko - Alkaline
figure("Name","Cherevko Alkaline Fitting")                                                 % Creates figure
scatter(Cherevko_E_alkaline, Cherevko_i_alkaline, 45, "filled", "red", "square")           % Scatter plot of the sampled values from Scohy
hold on
fig_cherevko_fit_alkaline = plot(Cherevko_curve_alkaline_reksten, "black");                % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_cherevko_fit_alkaline,'lineWidth',1);                                              % Changing the linewidth of the curve of the cfit
hold off
ax_cherevko_alkaline = gca; % current axes                                                 % Creating an ax with gca such that the fontsize can be changed
legend({'Data', 'Fitting'},...                                                             % Creating a legend for the graphs
    'Position', [.2 .65 .1 .1], 'Interpreter','latex', 'FontSize',15)
str_cherevko_alkaline = sprintf("$R^{2}$ = %.5f", round(Cherevko_gof_alkaline_reksten.rsquare, 5));% Creating a string element for the annotation
annotation('textbox', [.15 .8 .1 .1], 'String',str_cherevko_alkaline,...                   % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
ax_cherevko_alkaline.XAxis.FontSize = 12;                                                  % Changing the tick size on the x-axis
ax_cherevko_alkaline.YAxis.FontSize = 12;                                                  % Changing the tick size on the y-axis
xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize', 15)                 % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                               % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Cherevko_E_alkaline(1) Cherevko_E_alkaline(end)])
ylim([Cherevko_i_alkaline(1)*0 Cherevko_i_alkaline(end)])
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Schalenbach - Alkaline
figure("Name","Schalenbach fitting")                                                % Creates figure
scatter(Schalenbach_E_alkaline, Schalenbach_i_alkaline, 45, "filled", "magenta", "diamond")  % Scatter plot of the sampled values from Scohy
hold on
fig_schalenbach_fit_alkaline = plot(Schalenbach_curve_alkaline_reksten, "black");           % Creating a fig to stor the plot of the curve fit (cfit element)
set(fig_schalenbach_fit_alkaline,'lineWidth',1);                                    % Changing the linewidth of the curve of the cfit
hold off
ax_schalenbach_alkaline = gca; % current axes                                       % Creating an ax with gca such that the fontsize can be changed
legend({'Data', 'Fitting'},...                                                      % Creating a legend for the graphs
    'Position', [.2 .65 .1 .1], 'Interpreter','latex', 'FontSize',15)
str_schalenbach_alkaline = sprintf("$R^{2}$ = %.5f",...
    round(Schalenbach_gof_alkaline_reksten.rsquare, 5));                                    % Creating a string element for the annotation
annotation('textbox', [.15 .8 .1 .1], 'String',str_schalenbach_alkaline,...         % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
ax_schalenbach_alkaline.XAxis.FontSize = 12;                                        % Changing the tick size on the x-axis
ax_schalenbach_alkaline.YAxis.FontSize = 12;                                        % Changing the tick size on the y-axis
xlabel('Potential - E/[$V$] vs RHE','Interpreter','latex', 'FontSize', 15)          % Creating x-label
ylabel('Current density - i/[$Am^{-2}$]',...                                        % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([Schalenbach_E_alkaline(1) Schalenbach_E_alkaline(end)])
ylim([Schalenbach_i_alkaline(1)*0 Schalenbach_i_alkaline(end)])
%--------------------------------------------------------------------------

%% %%%%%%%%%% Marshall %%%%%%%%%%%% %%