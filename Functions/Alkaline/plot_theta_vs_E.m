% This will be a plotting script where I will plot theta_2 as a function of
% potential

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
gamma = 8.16*10^(-6);                                                       % mol/m^2 [concentration of active sites]
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
%% ################## Theta vs Potential #############################

% Transforming time to potential fot the interpolation
potential_interpol = CV_potential_alkaline(Schalenbach_time_CV_linear, "array");

% Creating a string element for the legends
string_array_1 = sprintf('$k^{0}_{4+}$ = %.1f', round(k_4_0_plus(1), 5));
string_array_2 = sprintf('$k^{0}_{4+}$ = %.2f', round(k_4_0_plus(2), 5));
string_array_3 = sprintf('$k^{0}_{4+}$ = %.3f', round(k_4_0_plus(3), 5));

%% Cherevko
[t_cherevko_1, theta_cherevko_1, potential_cherevko_1, theta_cherevko_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(1));
[t_cherevko_2, theta_cherevko_2, potential_cherevko_2, theta_cherevko_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(2));
[t_cherevko_3, theta_cherevko_3, potential_cherevko_3, theta_cherevko_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Cherevko_E_alkaline, Cherevko_i_alkaline, Cherevko_OH_alkaline, Cherevko_T_alkaline, "Linear", k_4_0_plus(3));

figure('Name', 'Cherevko: theta_2 vs potential')                                % Creating figure
%yyaxis left
plot(potential_cherevko_1, theta_cherevko_1, "Color", "red")                    % Plots the line for 1
hold on
scatter(potential_interpol, theta_cherevko_interpol_1,...                       % Scatter interpolated values for 1
    45,"red", 'o')                                                                      
plot(potential_cherevko_2, theta_cherevko_2, "Color", "blue")                   % Plots the line for 2
scatter(potential_interpol, theta_cherevko_interpol_2,...                       % Scatter interpolated values for 2
   45,"blue", 'square')     
plot(potential_cherevko_3, theta_cherevko_3, "Color", "green")                  % Plots the line for 3
scatter(potential_interpol, theta_cherevko_interpol_3,...                       % Scatter the interpolated values for 3
   45, "green", 'diamond')  
hold off
ax_cherevko_alkaline = gca; % current axes                                      % Creating an ax with gca such that the fontsize can be changed
ax_cherevko_alkaline.XAxis.FontSize = 12;                                       % Changing the tick size on the x-axis
ax_cherevko_alkaline.YAxis.FontSize = 12;                                       % Changing the tick size on the y-axis
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
xlim([min(potential_cherevko_3) max(potential_cherevko_3)])
ylim([min(theta_cherevko_3)*0 max(theta_cherevko_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol(22);                                                  % x_begin for arrow
x2_1 = (potential_interpol(23) + x1_1)/2;                                       % x_end for arrow
y1_1 = theta_cherevko_interpol_1(22);                                           % y_begin for arrow
y2_1 = (theta_cherevko_interpol_1(23)+y1_1)/2;                                  % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                                 % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                                 % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                                 % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                                 % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                      % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                       % Defines position
arh1.Color = 'red';                                                             % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol(end-2);                                               % x_begin for arrow
x2_2 = (potential_interpol(end-1) + x1_2)/2;                                    % x_end for arrow
y1_2 = theta_cherevko_interpol_2(end-2);                                        % y_begin for arrow
y2_2 = (theta_cherevko_interpol_2(end-1)+y1_2)/2;                               % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                                 % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                                 % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                                 % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                                 % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                      % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                       % Defines position
arh2.Color = 'blue';                                                            % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol(end-2);                                               % x_begin for arrow
x2_3 = (potential_interpol(end-1) + x1_3)/2;                                    % x_end for arrow
y1_3 = theta_cherevko_interpol_3(end-2);                                        % y_begin for arrow
y2_3 = (theta_cherevko_interpol_3(end-1)+y1_3)/2;                               % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                                 % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                                 % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                                 % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                                 % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                      % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                       % Defines position
arh3.Color = 'green';                                                           % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyaxis right
%ax_2_cherevko_alkaline = gca;
%ax_2_cherevko_alkaline.YAxis.FontSize = 12;
plot(potential_interpol, Schalenbach_dissolution_mole,...                         % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex')                                % Label for second y_axis
annotation('textbox', [.15 .676 .1 .1], 'String',["Cherevko -", "Alkaline"],...  % Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);

%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol(8);                                                   % x_begin for arrow
x2_d = (potential_interpol(9) + x1_d)/2;                                        % x_end for arrow
y1_d = Schalenbach_dissolution_mole(8);                                         % y_begin for arrow
y2_d = (Schalenbach_dissolution_mole(9)+y1_d)/2;                                % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                                 % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                                 % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                                 % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                                 % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                      % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                       % Defines position
arhd.Color = [.5 .5 .5];                                                        % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend({'',string_array_1, '', string_array_2, '',...
    string_array_3, "$\frac{d Ir}{d t}$"},...                                   % Creating a legend for the graphs
    'Position', [.2375 .45 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------

%% Damjanovic
[t_damj_1, theta_damj_1, potential_damj_1, theta_damj_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(1));
[t_damj_2, theta_damj_2, potential_damj_2, theta_damj_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(2));
[t_damj_3, theta_damj_3, potential_damj_3, theta_damj_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Linear", k_4_0_plus(3));


figure('Name', 'Damjanovic: theta_2 vs potential')                              % Creating figure
%yyaxis left
plot(potential_damj_1, theta_damj_1, "Color", "red")                            % Plots the line for 1
hold on
scatter(potential_interpol, theta_damj_interpol_1,...                           % Scatter interpolated values for 1
    45,"red", 'o')                                                                      
plot(potential_damj_2, theta_damj_2, "Color", "blue")                           % Plots the line for 2
scatter(potential_interpol, theta_damj_interpol_2,...                           % Scatter interpolated values for 2
    45,"blue", 'square')     
plot(potential_damj_3, theta_damj_3, "Color", "green")                          % Plots the line for 3
scatter(potential_interpol, theta_damj_interpol_3,...                           % Scatter the interpolated values for 3
    45, "green", 'diamond')  
%hold off
ax_damj_alkaline = gca; % current axes                                          % Creating an ax with gca such that the fontsize can be changed
ax_damj_alkaline.XAxis.FontSize = 12;                                           % Changing the tick size on the x-axis
ax_damj_alkaline.YAxis.FontSize = 12;                                           % Changing the tick size on the y-axis
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
xlim([min(potential_damj_3) max(potential_damj_3)])
ylim([min(theta_damj_3)*0 max(theta_damj_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol(22);                                                  % x_begin for arrow
x2_1 = (potential_interpol(23) + x1_1)/2;                                       % x_end for arrow
y1_1 = theta_damj_interpol_1(22);                                               % y_begin for arrow
y2_1 = (theta_damj_interpol_1(23)+y1_1)/2;                                      % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                                 % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                                 % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                                 % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                                 % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                      % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                       % Defines position
arh1.Color = 'red';                                                             % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol(end-2);                                               % x_begin for arrow
x2_2 = (potential_interpol(end-1) + x1_2)/2;                                    % x_end for arrow
y1_2 = theta_damj_interpol_2(end-2);                                            % y_begin for arrow
y2_2 = (theta_damj_interpol_2(end-1)+y1_2)/2;                                   % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                                 % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                                 % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                                 % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                                 % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                      % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                       % Defines position
arh2.Color = 'blue';                                                            % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol(end-2);                                               % x_begin for arrow
x2_3 = (potential_interpol(end-1) + x1_3)/2;                                    % x_end for arrow
y1_3 = theta_damj_interpol_3(end-2);                                            % y_begin for arrow
y2_3 = (theta_damj_interpol_3(end-1)+y1_3)/2;                                   % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                                 % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                                 % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                                 % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                                 % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                      % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                       % Defines position
arh3.Color = 'green';                                                           % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyaxis right
%ax_2_damj_alkaline = gca;
%ax_2_damj_alkaline.YAxis.FontSize = 12;
plot(potential_interpol, Schalenbach_dissolution_mole,...                        % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex')      % Label for second y_axis
annotation('textbox', [.15 .676 .1 .1], 'String',["Damjanovic -", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);

%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol(8);                                                   % x_begin for arrow
x2_d = (potential_interpol(9) + x1_d)/2;                                        % x_end for arrow
y1_d = Schalenbach_dissolution_mole(8);                                         % y_begin for arrow
y2_d = (Schalenbach_dissolution_mole(9)+y1_d)/2;                                % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                                 % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                                 % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                                 % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                                 % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                      % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                       % Defines position
arhd.Color = [.5 .5 .5];                                                        % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend({'',string_array_1, '', string_array_2, '',...
    string_array_3, "$\frac{d Ir}{d t}$"},...                                   % Creating a legend for the graphs
    'Position', [.2375 .45 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------

%% Damjanovic log
[t_damj_log_1, theta_damj_log_1, potential_damj_log_1, theta_damj_log_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(1));
[t_damj_log_2, theta_damj_log_2, potential_damj_log_2, theta_damj_log_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(2));
[t_damj_log_3, theta_damj_log_3, potential_damj_log_3, theta_damj_log_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Damjanovic_E_alkaline, Damjanovic_i_alkaline, Damjanovic_OH_alkaline, Damjanovic_T_alkaline, "Logarithmic", k_4_0_plus(3));


figure('Name', 'Damjanovic log: theta_2 vs potential')                          % Creating figure
%yyaxis left
plot(potential_damj_log_1, theta_damj_log_1, "Color", "red")                    % Plots the line for 1
hold on
scatter(potential_interpol, theta_damj_log_interpol_1,...                       % Scatter interpolated values for 1
    45,"red", 'o')                                                                      
plot(potential_damj_log_2, theta_damj_log_2, "Color", "blue")                   % Plots the line for 2
scatter(potential_interpol, theta_damj_log_interpol_2,...                       % Scatter interpolated values for 2
    45,"blue", 'square')     
plot(potential_damj_log_3, theta_damj_log_3, "Color", "green")                  % Plots the line for 3
scatter(potential_interpol, theta_damj_log_interpol_3,...                       % Scatter the interpolated values for 3
    45, "green", 'diamond')  
%hold off
ax_damj_log_alkaline = gca; % current axes                                      % Creating an ax with gca such that the fontsize can be changed
ax_damj_log_alkaline.XAxis.FontSize = 12;                                       % Changing the tick size on the x-axis
ax_damj_log_alkaline.YAxis.FontSize = 12;                                       % Changing the tick size on the y-axis
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
xlim([min(potential_damj_log_3) max(potential_damj_log_3)])
ylim([min(theta_damj_log_3)*0 max(theta_damj_log_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol(22);                                                  % x_begin for arrow
x2_1 = (potential_interpol(23) + x1_1)/2;                                       % x_end for arrow
y1_1 = theta_damj_log_interpol_1(22);                                           % y_begin for arrow
y2_1 = (theta_damj_log_interpol_1(23)+y1_1)/2;                                  % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                                 % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                                 % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                                 % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                                 % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                      % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                       % Defines position
arh1.Color = 'red';                                                             % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol(end-2);                                               % x_begin for arrow
x2_2 = (potential_interpol(end-1) + x1_2)/2;                                    % x_end for arrow
y1_2 = theta_damj_log_interpol_2(end-2);                                        % y_begin for arrow
y2_2 = (theta_damj_log_interpol_2(end-1)+y1_2)/2;                               % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                                 % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                                 % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                                 % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                                 % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                      % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                       % Defines position
arh2.Color = 'blue';                                                            % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol(end-2);                                               % x_begin for arrow
x2_3 = (potential_interpol(end-1) + x1_3)/2;                                    % x_end for arrow
y1_3 = theta_damj_log_interpol_3(end-2);                                        % y_begin for arrow
y2_3 = (theta_damj_log_interpol_3(end-1)+y1_3)/2;                               % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                                 % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                                 % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                                 % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                                 % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                      % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                       % Defines position
arh3.Color = 'green';                                                           % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyaxis right
%ax_2_damj_alkaline = gca;
%ax_2_damj_alkaline.YAxis.FontSize = 12;
plot(potential_interpol, Schalenbach_dissolution_mole,...                           % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex')         % Label for second y_axis
annotation('textbox', [.15 .676 .1 .1], 'String',["Damjanovic log-", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);

%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol(8);                                                   % x_begin for arrow
x2_d = (potential_interpol(9) + x1_d)/2;                                        % x_end for arrow
y1_d = Schalenbach_dissolution_mole(8);                                         % y_begin for arrow
y2_d = (Schalenbach_dissolution_mole(9)+y1_d)/2;                                % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                                 % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                                 % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                                 % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                                 % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                      % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                       % Defines position
arhd.Color = [.5 .5 .5];                                                        % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend({'',string_array_1, '', string_array_2, '',...
    string_array_3, "$\frac{d Ir}{d t}$"},...                                   % Creating a legend for the graphs
    'Position', [.2375 .45 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------
%
%% Schalenbach
[t_schalenbach_1, theta_schalenbach_1, potential_schalenbach_1, theta_schalenbach_interpol_1] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(1));
[t_schalenbach_2, theta_schalenbach_2, potential_schalenbach_2, theta_schalenbach_interpol_2] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(2));
[t_schalenbach_3, theta_schalenbach_3, potential_schalenbach_3, theta_schalenbach_interpol_3] =...
    time_theta_potential_ode15s_alkaline(Schalenbach_E_alkaline, Schalenbach_i_alkaline, Schalenbach_OH_alkaline, Schalenbach_T_alkaline, "Linear", k_4_0_plus(3));

figure('Name', 'Schalenbach: theta_2 vs potential')                              % Creating figure
%yyaxis left
plot(potential_schalenbach_1, theta_schalenbach_1, "Color", "red")               % Plots the line for 1
hold on
scatter(potential_interpol, theta_schalenbach_interpol_1,...                     % Scatter interpolated values for 1
    45,"red", 'o')                                                                      
plot(potential_schalenbach_2, theta_schalenbach_2, "Color", "blue")              % Plots the line for 2
scatter(potential_interpol, theta_schalenbach_interpol_2,...                     % Scatter interpolated values for 2
    45,"blue", 'square')     
plot(potential_schalenbach_3, theta_schalenbach_3, "Color", "green")             % Plots the line for 3
scatter(potential_interpol, theta_schalenbach_interpol_3,...                     % Scatter the interpolated values for 3
    45, "green", 'diamond')  
%hold off
ax_schalenbach_alkaline = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_schalenbach_alkaline.XAxis.FontSize = 12;                                     % Changing the tick size on the x-axis
ax_schalenbach_alkaline.YAxis.FontSize = 12;                                     % Changing the tick size on the y-axis
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex')
ylabel('$\theta_{2}(t)$ - [$-$]','Interpreter','latex')
xlim([min(potential_schalenbach_1) max(potential_schalenbach_1)])
ylim([min(theta_schalenbach_3)*0 max(theta_schalenbach_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol(22);                                                  % x_begin for arrow
x2_1 = (potential_interpol(23) + x1_1)/2;                                       % x_end for arrow
y1_1 = theta_schalenbach_interpol_1(22);                                        % y_begin for arrow
y2_1 = (theta_schalenbach_interpol_1(23)+y1_1)/2;                               % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                                 % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                                 % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                                 % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                                 % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                      % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                       % Defines position
arh1.Color = 'red';                                                             % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol(end-2);                                               % x_begin for arrow
x2_2 = (potential_interpol(end-1) + x1_2)/2;                                    % x_end for arrow
y1_2 = theta_schalenbach_interpol_2(end-2);                                     % y_begin for arrow
y2_2 = (theta_schalenbach_interpol_2(end-1)+y1_2)/2;                            % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                                 % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                                 % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                                 % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                                 % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                      % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                       % Defines position
arh2.Color = 'blue';                                                            % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol(end-2);                                               % x_begin for arrow
x2_3 = (potential_interpol(end-1) + x1_3)/2;                                    % x_end for arrow
y1_3 = theta_schalenbach_interpol_3(end-2);                                     % y_begin for arrow
y2_3 = (theta_schalenbach_interpol_3(end-1)+y1_3)/2;                            % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                                 % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                                 % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                                 % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                                 % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                      % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                       % Defines position
arh3.Color = 'green';                                                           % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis right
%ax_2_damj_alkaline = gca;
%ax_2_damj_alkaline.YAxis.FontSize = 12;
plot(potential_interpol, Schalenbach_dissolution_mole,...                        % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex')      % Label for second y_axis
annotation('textbox', [.15 .676 .1 .1], 'String',["Schalenbach-", "Alkaline"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol(8);                                                   % x_begin for arrow
x2_d = (potential_interpol(9) + x1_d)/2;                                        % x_end for arrow
y1_d = Schalenbach_dissolution_mole(8);                                         % y_begin for arrow
y2_d = (Schalenbach_dissolution_mole(9)+y1_d)/2;                                % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                                 % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                                 % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                                 % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                                 % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',12, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                      % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                       % Defines position
arhd.Color = [.5 .5 .5];                                                        % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend({'',string_array_1, '', string_array_2, '',...
    string_array_3, "$\frac{d Ir}{d t}$"},...                                   % Creating a legend for the graphs
    'Position', [.2375 .45 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------