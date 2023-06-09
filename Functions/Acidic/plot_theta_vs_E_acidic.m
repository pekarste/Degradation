% This will be a plotting script where I will plot theta_2 as a function of
% potential

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
Damjanovic_acidic = readmatrix("Data\Acidic\Damjanovic_Ir_E_vs_log_i_acidic.xlsx"); % Current density/potential data from Damjanovic
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
string_array_1 = sprintf('$k^{0}_{4+}$ = %.2f $s^{-1}$', round(k_4_0_plus(1), 5));
string_array_2 = sprintf('$k^{0}_{4+}$ = %.3f $s^{-1}$', round(k_4_0_plus(2), 5));
string_array_3 = sprintf('$k^{0}_{4+}$ = %.4f $s^{-1}$', round(k_4_0_plus(3), 5));
string_array_4 = '$\frac{d Ir}{d t}$';

% Colour blind pallette
Orange          = [.90 .60 .0];                                        % Orange                                        
Reddish_purple  = [.80 .60 .70];                                       % Reddish purple
Sky_blue        = [.35 .70 .90];                                       % Sky blue

%% ####################       Scohy        ################################
[t_scohy_1, gamma_theta_scohy_1, potential_scohy_1, theta_scohy_interpol_1] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(1));
[t_scohy_2, gamma_theta_scohy_2, potential_scohy_2, theta_scohy_interpol_2] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(2));
[t_scohy_3, gamma_theta_scohy_3, potential_scohy_3, theta_scohy_interpol_3] =...
    time_theta_potential_ode15s_acidic(Scohy_E_acidic, Scohy_i_acidic, Scohy_a_H_plus, Scohy_T_acidic, "Linear", k_4_0_plus(3));

figure('Name', 'Scohy: theta_2 vs potential')                               % Creating figure
%yyaxis left
plot(potential_scohy_1, gamma_theta_scohy_1, "Color", Orange)                % Plots the line for 1
hold on
scatter(potential_interpol_acidic, theta_scohy_interpol_1,...               % Scatter interpolated values for 1
   45,Orange, 'o', 'filled')                                                                      
plot(potential_scohy_2, gamma_theta_scohy_2, "Color", Reddish_purple)               % Plots the line for 2
scatter(potential_interpol_acidic, theta_scohy_interpol_2,...               % Scatter interpolated values for 2
   45,Reddish_purple, 'square', 'filled')     
plot(potential_scohy_3, gamma_theta_scohy_3, "Color", Sky_blue)              % Plots the line for 3
scatter(potential_interpol_acidic, theta_scohy_interpol_3,...               % Scatter the interpolated values for 3
   45, Sky_blue, 'diamond', 'filled')  
hold off
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex',...
    'FontSize',15)
ylabel('$\Gamma\theta_{2}(t)$ - [$\frac{mol}{m^{2}}$]',...
    'FontSize', 15,'Interpreter','latex')
ax_scohy_acidic = gca; % current axes                                       % Creating an ax with gca such that the fontsize can be changed
ax_scohy_acidic.XAxis.FontSize = 15;                                        % Changing the tick size on the x-axis
ax_scohy_acidic.YAxis(1).FontSize = 15;                                        % Changing the tick size on the y-axis
xlim([min(potential_interpol_acidic) 1.6])
%ylim([min(gamma_theta_scohy_3)*0 max(gamma_theta_scohy_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                  % x_lim for normalising position
yL = ylim;                                                                  % y_lim for normalasing position
ah = gca;                                                                   % gives current axis handle (left axis)
aPos = ah.Position;                                                         % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                           % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                           % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol_acidic(12);                                       % x_begin for arrow
x2_1 = (potential_interpol_acidic(13) + x1_1)/2;                            % x_end for arrow
y1_1 = theta_scohy_interpol_1(12);                                          % y_begin for arrow
y2_1 = (theta_scohy_interpol_1(13)+y1_1)/2;                                 % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                             % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                             % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                             % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                             % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                  % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                   % Defines position
arh1.Color = Orange;                                                         % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol_acidic(end-2);                                    % x_begin for arrow
x2_2 = (potential_interpol_acidic(end-1) + x1_2)/2;                         % x_end for arrow
y1_2 = theta_scohy_interpol_2(end-2);                                       % y_begin for arrow
y2_2 = (theta_scohy_interpol_2(end-1)+y1_2)/2;                              % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                             % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                             % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                             % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                             % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                  % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                   % Defines position
arh2.Color = Reddish_purple;                                                        % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol_acidic(end-2);                                    % x_begin for arrow
x2_3 = (potential_interpol_acidic(end-1) + x1_3)/2;                         % x_end for arrow
y1_3 = theta_scohy_interpol_3(end-2);                                       % y_begin for arrow
y2_3 = (theta_scohy_interpol_3(end-1)+y1_3)/2;                              % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                             % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                             % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                             % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                             % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                  % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                   % Defines position
arh3.Color = Sky_blue;                                                       % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyaxis right
plot(potential_interpol_acidic, Mayrhofer_dissolution_mole,...            % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ax_scohy_acidic.YAxis(2).FontSize = 15;                                        % Changing the tick size on the y-axis
ax_scohy_acidic.YAxis(2).Color = 'black';
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]',...                    % Label for second y_axis
    'Interpreter','latex', 'FontSize',20) 

%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                  % x_lim for normalising position
yL = ylim;                                                                  % y_lim for normalasing position
ah = gca;                                                                   % gives current axis handle (left axis)
aPos = ah.Position;                                                         % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                           % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                           % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol_acidic(8);                                        % x_begin for arrow
x2_d = (potential_interpol_acidic(9) + x1_d)/2;                             % x_end for arrow
y1_d = Mayrhofer_dissolution_mole(8);                                     % y_begin for arrow
y2_d = (Mayrhofer_dissolution_mole(9)+y1_d)/2;                            % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                             % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                             % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                             % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                             % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                  % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                   % Defines position
arhd.Color = [.5 .5 .5];                                                    % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% legend({'',string_array_1, '', string_array_2, '',...
%     string_array_3, "$\frac{d Ir}{d t}$"},...                                % Creating a legend for the graphs
%     'Position', [.25 .675 .1 .1],'Interpreter','latex', 'FontSize',15)

annotation('textbox', [.15 .67 .1 .1], 'String',["Scohy -", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.60 .14 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange, 'Rotation',45);
annotation('textbox', [.25 .45 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple, 'Rotation',25);
annotation('textbox', [.25 .75 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue, 'Rotation',5);
annotation('textbox', [.70 .49 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',20, 'Color',[.5 .5 .5]);
%--------------------------------------------------------------------------

%% ##################         Damjanovic          #########################
[t_damj_1, gamma_theta_damj_1, potential_damj_1, theta_damj_interpol_1] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(1));
[t_damj_2, gamma_theta_damj_2, potential_damj_2, theta_damj_interpol_2] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(2));
[t_damj_3, gamma_theta_damj_3, potential_damj_3, theta_damj_interpol_3] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Linear", k_4_0_plus(3));


figure('Name', 'Damjanovic: theta_2 vs potential')                          % Creating figure
%yyaxis left
plot(potential_damj_1, gamma_theta_damj_1, "Color", Orange)                  % Plots the line for 1
hold on
scatter(potential_interpol_acidic, theta_damj_interpol_1,...                % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(potential_damj_2, gamma_theta_damj_2, "Color", Reddish_purple)                 % Plots the line for 2
scatter(potential_interpol_acidic, theta_damj_interpol_2,...                % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(potential_damj_3, gamma_theta_damj_3, "Color", Sky_blue)                % Plots the line for 3
scatter(potential_interpol_acidic, theta_damj_interpol_3,...                % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_damj_acidic = gca; % current axes                                        % Creating an ax with gca such that the fontsize can be changed
ax_damj_acidic.XAxis.FontSize = 15;                                         % Changing the tick size on the x-axis
ax_damj_acidic.YAxis(1).FontSize = 15;                                         % Changing the tick size on the y-axis
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex')
ylabel('$\Gamma\theta_{2}(t)$ - [$\frac{mol}{m^{2}}$]','Interpreter','latex')
xlim([min(potential_interpol_acidic) 1.6])
%ylim([min(gamma_theta_damj_3)*0 max(gamma_theta_damj_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                  % x_lim for normalising position
yL = ylim;                                                                  % y_lim for normalasing position
ah = gca;                                                                   % gives current axis handle (left axis)
aPos = ah.Position;                                                         % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                           % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                           % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol_acidic(11);                                       % x_begin for arrow
x2_1 = (potential_interpol_acidic(12) + x1_1)/2;                            % x_end for arrow
y1_1 = theta_damj_interpol_1(11);                                           % y_begin for arrow
y2_1 = (theta_damj_interpol_1(12)+y1_1)/2;                                  % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                             % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                             % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                             % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                             % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                  % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                   % Defines position
arh1.Color = Orange;                                                         % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol_acidic(end-2);                                    % x_begin for arrow
x2_2 = (potential_interpol_acidic(end-1) + x1_2)/2;                         % x_end for arrow
y1_2 = theta_damj_interpol_2(end-2);                                        % y_begin for arrow
y2_2 = (theta_damj_interpol_2(end-1)+y1_2)/2;                               % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                             % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                             % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                             % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                             % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                  % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                   % Defines position
arh2.Color = Reddish_purple;                                                        % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol_acidic(end-2);                                    % x_begin for arrow
x2_3 = (potential_interpol_acidic(end-1) + x1_3)/2;                         % x_end for arrow
y1_3 = theta_damj_interpol_3(end-2);                                        % y_begin for arrow
y2_3 = (theta_damj_interpol_3(end-1)+y1_3)/2;                               % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                             % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                             % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                             % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                             % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                  % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                   % Defines position
arh3.Color = Sky_blue;                                                       % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyaxis right
%ax_2_damj_acidic = gca;
ax_damj_acidic.YAxis(2).FontSize = 15;
ax_damj_acidic.YAxis(2).Color = 'black';
plot(potential_interpol_acidic, Mayrhofer_dissolution_mole,...            % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex') % Label for second y_axis

%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                  % x_lim for normalising position
yL = ylim;                                                                  % y_lim for normalasing position
ah = gca;                                                                   % gives current axis handle (left axis)
aPos = ah.Position;                                                         % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                           % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                           % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol_acidic(8);                                        % x_begin for arrow
x2_d = (potential_interpol_acidic(9) + x1_d)/2;                             % x_end for arrow
y1_d = Mayrhofer_dissolution_mole(8);                                     % y_begin for arrow
y2_d = (Mayrhofer_dissolution_mole(9)+y1_d)/2;                            % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                             % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                             % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                             % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                             % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                  % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                   % Defines position
arhd.Color = [.5 .5 .5];                                                    % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% legend({'',string_array_1, '', string_array_2, '',...
%     string_array_3, "$\frac{d Ir}{d t}$"},...                               % Creating a legend for the graphs
%     'Position', [.2375 .675 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------
annotation('textbox', [.20 .6 .1 .1], 'String',["Damjanovic -", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.70 .15 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange, 'Rotation',55);
annotation('textbox', [.27 .30 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple, 'Rotation',25);
annotation('textbox', [.24 .7 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue, 'Rotation',5);
annotation('textbox', [.725 .46 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',20, 'Color',[.5 .5 .5]);
 
%% ################        Damjanovic log        #########################
[t_damj_log_1, gamma_theta_damj_log_1, potential_damj_log_1, theta_damj_log_interpol_1] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(1));
[t_damj_log_2, gamma_theta_damj_log_2, potential_damj_log_2, theta_damj_log_interpol_2] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(2));
[t_damj_log_3, gamma_theta_damj_log_3, potential_damj_log_3, theta_damj_log_interpol_3] =...
    time_theta_potential_ode15s_acidic(Damjanovic_E_acidic, Damjanovic_i_acidic, Damjanovic_a_H_plus, Damjanovic_T_acidic, "Logarithmic", k_4_0_plus(3));


figure('Name', 'Damjanovic log: theta_2 vs potential')                      % Creating figure
%yyaxis left
plot(potential_damj_log_1, gamma_theta_damj_log_1, "Color", Orange)          % Plots the line for 1
hold on
scatter(potential_interpol_acidic, theta_damj_log_interpol_1,...            % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(potential_damj_log_2, gamma_theta_damj_log_2, "Color", Reddish_purple)         % Plots the line for 2
scatter(potential_interpol_acidic, theta_damj_log_interpol_2,...            % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(potential_damj_log_3, gamma_theta_damj_log_3, "Color", Sky_blue)        % Plots the line for 3
scatter(potential_interpol_acidic, theta_damj_log_interpol_3,...            % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_damj_log_acidic = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_damj_log_acidic.XAxis.FontSize = 15;                                     % Changing the tick size on the x-axis
ax_damj_log_acidic.YAxis.FontSize = 15;                                     % Changing the tick size on the y-axis
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex')
ylabel('$\Gamma\theta_{2}(t)$ - [$\frac{mol}{m^{2}}$]','Interpreter','latex')
xlim([min(potential_interpol_acidic) 1.6])
%ylim([min(gamma_theta_damj_log_3)*0 max(gamma_theta_damj_log_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                  % x_lim for normalising position
yL = ylim;                                                                  % y_lim for normalasing position
ah = gca;                                                                   % gives current axis handle (left axis)
aPos = ah.Position;                                                         % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                           % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                           % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol_acidic(11);                                       % x_begin for arrow
x2_1 = (potential_interpol_acidic(12) + x1_1)/2;                            % x_end for arrow
y1_1 = theta_damj_log_interpol_1(11);                                       % y_begin for arrow
y2_1 = (theta_damj_log_interpol_1(12)+y1_1)/2;                              % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                             % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                             % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                             % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                             % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                  % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                   % Defines position
arh1.Color = Orange;                                                         % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol_acidic(end-2);                                    % x_begin for arrow
x2_2 = (potential_interpol_acidic(end-1) + x1_2)/2;                         % x_end for arrow
y1_2 = theta_damj_log_interpol_2(end-2);                                    % y_begin for arrow
y2_2 = (theta_damj_log_interpol_2(end-1)+y1_2)/2;                           % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                             % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                             % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                             % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                             % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                  % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                   % Defines position
arh2.Color = Reddish_purple;                                                        % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol_acidic(end-2);                                    % x_begin for arrow
x2_3 = (potential_interpol_acidic(end-1) + x1_3)/2;                         % x_end for arrow
y1_3 = theta_damj_log_interpol_3(end-2);                                    % y_begin for arrow
y2_3 = (theta_damj_log_interpol_3(end-1)+y1_3)/2;                           % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                             % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                             % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                             % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                             % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                  % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                   % Defines position
arh3.Color = Sky_blue;                                                       % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyaxis right
ax_damj_log_acidic.YAxis(2).FontSize = 15;
ax_damj_log_acidic.YAxis(2).Color = 'black';
plot(potential_interpol_acidic, Mayrhofer_dissolution_mole,...            % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex') % Label for second y_axis


%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                  % x_lim for normalising position
yL = ylim;                                                                  % y_lim for normalasing position
ah = gca;                                                                   % gives current axis handle (left axis)
aPos = ah.Position;                                                         % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                           % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                           % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol_acidic(8);                                        % x_begin for arrow
x2_d = (potential_interpol_acidic(9) + x1_d)/2;                             % x_end for arrow
y1_d = Mayrhofer_dissolution_mole(8);                                     % y_begin for arrow
y2_d = (Mayrhofer_dissolution_mole(9)+y1_d)/2;                            % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                             % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                             % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                             % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                             % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                           % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                  % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                   % Defines position
arhd.Color = [.5 .5 .5];                                                    % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%legend({'',string_array_1, '', string_array_2, '',...
%    string_array_3, "$\frac{d Ir}{d t}$"},...                               % Creating a legend for the graphs
%    'Position', [.2375 .45 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------
annotation('textbox', [.15 .60 .1 .1], 'String',["Damjanovic log-", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.70 .15 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange, 'Rotation',55);
annotation('textbox', [.27 .30 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple, 'Rotation',25);
annotation('textbox', [.25 .69 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue, 'Rotation',5);
annotation('textbox', [.725 .46 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',20, 'Color',[.5 .5 .5]);

%% ####################       Cherevko        #############################
[t_cherevko_1, gamma_theta_cherevko_1, potential_cherevko_1, theta_cherevko_interpol_1] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(1));
[t_cherevko_2, gamma_theta_cherevko_2, potential_cherevko_2, theta_cherevko_interpol_2] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(2));
[t_cherevko_3, gamma_theta_cherevko_3, potential_cherevko_3, theta_cherevko_interpol_3] =...
    time_theta_potential_ode15s_acidic(Cherevko_E_acidic, Cherevko_i_acidic, Cherevko_a_H_plus, Cherevko_T_acidic, "Linear", k_4_0_plus(3));

figure('Name', 'Cherevko: theta_2 vs potential')                              % Creating figure
%yyaxis left
plot(potential_cherevko_1, gamma_theta_cherevko_1, "Color", Orange)               % Plots the line for 1
hold on
scatter(potential_interpol_acidic, theta_cherevko_interpol_1,...                     % Scatter interpolated values for 1
    45,Orange, 'o', 'filled')                                                                      
plot(potential_cherevko_2, gamma_theta_cherevko_2, "Color", Reddish_purple)              % Plots the line for 2
scatter(potential_interpol_acidic, theta_cherevko_interpol_2,...                     % Scatter interpolated values for 2
    45,Reddish_purple, 'square', 'filled')     
plot(potential_cherevko_3, gamma_theta_cherevko_3, "Color", Sky_blue)             % Plots the line for 3
scatter(potential_interpol_acidic, theta_cherevko_interpol_3,...                     % Scatter the interpolated values for 3
    45, Sky_blue, 'diamond', 'filled')  
%hold off
ax_cherevko_acidic = gca; % current axes                                    % Creating an ax with gca such that the fontsize can be changed
ax_cherevko_acidic.XAxis.FontSize = 15;                                     % Changing the tick size on the x-axis
ax_cherevko_acidic.YAxis(1).FontSize = 15;                                     % Changing the tick size on the y-axis
xlabel('Potential -E vs RHE [$V$]','Interpreter','latex')
ylabel('$\Gamma\theta_{2}(t)$ - [$\frac{mol}{m^{2}}$]','Interpreter','latex')
xlim([min(potential_interpol_acidic) 1.6])
%ylim([min(gamma_theta_cherevko_3)*0 max(gamma_theta_cherevko_3)])

%%%%%%%%%%%%%%%%%%%  Creating arrowheads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead 1 -------------------------------------------------------------
x1_1 = potential_interpol_acidic(11);                                                  % x_begin for arrow
x2_1 = (potential_interpol_acidic(12) + x1_1)/2;                                       % x_end for arrow
y1_1 = theta_cherevko_interpol_1(11);                                        % y_begin for arrow
y2_1 = (theta_cherevko_interpol_1(12)+y1_1)/2;                               % y_end for arrow

x1p_1 = interp1(xL, ahx, x1_1);                                                 % These lines gives the interpolated values
x2p_1 = interp1(xL, ahx, x2_1);                                                 % I belive this has something to do with the
y1p_1 = interp1(yL, ahy, y1_1);                                                 % definition of the angle of the arrowhead and
y2p_1 = interp1(yL, ahy, y2_1);                                                 % the use of normalised coordinates

arh1 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh1.Units = 'normalized';                                                      % Normalaizing the units
arh1.Position = [x1p_1, y1p_1, x2p_1-x1p_1, y2p_1-y1p_1];                       % Defines position
arh1.Color = Orange;                                                             % Defines colour for arrowhead

% Arrowhead 2 -------------------------------------------------------------
x1_2 = potential_interpol_acidic(end-2);                                               % x_begin for arrow
x2_2 = (potential_interpol_acidic(end-1) + x1_2)/2;                                    % x_end for arrow
y1_2 = theta_cherevko_interpol_2(end-2);                                     % y_begin for arrow
y2_2 = (theta_cherevko_interpol_2(end-1)+y1_2)/2;                            % y_end for arrow

x1p_2 = interp1(xL, ahx, x1_2);                                                 % These lines gives the interpolated values
x2p_2 = interp1(xL, ahx, x2_2);                                                 % I belive this has something to do with the
y1p_2 = interp1(yL, ahy, y1_2);                                                 % definition of the angle of the arrowhead and
y2p_2 = interp1(yL, ahy, y2_2);                                                 % the use of normalised coordinates

arh2 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh2.Units = 'normalized';                                                      % Normalaizing the units
arh2.Position = [x1p_2, y1p_2, x2p_2-x1p_2, y2p_2-y1p_2];                       % Defines position
arh2.Color = Reddish_purple;                                                            % Defines colour for arrowhead

% Arrowhead 3 -------------------------------------------------------------
x1_3 = potential_interpol_acidic(end-2);                                               % x_begin for arrow
x2_3 = (potential_interpol_acidic(end-1) + x1_3)/2;                                    % x_end for arrow
y1_3 = theta_cherevko_interpol_3(end-2);                                     % y_begin for arrow
y2_3 = (theta_cherevko_interpol_3(end-1)+y1_3)/2;                            % y_end for arrow

x1p_3 = interp1(xL, ahx, x1_3);                                                 % These lines gives the interpolated values
x2p_3 = interp1(xL, ahx, x2_3);                                                 % I belive this has something to do with the
y1p_3 = interp1(yL, ahy, y1_3);                                                 % definition of the angle of the arrowhead and
y2p_3 = interp1(yL, ahy, y2_3);                                                 % the use of normalised coordinates

arh3 = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arh3.Units = 'normalized';                                                      % Normalaizing the units
arh3.Position = [x1p_3, y1p_3, x2p_3-x1p_3, y2p_3-y1p_3];                       % Defines position
arh3.Color = Sky_blue;                                                           % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis right
ax_cherevko_acidic.YAxis(2).FontSize = 15;
ax_cherevko_acidic.YAxis(2).Color = 'black';
plot(potential_interpol_acidic, Mayrhofer_dissolution_mole,...                        % Plots the potential regime
    'color', [.5 .5 .5], 'LineWidth', 1.5, 'LineStyle', '--')
ylabel('$\frac{d Ir}{d t}$ - [$\frac{mol}{m^{2}s}$]','Interpreter','latex')      % Label for second y_axis

%%%%%%%%%%%%%%%%%%%%%%%%% creating arrowhead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = xlim;                                                                      % x_lim for normalising position
yL = ylim;                                                                      % y_lim for normalasing position
ah = gca;                                                                       % gives current axis handle (left axis)
aPos = ah.Position;                                                             % These three lines gives the 
ahx = [aPos(1), aPos(1)+aPos(3)];                                               % position in some way
ahy = [aPos(2), aPos(2)+aPos(4)];                                               % I don-t really know how they work

% Arrowhead diss-----------------------------------------------------------
x1_d = potential_interpol_acidic(8);                                                   % x_begin for arrow
x2_d = (potential_interpol_acidic(9) + x1_d)/2;                                        % x_end for arrow
y1_d = Mayrhofer_dissolution_mole(8);                                         % y_begin for arrow
y2_d = (Mayrhofer_dissolution_mole(9)+y1_d)/2;                                % y_end for arrow

x1p_d = interp1(xL, ahx, x1_d);                                                 % These lines gives the interpolated values
x2p_d = interp1(xL, ahx, x2_d);                                                 % I belive this has something to do with the
y1p_d = interp1(yL, ahy, y1_d);                                                 % definition of the angle of the arrowhead and
y2p_d = interp1(yL, ahy, y2_d);                                                 % the use of normalised coordinates

arhd = annotation('arrow', 'LineStyle','none',...                               % Creates the arrowhead annotaion and removing the line segment
    'HeadWidth',15, 'HeadStyle','vback2');
arhd.Units = 'normalized';                                                      % Normalaizing the units
arhd.Position = [x1p_d, y1p_d, x2p_d-x1p_d, y2p_d-y1p_d];                       % Defines position
arhd.Color = [.5 .5 .5];                                                        % Defines colour for arrowhead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend({'',string_array_1, '', string_array_2, '',...
%     string_array_3, "$\frac{d Ir}{d t}$"},...                                   % Creating a legend for the graphs
%     'Position', [.2375 .45 .1 .1],'Interpreter','latex', 'FontSize',15)
%--------------------------------------------------------------------------
annotation('textbox', [.15 .60 .1 .1], 'String',["Cherevko-", "Acidic"],...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText','on', 'FontSize',15);
annotation('textbox', [.70 .14 .1 .1], 'String',string_array_1,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Orange, 'Rotation',55);
annotation('textbox', [.27 .30 .1 .1], 'String',string_array_2,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Reddish_purple, 'Rotation',23);
annotation('textbox', [.25 .69 .1 .1], 'String',string_array_3,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',17, 'Color',Sky_blue, 'Rotation',5);
annotation('textbox', [.725 .46 .1 .1], 'String',string_array_4,...% Creating an annotation, textbox, with the rsquare value from the cfit
    'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor','none' ,'FontSize',20, 'Color',[.5 .5 .5]);
