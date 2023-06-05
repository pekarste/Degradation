% This script will make a Pourbaix diagram of Ir in contact with water

%% Physical constants
R = 8.31446261815324;                                                       % J mol^-1 K^-1
F = 96485.3329;                                                             % A s mol^-1
T = 298.15;                                                                 % K

%% Thermodynamic values
E_1_0 = 0.926;                                                              % Standard electrode potential for the first equilibrium -  [V vs SHE] - Milazzo/Pourbaix
E_2_0 = 1.369;                                                              % Standard electrode potential for the second equilibrium - [V vs SHE] - Calculated by Materials values
Delta_G_3_0 = -134798.5;                                                    % Standard gibbs energy for reaction three at equilibrium - [J/mol]    - Calculated by Materials values
E_4_0 = 2.057;                                                              % Standard electrode potential for the fourth equilibrium - [V vs SHE} - Milazzo/Pourbaix
E_HER_0 = 0;                                                                % Standard electrode potential for HER                    - [V vs SHE] - def
E_OER_0 = 1.229;                                                            % Standard electrode potential for OER                    - [V vs SHE] - SI


% Number of electrons transferred
n_1 = 4;
n_2 = 2;
n_4 = 2;
n_HER = 2;
n_OER = 4;
%% Lines in Pourbaix diagram
pH = linspace(0,14,100);                                                    % Linspace of the pH range
pIr = [0 2 4 6];                                                                    % pIr = -log(a_IrO4^2-)

E_1 = E_1_0 - 4*((R*T)/(n_1*F))*log(10)*pH;                                 % The potential of reaction 1 as a function of pH
E_2 = E_2_0 - 2*((R*T)/(n_2*F))*log(10)*pH;                                 % The potential of reaction 2 as a function of pH

%pH_3 = -pIr/2 - (Delta_G_3_0)/(2*log(10)*R*T);                              % The equilibrium pH for reaction 3 --- NOT USED

E_HER = E_HER_0 - ((2*R*T)/(n_HER*F))*log(10)*pH;                           % The potential for the HER as a function of pH
E_OER = E_OER_0 - ((4*R*T)/(n_OER*F))*log(10)*pH;                           % The potential for the OER as a function of pH

%% Labels for the different regions and curves
Label_1 = "Ir(s)";
Label_2 = "IrO$_{2}$(s)";
Label_3 = "IrO$_{3}$(s)";
Label_4 = "IrO$^{2-}_{4}$(aq)";
Label_HER = "HER";
Label_OER = "OER";
Label_activities = ["pIr = 0" "pIr = 2" "pIr = 4" "pIr = 6"];

%% Finding where the lines intersect



f = figure('Name', 'Pourbaix Diagram of Ir');
fig_E_1 = plot(pH, E_1, 'Color','black');
hold on
for i = 1:length(pIr)
    E_4 = E_4_0 - ((R*T)/(n_4*F))*(log(10))*pIr(i) - 4*((R*T)/(n_4*F))*log(10)*pH; % The potential of reaction 4 as a function of pH
    pH_intersect_find = find(E_4 >= E_2);
    pH_intersect = pH_intersect_find(end);
    fig_E_2 = plot(pH(1:pH_intersect), E_2(1:pH_intersect), 'Color', 'black');
    fig_pH_3 = plot([pH(pH_intersect) pH(pH_intersect)], [E_2(pH_intersect) max(E_4)], 'Color', 'black');
    fig_E_4 = plot(pH(pH_intersect:end), E_4(pH_intersect:end), 'Color', 'black');
end
fig_HER = plot(pH, E_HER, 'LineStyle', '--', 'Color', [0.8500 0.3250 0.0980]);
fig_OER = plot(pH, E_OER, 'LineStyle', '--', 'Color', [0.8500 0.3250 0.0980]);

% Shade area between E_1 and pH-axis
%patch([pH fliplr(pH)], [E_1 ones(size(pH))*min(E_HER)], [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Shade area between E_HER and E_2
%patch([pH fliplr(pH(1:pH_intersect))], [E_HER fliplr(E_2(1:pH_intersect))], [0.4940 0.1840 0.5560], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Shade area between E_2 and pH_3
%patch([pH(pH_intersect) pH(pH_intersect:end)], [E_2(pH_intersect) fliplr(E_4(pH_intersect:end))], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Shade area between pH_3 and E_4
%patch([pH(pH_intersect:end) pH(pH_intersect)], [E_4(pH_intersect:end) fliplr(E_2(pH_intersect:end))], [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Shade area between E_4 and E_OER
%patch([pH(pH_intersect:end) pH], [E_4(pH_intersect:end) fliplr(E_OER)], [0.4940 0.1840 0.5560], 'FaceAlpha', 0.2, 'EdgeColor', 'none');


hold off
% Add labels with rotation
text(0.5, 0.3, Label_1, 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 15, 'Interpreter','latex');
text(0.5, 0.54, Label_2, 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 15, 'Interpreter','latex', 'Rotation', -12);
text(0.35, 0.85, Label_3, 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 15, 'Interpreter','latex');
text(0.9, 0.8, Label_4, 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 15, 'Interpreter','latex', 'Rotation', -90);
text(0.10, 0.225, Label_HER, 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 15, 'Interpreter','latex', 'Rotation', -12);
text(0.10, 0.675, Label_OER, 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 15, 'Interpreter','latex', 'Rotation', -12);
text(0.805, 0.9, Label_activities(1), 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12, 'Interpreter','latex', 'Rotation', -90);
text(0.735, 0.9, Label_activities(2), 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12, 'Interpreter','latex', 'Rotation', -90);
text(0.665, 0.9, Label_activities(3), 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12, 'Interpreter','latex', 'Rotation', -90);
text(0.595, 0.9, Label_activities(4), 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12, 'Interpreter','latex', 'Rotation', -90);

ax = gca;
ax.XAxis.FontSize = 12;                                                     % Changing the tick size on the x-axis
ax.YAxis.FontSize = 12;                                                     % Changing the tick size on the y-axis
xlabel('pH','Interpreter','latex', 'FontSize', 15)                          % Creating x-label
ylabel('Potential vs SHEH - E/[$V$]',...                                    % Creating y-label
    'Interpreter','latex', 'FontSize', 15)
xlim([pH(1) pH(end)])
ylim([min(E_HER) max(E_4)])