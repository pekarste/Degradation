function [t, gamma_theta] = diff_equation_solver_alkaline(time_array,data_type, curve, a_OH, Temperature, k_4_0_plus, theta_0)
%diff_equation_solver solves the differential equation describing the rate of change of site coverage
   


%% Solving the differential equation

t_span = [time_array(1) time_array(end)];                                   % The start en end time of the integration fot the solver
opts = odeset('RelTol',1e-10, 'AbsTol',1e-14);                              % Defining the difference between the solution and something, usded to not get so ragged curves 

[t,gamma_theta] = ode15s(@(t,gamma_theta)...                                % using ode15s to solve the diff equation described in diff_equation(...)
    diff_equation_alkaline(t,data_type, gamma_theta,curve,...
    a_OH, Temperature, k_4_0_plus), t_span, theta_0,opts);                  % This would then solve the differential equation based on a value for k_4_0_plus and gamma_theta_0


end