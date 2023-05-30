function [t, theta] = diff_equation_solver_alkaline(time_array,data_type, curve, a_OH, Temperature, k_4_0_plus, theta_0, frac)
%diff_equation_solver solves the differential equation describing the rate of change of site coverage
   


%% Solving the differential equation

t_span = [time_array(1) time_array(end)];                                   % The start en end time of the integration fot the solver

[t,theta] = ode15s(@(t,theta)...                                            % using ode15s to solve the diff equation described in diff_equation(...)
    diff_equation_alkaline(t,data_type, theta,curve,...
    a_OH, Temperature, k_4_0_plus, frac), t_span, theta_0);                       % This would then solve the differential equation based on a value for k_4_0_plus and gamma_theta_0


end