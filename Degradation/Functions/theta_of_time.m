function [thet_interpol] = theta_of_time(time_solver, theta_solver, time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
thet_interpol = interp1(time_solver,theta_solver,time);

end