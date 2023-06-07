function potential = CV_potential_alkaline_2(time,data_type)
%CV_potential_alkaline will calculate the potential as a function of time
%   Schalenbach is using cyclic sweep voltammetry, in order to express the
%   potential as a function of time in r_2, we need to adress that the
%   sweep rate changes direction at a certain time

% used fot the second peaks from Schalenbach et al

%% Extracted Schalenbach data that needs to be used
% have used the Schalenbach_dissolution_CV_linear picture as source in
% plotdigitizer
sweep_rate = 2*10^(-3);                                                     % Schalenbach Sweep rate [V/s]
%t_max = 730;                                                               % Time at max potential - read off by plot digitizer [s]

t_0 = 5;                                                                    % The lowest value is recorded at this time
E_0 = 0.05;
E_max = 1.5;                                                                % Highest potential on the CV             
t_min = t_0 + 2.*(E_max - E_0)/sweep_rate;
t_max = t_min + (E_max - E_0)/sweep_rate;%


%% Making the potential array as a functionof time


if data_type == "array"
    potential = ones(size(time))*E_0;                                       % It might be that this actually does not take in an array when used with ode15s...
    for  i = 1:length(time)                                 
        t = time(i);
    
        if t >= t_min && t < t_max
            potential(i) = E_0 + sweep_rate*(t-t_min);
        elseif t>t_max
            potential(i) = E_0 + sweep_rate*(2*t_max-t - t_min);
        elseif t<t_min
            potential(i) = E_0 - sweep_rate*(t-t_min);
        end
    end

% This needs to take in a single value because ode15s evaluates one step at a time 
elseif data_type == "value"
    t = time;

    if t >= t_min && t < t_max
        potential = E_0 + sweep_rate*(t-t_min);
    elseif t>=t_max
        potential = E_0 + sweep_rate*(2*t_max-t - t_min);
    elseif t<t_min
        potential = E_0 - sweep_rate*(t - t_min);
    end

end