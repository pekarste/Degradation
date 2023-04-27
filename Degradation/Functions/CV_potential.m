function potential = CV_potential(time, data_type)
%CV_potential will alculate the potential as a function of time
%   Mayrhofer is using linear sweep voltammetry, in order to express the
%   potential as a function of time in r_2, we need to adress that the
%   sweep rate changes direction at a certain time

%% Extracted Mayrhofer data that needs to be used

sweep_rate = 10*10^(-3);                                                    % Mayrhofer Sweep rate [V/s]
%t_max = 1630.9059738019;                                                   % Time at max potential - read off by plot digitizer [s]

t_0 = 16.826;                                                               % The lowest value is recorded at this time
E_0 = 0.040;                                                                % Lowest value for the potential - different from zero
t_min = t_0 + 2*((1.0 - E_0)/sweep_rate + (1.1 - E_0)/sweep_rate + ...
    (1.2 - E_0)/sweep_rate + (1.3 - E_0)/sweep_rate + ...
    (1.4 - E_0)/sweep_rate + (1.5 - E_0)/sweep_rate);% = 1468.826           % Time at the start of the positive sweep for our interest (E_0 --> 1.6)
t_max = t_min + (1.6 - E_0)/sweep_rate;% = 1624.826


%% Making the potential array as a functionof time


if data_type == "array"
    potential = ones(size(time))*E_0;                                       % It might be that this actually does not take in an array when used with ode15s...
    for  i = 1:length(time)                                 
        t = time(i);
    
        if t >= t_min && t < t_max
            potential(i) = E_0 + sweep_rate*(t-t_min);
        elseif t>t_max
            potential(i) = E_0 + sweep_rate*(2*t_max-t - t_min);
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