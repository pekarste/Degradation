% This script is just made to run k_3_0_extraction

k_4_0_plus_array = [10^(-1) 10^(-2) 10^(-3) 10^(-4)];                       % An inital array of k_4_0_plus values

k_3_0_plus_scohy_array = ones(size(k_4_0_plus_array));                      % An array where the values of k_3_0_plus based on Scohy will be stored
k_3_0_plus_damjanovic_array = ones(size(k_4_0_plus_array));                 % An array where the values of k_3_0_plus based on Damjanovic will be stored
k_3_0_plus_damjanovic_log_array = ones(size(k_4_0_plus_array));             % An array where the values of k_3_0_plus based on Damjanovic_log will be stored
k_3_0_plus_cherevko_array = ones(size(k_4_0_plus_array));

for i = 1:length(k_4_0_plus_array)                                          % Looping through the number of elements corresponding to the lenght of k_4_0_plus_array

    [k_3_0_plus_scohy, k_3_0_plus_damjanovic, k_3_0_plus_damjanovic_log, k_3_0_plus_cherevko] = ...
        k_3_0_plus_extraction(k_4_0_plus_array(i));                         % Calculating the k_3_0_plus values based on the htree parameter fits for ever value of k_4_0_plus

    k_3_0_plus_scohy_array(i) = k_3_0_plus_scohy;                           % Appending the k_3_0_plus value to its correct array
    k_3_0_plus_damjanovic_array(i) = k_3_0_plus_damjanovic;                 % Appending the k_3_0_plus value to its correct array
    k_3_0_plus_damjanovic_log_array(i) = k_3_0_plus_damjanovic_log;         % Appending the k_3_0_plus value to its correct array
    k_3_0_plus_cherevko_array(i) = k_3_0_plus_cherevko;
end