%% Plot feasibility test evaluations.
clear;
clc;
close all;
% To generate plot run test_feasibility first.

addpath(genpath('functions'));
result_path = '~/catkin_ws/build/mav_trajectory_generation_ros/';

%% Result plots.
% Get data.
result_sampling_01_file = 'result_sampling_01.txt';
result_sampling_05_file = 'result_sampling_05.txt';
result_sampling_10_file = 'result_sampling_10.txt';

result_recursive_01_file = 'result_recursive_01.txt';
result_recursive_05_file = 'result_recursive_05.txt';
result_recursive_10_file = 'result_recursive_10.txt';

file_id = fopen(strcat(result_path, result_sampling_01_file));
data_sampling_01 = fscanf(file_id, '%f');
fclose(file_id);

file_id = fopen(strcat(result_path, result_sampling_05_file));
data_sampling_05 = fscanf(file_id, '%f');
fclose(file_id);

file_id = fopen(strcat(result_path, result_sampling_10_file));
data_sampling_10 = fscanf(file_id, '%f');
fclose(file_id);

file_id = fopen(strcat(result_path, result_recursive_01_file));
data_recursive_01 = fscanf(file_id, '%f');
fclose(file_id);

file_id = fopen(strcat(result_path, result_recursive_05_file));
data_recursive_05 = fscanf(file_id, '%f');
fclose(file_id);

file_id = fopen(strcat(result_path, result_recursive_10_file));
data_recursive_10 = fscanf(file_id, '%f');
fclose(file_id);

% Create plots.
num_segments = length(data_recursive_10);
possible_outcomes = 0:8;
for i = 1 : length(possible_outcomes)
    possible_results{i} = feasibilityResultToString(possible_outcomes(i));
end
% Result summary
figure();
result_sampling_01_cnt = hist(data_sampling_01, possible_outcomes);
result_sampling_05_cnt = hist(data_sampling_05, possible_outcomes);
result_sampling_10_cnt = hist(data_sampling_10, possible_outcomes);
bar([result_sampling_01_cnt', result_sampling_05_cnt', result_sampling_10_cnt'])
set(gca, 'XTickLabel', possible_results)
legend('\Delta t = 0.01s', '\Delta t = 0.05s', '\Delta t = 0.10s')
ylim([0; num_segments])
grid
xlabel('Result')
ylabel('Count')
title('Results Sampling Test')

figure()
result_recursive_01_cnt = hist(data_recursive_01, possible_outcomes);
result_recursive_05_cnt = hist(data_recursive_05, possible_outcomes);
result_recursive_10_cnt = hist(data_recursive_10, possible_outcomes);
bar([result_recursive_01_cnt', result_recursive_05_cnt', result_recursive_10_cnt'])
set(gca, 'XTickLabel', possible_results)
legend('t_{min} = 0.01s', 't_{min} = 0.05s', 't_{min} = 0.10s')
ylim([0; num_segments])
grid
xlabel('Result')
ylabel('Count')
title('Results Recursive Test')

%% Timing plots.
times_file = 'feasibility_times.txt';
str = fileread(strcat(result_path, times_file));
str = strrep(str, '(', ' ');
str = strrep(str, '+-', '');
str = strrep(str, ')', '');
str = strrep(str, '[', '');
str = strrep(str, ']', '');
str = strrep(str, ',', ' ');
times_data = textscan(str, '%s %d %f %f %f %f %f', 'HeaderLines', 2);

mean = times_data{4};
stddev = times_data{5};

min_time = min(mean) - max(stddev);
max_time = max(mean) + max(stddev);

% Computation time mean and standard deviation.
% Sampling
figure()
sampling_interval = [0.01 0.05 0.1];
errorbar(sampling_interval, mean(4:end), stddev(4:end))
xlabel('Sampling Interval \Delta t');
ylabel('Time [s]');
ylim([min_time max_time]);
title('Sampling: mean and stdev of computation times')

% Recursive
figure()
min_section_time = [0.01 0.05 0.1];
errorbar(min_section_time, mean(1:3), stddev(1:3))
xlabel('Minimum section time t_{min}');
ylabel('Time [s]');
ylim([min_time max_time]);
title('Recursive: mean and stdev of computation times')


