%% Plot feasibility test evaluations.
% To generate plot run test_feasibility first.

addpath(genpath('functions'));
result_path = '~/catkin_ws/build/mav_trajectory_generation_ros/';

file_id = fopen(strcat(result_path, 'result_recursive_10.txt'));
data = fscanf(file_id, '%f');
fclose(file_id);

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

%% Timing plots.