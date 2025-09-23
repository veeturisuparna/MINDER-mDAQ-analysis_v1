clc;
clear all;

% Prompt user to select folder
folderPath = uigetdir('', 'Select the folder containing CSV files');
if folderPath == 0
    error('No folder selected. Exiting script.');
end

% Get list of all CSV files in selected folder
files = dir(fullfile(folderPath, '*.csv'));

% Sort files by name (timestamp order)
fileNames = {files.name};
[~, idx] = sort(fileNames);
files = files(idx);

% Initialize empty table for concatenated data
allData = [];

% Loop through each file and concatenate
for i = 1:length(files)
    filePath = fullfile(folderPath, files(i).name);
    currentData = readtable(filePath);
    
    if isempty(allData)
        allData = currentData;
    else
        allData = [allData; currentData];
    end
    
    fprintf('Processed file %d of %d: %s\n', i, length(files), files(i).name);
end

% Sampling rate
Fs = 100;
%%
%writetable(allData, 'C:\Users\Kaya\Desktop\MINDER data collection\Current data\concatenated_data.csv');
%%
% Sampling frequency
Fs = 100; 
t = (0:height(allData)-1)' / Fs; % Time vector in seconds

% Extract signals
ECG = allData.Var1;
EDA = allData.Var2;
PPG = allData.Var3;

% Create figure and layout
figure('Name','Signal Visualizer','NumberTitle','off','Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(3, 1);

% Plot ECG
ax1 = nexttile;
plot(t, ECG, 'b');
ylabel('ECG');
title('ECG Signal');
grid on;

% Plot EDA
ax2 = nexttile;
plot(t, EDA, 'r');
ylabel('EDA');
title('EDA Signal');
grid on;

% Plot PPG
ax3 = nexttile;
plot(t, PPG, 'g');
xlabel('Time (s)');
ylabel('PPG IR');
title('PPG IR Signal');
grid on;

% Link x-axes
linkaxes([ax1, ax2, ax3], 'x');

