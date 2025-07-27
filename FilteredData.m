% Filter ATF Data: Lowpass and Notch Options for All Traces

clear; clc;


% --- Prompt for Cell ID, Mutation, Drug, number of sweeps, and concentrations ---
prompt = {'Enter Cell ID:', 'Enter Mutation:', 'Enter Drug:', 'Number of sweeps:', 'Concentrations (comma-separated, lowest to highest):'};
dlg_title = 'Input Cell Info';
num_lines = 1;
defaultans = {'', '', '', '5', '1.111,3.333,11.111,33.33,100'};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
if isempty(answer)
    error('No input provided. Script terminated.');
end
cell_id = answer{1};
mutation = answer{2};
drug = answer{3};
num_sweeps = str2double(answer{4});
concs = str2num(answer{5}); %#ok<ST2NM>
if length(concs) ~= num_sweeps
    error('Number of concentrations must match number of sweeps.');
end

disp('Select .atf file');
[filename, pathname] = uigetfile({'*.atf'}, 'Select ATF trace file');
pathFile = fullfile(pathname, filename);

a = importdata(pathFile, '\t', 11);
data = a.data(:,2:(1+num_sweeps));  % dynamically select sweeps
[r, c] = size(data);
sr = 10000;              % sampling rate (Hz)
t = (0:r-1)/sr;          % time vector
Fs = sr;


Fs = sr;  % 'sr' must be defined earlier in your script

% --- Lowpass filter setup ---
lowpass_cutoff = 2;  % Hz
[b_lp, a_lp] = butter(4, lowpass_cutoff / (Fs/2), 'low');

% --- Notch filter setup (manual, compatible with older MATLAB) ---
notch_freq = 3.5;     % Center frequency of the notch (Hz)
notch_bw   = 1.5;     % Bandwidth of the notch (Hz)

f0 = notch_freq;               % Notch center frequency
bw_hz = notch_bw;              % Bandwidth in Hz
Q = f0 / bw_hz;                % Quality factor
w0 = f0 / (Fs/2);              % Normalized center frequency (0–1)
bw = w0 / Q;                   % Normalized bandwidth
[b_notch, a_notch] = butter(2, [w0 - bw/2, w0 + bw/2], 'stop');  % 2nd-order band-stop


filtered_data = zeros(size(data));

for i = 1:size(data, 2)
    raw_trace = data(:, i);

    % Apply lowpass filter
    trace_lp = filtfilt(b_lp, a_lp, raw_trace);

    % Apply notch filter
    trace_filtered = filtfilt(b_notch, a_notch, trace_lp);

    % Store filtered trace
    filtered_data(:, i) = trace_filtered;
end



% --- Dynamic titles based on concentrations ---
titles = cell(1, num_sweeps);
for i = 1:num_sweeps
    titles{i} = sprintf('METH %.3f \muM', concs(i));
end

figure;
for i = 1:num_sweeps
    subplot(num_sweeps, 1, i);
    plot(t, data(:, i), 'Color', [0.6 0.6 0.6]); hold on;
    plot(t, filtered_data(:, i), 'r');
    % Vertical lines for drug on/off
    xline(35.5, '--k', 'METH on', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
    xline(45.5, '--k', 'METH off', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    % Add cocaine start and end lines
    xline(6, '--b', 'COC on', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
    xline(21, '--b', 'COC off', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    title(titles{i}, 'Interpreter', 'tex');
    ylabel('pA');
    if i == num_sweeps
        xlabel('Time (s)');
    end
    legend('Raw', 'Filtered', 'Location', 'northeastoutside');
end
% Set figure title as Cell ID, Mutation, and Drug
fig_title = sprintf('%s | %s | %s', cell_id, mutation, drug);
sgtitle(fig_title, 'FontWeight', 'bold');

% --- Save figure in same directory as ATF file ---
[~, baseName, ~] = fileparts(filename);
figSavePath = fullfile(pathname, [baseName '_filtered_plot.png']);
saveas(gcf, figSavePath);

% Display the saved figure path in the command window
fprintf('Filtered figure saved to: %s\n', figSavePath);

