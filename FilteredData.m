% Filter ATF Data: Lowpass and Notch Options for All Traces

clear; clc;
disp('Select .atf file');
[filename, pathname] = uigetfile({'*.atf'}, 'Select ATF trace file');
pathFile = fullfile(pathname, filename);


a = importdata(pathFile, '\t', 11);
data = a.data(:,2:end);  % exclude time column
[r, c] = size(data);
sr = 10000;              % sampling rate (Hz)
t = (0:r-1)/sr;          % time vector
Fs = sr;

%%
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


%% --- Plot all 4 raw vs filtered traces with custom titles and annotations
titles = {
    'METH 11.11 \muM', 
    'METH 33.33 \muM', 
    'METH 100 \muM', 
    'METH 300 \muM'
};

figure;
for i = 1:4
    subplot(4, 1, i);
    plot(t, data(:, i), 'Color', [0.6 0.6 0.6]); hold on;
    plot(t, filtered_data(:, i), 'r');
    
    % Vertical lines for drug on/off
    xline(5.5, '--k', 'drug on', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
    xline(20.5, '--k', 'drug off', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    
    title(titles{i}, 'Interpreter', 'tex');
    ylabel('pA');
    if i == 4
        xlabel('Time (s)');
    end
    legend('Raw', 'Filtered');
end
