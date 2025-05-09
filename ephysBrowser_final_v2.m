function ephysBrowser_final
% Final version with 6-panel figure and EC50 outputs to workspace (no normalization or Excel)

clear;
clc;
disp('Select .atf file');

[filename, pathname] = uigetfile({'*.atf'}, 'Select trace file');

pathFile = fullfile(pathname, filename);
a = importdata(pathFile, '\t', 11);
b = a.data(:,2:end); 
[r, ~] = size(b);
sr = 10000; 
t = (0:r-1)/sr;

%%%%%%%%% RAW PLOT FOR VISUAL TIMING %%%%%%%%%
figure;
plot(t, b);
title('Raw Data Preview - Use this to select baseline and drug timing');
numTracesPreview = size(b, 2);
text(0.95, 0.95, sprintf('Num Traces: %d', numTracesPreview), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (s)');
ylabel('pA');
yline(0, '--k');
grid on;

% User inputs
prompt = {'Baseline start (sec):','Baseline end (sec):',... 
    'Drug on (sec):','Drug off (sec):',...
    'Mutation name:','Drug + dosage:'};
dlgtitle = 'Experiment Info';
input_info = inputdlg(prompt, dlgtitle, [1 50]);

base1 = str2double(input_info{1});
base2 = str2double(input_info{2});
drugStart = str2double(input_info{3});
drugEnd = str2double(input_info{4});
mutationName = strrep(input_info{5}, '_', '\_');
drugInfo = strrep(input_info{6}, '_', '\_');

prompt_numSweeps = 'How many sweeps?';
numSweeps_str = inputdlg(prompt_numSweeps, 'Input', [1 50]);
numSweeps = str2double(numSweeps_str{1});
b = b(:, 1:numSweeps);

prompt_conc = sprintf('Enter %d concentrations (uM), separated by spaces:', numSweeps);
user_conc = inputdlg(prompt_conc, 'Input Concentrations', [1 70]);
conc = str2num(user_conc{1});
if length(conc) ~= numSweeps
    error('Mismatch between # of sweeps and # of concentrations.');
end

prompt_exclude = 'Sweeps to exclude (space-separated), or leave blank:';
user_exclude = inputdlg(prompt_exclude, 'Exclude Sweeps', [1 70]);
exclude_idx = [];
if ~isempty(user_exclude{1})
    exclude_idx = str2num(user_exclude{1});
end

keep_idx = setdiff(1:numSweeps, exclude_idx);
b = b(:, keep_idx);
conc = conc(keep_idx);
numTraces = length(keep_idx);

base = mean(b(base1*sr : base2*sr, :));
bs = b - repmat(base, r, 1);

ind1p = drugStart;  
ind2p = drugEnd;
ind1ss = drugEnd - 5;
ind2ss = drugEnd;

for i = 1:numTraces
    [~, I] = min(bs(ind1p*sr:ind2p*sr, i));
    medMin = median(bs((I-5)+ind1p*sr:I+ind1p*sr,i));
    minB(i) = medMin;
    minBss(i) = mean(bs(ind1ss*sr:ind2ss*sr, i));
end

Z = abs(minB);
Zss = abs(minBss);

keep_idx_fit = 2:numTraces; % skip 0 uM
c = conc(keep_idx_fit);
zr = Z(keep_idx_fit);
zssr = Zss(keep_idx_fit);

results_raw = ec50(c', zr');
resultsSS_raw = ec50(c', zssr');

% --- Plotting ---
hfig = figure('Position', [100, 100, 1200, 1000]);

subplot(3,2,1); 
plot(t, b); 
title(['Raw traces - ' mutationName]); 
xlabel('Time (s)'); ylabel('pA');
xline(drugStart, 'r--'); xline(drugEnd, 'r--');

subplot(3,2,2); 
plot(t, bs); 
title(['Baseline-subtracted - ' mutationName]); 
xlabel('Time (s)'); ylabel('pA');
xline(drugStart, 'r--'); xline(drugEnd, 'r--');

subplot(3,2,3); 
plot(log10(c), zr, 'ko-','MarkerFaceColor','k');
xlabel('log[Conc] (uM)'); ylabel('Peak pA');
title(['Peak Dose Response - ' drugInfo]);

subplot(3,2,4); hold on;
xFit = logspace(log10(min(c)), log10(max(c)), 100);
semilogx(c, zr, 'ko-','MarkerFaceColor','k');
semilogx(xFit, results_raw(1)+(results_raw(2)-results_raw(1))./(1+(results_raw(3)./xFit).^results_raw(4)), '-k');
title(sprintf('Peak EC50 = %.1f uM', results_raw(3)));
xlabel('Conc (uM)'); ylabel('Response (pA)');

subplot(3,2,5); 
plot(log10(c), zssr, 'ko-','MarkerFaceColor','k');
xlabel('log[Conc] (uM)'); ylabel('Steady-State pA');
title(['Steady-State Dose Response - ' drugInfo]);

subplot(3,2,6); hold on;
semilogx(c, zssr, 'ko-','MarkerFaceColor','k');
semilogx(xFit, resultsSS_raw(1)+(resultsSS_raw(2)-resultsSS_raw(1))./(1+(resultsSS_raw(3)./xFit).^resultsSS_raw(4)), '-k');
title(sprintf('SS EC50 = %.1f uM', resultsSS_raw(3)));
xlabel('Conc (uM)'); ylabel('Response (pA)');

sgtitle([mutationName ' - ' drugInfo]);

% --- Output tables to MATLAB workspace ---
T = table(conc(:), Z(:), Zss(:), ...
    'VariableNames', {'Conc_uM', 'Peak_pA', 'SS_pA'});

U = array2table([results_raw(:), resultsSS_raw(:)], ...
    'RowNames', {'Bottom', 'Top', 'EC50', 'Hill'}, ...
    'VariableNames', {'PeakFit', 'SSFit'});

assignin('base', 'T', T);
assignin('base', 'U', U);
disp('📊 Variables T (response) and U (fit parameters) are available in your workspace.');
