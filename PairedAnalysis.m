%% === SETTINGS & FILE LOADING ===
samplingRate = 10000;  % Hz
[file, path] = uigetfile('*.atf', 'Select ATF file', 'C:\Users\jew052\Desktop\Jesse2');
if isequal(file, 0)
    error('No file selected');
end
fullpath = fullfile(path, file);
T = readtable(fullpath, 'FileType', 'text', 'Delimiter', '\t', 'HeaderLines', 10);

%% === METADATA INPUT ===
% Prompt for genotype using a GUI dialog
genotype = questdlg('Select Genotype:', 'Genotype Selection', 'WT', 'DKO', 'WT');
if isempty(genotype)
    error('No genotype selected');
end
pairID   = input('Enter Pair ID (e.g., R1, R2): ', 's');
cellID_X = input('Enter Cell X ID (e.g., A, B, etc.): ', 's');
cellID_Y = input('Enter Cell Y ID (e.g., A, B, etc.): ', 's');

%% === DATA WINDOW & COLUMN SETUP ===
timeVec = T{:,1};  % First column is time (in seconds)

% Define analysis windows (in seconds)
y_base_idx = find(timeVec >= 0.331   & timeVec <= 0.531); %changed baseline start from 0.031 to 0.331
y_step_idx = find(timeVec >= 0.681   & timeVec <= 0.731);
x_base_idx = find(timeVec >= 1.031   & timeVec <= 1.231); %changed baseline from 0.732 to 1.031
x_step_idx = find(timeVec >= 1.381   & timeVec <= 1.431); 

% Get sweep columns (interleaved format)
allCols = T.Properties.VariableNames(2:end);  % Skip time
cellX_cols = allCols(1:2:end);  % Columns 2, 4, ..., 48
cellY_cols = allCols(2:2:end);  % Columns 3, 5, ..., 49
nSweeps = length(cellX_cols);  % Should be 24
deltaV = -115:10:115;
assert(nSweeps == length(deltaV), 'Mismatch between sweeps and voltage steps');

%% === ANALYSIS & OUTPUT ===
deltaI_X = zeros(1, nSweeps);
deltaI_Y = zeros(1, nSweeps);
for i = 1:nSweeps
    traceX = T{:, cellX_cols{i}};
    traceY = T{:, cellY_cols{i}};
    baselineX = mean(traceX(x_base_idx));
    steadyX   = mean(traceX(x_step_idx));
    deltaI_X(i) = steadyX - baselineX;
    baselineY = mean(traceY(y_base_idx));
    steadyY   = mean(traceY(y_step_idx));
    deltaI_Y(i) = steadyY - baselineY;
end

% Fit lines & compute R²
pX = polyfit(deltaV, deltaI_X, 1);
fitX = polyval(pX, deltaV);
R2_X = 1 - sum((deltaI_X - fitX).^2) / sum((deltaI_X - mean(deltaI_X)).^2);
pY = polyfit(deltaV, deltaI_Y, 1);
fitY = polyval(pY, deltaV);
R2_Y = 1 - sum((deltaI_Y - fitY).^2) / sum((deltaI_Y - mean(deltaI_Y)).^2);

% Display results
idx_minus115 = find(deltaV == -115);
fprintf('\n✅ Analysis complete for %s (%s)\n', file, pairID);
fprintf('Cell X: Steady-State ΔI at -115 mV = %.3f pA, R² = %.3f\n', deltaI_X(idx_minus115), R2_X);
fprintf('Cell Y: Steady-State ΔI at -115 mV = %.3f pA, R² = %.3f\n', deltaI_Y(idx_minus115), R2_Y);

% Save summary to Excel
if strcmpi(genotype, 'WT')
    summaryDir = 'C:\Users\jew052\Desktop\Jesse2\WT';
    if ~exist(summaryDir, 'dir')
        mkdir(summaryDir);
    end
    summaryFile = fullfile(summaryDir, 'paired_cells_summary_WT.xlsx');
elseif strcmpi(genotype, 'DKO')
    summaryDir = 'C:\Users\jew052\Desktop\Jesse2\DKO';
    if ~exist(summaryDir, 'dir')
        mkdir(summaryDir);
    end
    summaryFile = fullfile(summaryDir, 'paired_cells_summary_DKO.xlsx');
else
    summaryFile = fullfile(path, 'paired_cells_summary_OTHER.xlsx');
end
rowX = table({genotype}, {pairID}, {file}, {fullpath}, {cellID_X}, deltaI_X(idx_minus115), R2_X, ...
    'VariableNames', {'Genotype', 'Pair', 'File', 'FilePath', 'CellID', 'MaxDeltaI', 'R2'});
rowY = table({genotype}, {pairID}, {file}, {fullpath}, {cellID_Y}, deltaI_Y(idx_minus115), R2_Y, ...
    'VariableNames', {'Genotype', 'Pair', 'File', 'FilePath', 'CellID', 'MaxDeltaI', 'R2'});
if isfile(summaryFile)
    existing = readtable(summaryFile);
    summary = [existing; rowX; rowY];
else
    summary = [rowX; rowY];
end
writetable(summary, summaryFile);

% Save full ΔI to CSV
T_DeltaI = table(deltaV', deltaI_X', deltaI_Y', 'VariableNames', {'Voltage_mV', ['DeltaI_' cellID_X], ['DeltaI_' cellID_Y]});
safePairID = regexprep(pairID, '[^a-zA-Z0-9_]', '_');
writetable(T_DeltaI, fullfile(path, [safePairID '_DeltaI.csv']));

%% === PLOTTING ===
sweepNum = 1;  % <-- Change this to view a different sweep
traceX = T{:, cellX_cols{sweepNum}};
traceY = T{:, cellY_cols{sweepNum}};
figure;
subplot(2,1,1)
plot(timeVec, traceX, 'b'); hold on;
yline(mean(traceX(x_base_idx)), '--k', 'Baseline');
yline(mean(traceX(x_step_idx)), '--r', 'Steady-State');
area(timeVec(x_base_idx), traceX(x_base_idx), 'FaceAlpha', 0.2, 'FaceColor', 'k', 'EdgeColor', 'none');
area(timeVec(x_step_idx), traceX(x_step_idx), 'FaceAlpha', 0.2, 'FaceColor', 'r', 'EdgeColor', 'none');
title(sprintf('Cell X – Sweep %d', sweepNum));
xlabel('Time (s)');
ylabel('Current (pA)');
subplot(2,1,2)
plot(timeVec, traceY, 'g'); hold on;
yline(mean(traceY(y_base_idx)), '--k', 'Baseline');
yline(mean(traceY(y_step_idx)), '--r', 'Steady-State');
area(timeVec(y_base_idx), traceY(y_base_idx), 'FaceAlpha', 0.2, 'FaceColor', 'k', 'EdgeColor', 'none');
area(timeVec(y_step_idx), traceY(y_step_idx), 'FaceAlpha', 0.2, 'FaceColor', 'r', 'EdgeColor', 'none');
title(sprintf('Cell Y – Sweep %d', sweepNum));
xlabel('Time (s)');
ylabel('Current (pA)');
