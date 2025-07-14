%% === SETTINGS ===
samplingRate = 10000;  % Hz


[file, path] = uigetfile('*.atf', 'Select ATF file');
if isequal(file, 0)
    error('No file selected');
end
fullpath = fullfile(path, file);
T = readtable(fullpath, 'FileType', 'text', 'Delimiter', '\t', 'HeaderLines', 10);

%% === GET METADATA ===
genotype = input('Enter Genotype (e.g., WT, DKO): ', 's');
pairID   = input('Enter Pair ID (e.g., R1, R2): ', 's');



%% === SET TIME WINDOWS (in ms) ===
timeVec = T{:,1};  % First column is time (in seconds)

y_base_idx = find(timeVec >= 0.031   & timeVec <= 0.531);
y_step_idx = find(timeVec >= 0.681   & timeVec <= 0.731);
x_base_idx = find(timeVec >= 0.731   & timeVec <= 1.231);
x_step_idx = find(timeVec >= 1.381   & timeVec <= 1.431);




%% === GET SWEEP COLUMNS (INTERLEAVED FORMAT) ===
% Column 1 = Time, columns 2–49 = alternating X1, Y1, X2, Y2, ...
allCols = T.Properties.VariableNames(2:end);  % Skip time

cellX_cols = allCols(1:2:end);  % Columns 2, 4, ..., 48
cellY_cols = allCols(2:2:end);  % Columns 3, 5, ..., 49

nSweeps = length(cellX_cols);  % Should be 24
deltaV = -115:10:115;
assert(nSweeps == length(deltaV), 'Mismatch between sweeps and voltage steps');

%% === EXTRACT ΔI FOR EACH SWEEP ===
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

%% === FIT LINES & COMPUTE R² ===
pX = polyfit(deltaV, deltaI_X, 1);
fitX = polyval(pX, deltaV);
R2_X = 1 - sum((deltaI_X - fitX).^2) / sum((deltaI_X - mean(deltaI_X)).^2);

pY = polyfit(deltaV, deltaI_Y, 1);
fitY = polyval(pY, deltaV);
R2_Y = 1 - sum((deltaI_Y - fitY).^2) / sum((deltaI_Y - mean(deltaI_Y)).^2);

% Display results
fprintf('\n✅ Analysis complete for %s (%s)\n', file, pairID);
fprintf('Cell X: Max ΔI = %.3f nA, R² = %.3f\n', min(deltaI_X), R2_X);
fprintf('Cell Y: Max ΔI = %.3f nA, R² = %.3f\n', min(deltaI_Y), R2_Y);

%% === SAVE SUMMARY TO EXCEL ===
summaryFile = fullfile(path, 'paired_cells_summary.xlsx');
rowX = table({genotype}, {pairID}, {file}, {'Cell X'}, min(deltaI_X), R2_X, ...
    'VariableNames', {'Genotype', 'Pair', 'File', 'CellID', 'MaxDeltaI', 'R2'});
rowY = table({genotype}, {pairID}, {file}, {'Cell Y'}, min(deltaI_Y), R2_Y, ...
    'VariableNames', {'Genotype', 'Pair', 'File', 'CellID', 'MaxDeltaI', 'R2'});

if isfile(summaryFile)
    existing = readtable(summaryFile);
    summary = [existing; rowX; rowY];
else
    summary = [rowX; rowY];
end

writetable(summary, summaryFile);

%% === SAVE FULL ΔI TO CSV ===
T_X = table(deltaV', deltaI_X', 'VariableNames', {'Voltage_mV', 'SteadyState_DeltaI_nA'});
T_Y = table(deltaV', deltaI_Y', 'VariableNames', {'Voltage_mV', 'SteadyState_DeltaI_nA'});

safePairID = regexprep(pairID, '[^a-zA-Z0-9_]', '_');
writetable(T_X, fullfile(path, [safePairID '_CellX_DeltaI.csv']));
writetable(T_Y, fullfile(path, [safePairID '_CellY_DeltaI.csv']));


%% === PLOT TRACE WITH BASELINE AND STEADY-STATE WINDOWS ===
sweepNum = 1;  % <-- Change this to view a different sweep

traceX = T{:, cellX_cols{sweepNum}};
traceY = T{:, cellY_cols{sweepNum}};

figure;
subplot(2,1,1)
plot(timeVec, traceX, 'b'); hold on;
yline(mean(traceX(x_base_idx)), '--k', 'Baseline');
yline(mean(traceX(x_step_idx)), '--r', 'Steady-State'); % <-- Add missing steady-state for Cell X
area(timeVec(x_base_idx), traceX(x_base_idx), 'FaceAlpha', 0.2, 'FaceColor', 'k', 'EdgeColor', 'none');
area(timeVec(x_step_idx), traceX(x_step_idx), 'FaceAlpha', 0.2, 'FaceColor', 'r', 'EdgeColor', 'none');
title(sprintf('Cell X – Sweep %d', sweepNum));
xlabel('Time (s)');
ylabel('Current (pA)');

subplot(2,1,2)
plot(timeVec, traceY, 'g'); hold on;
yline(mean(traceY(y_base_idx)), '--k', 'Baseline'); % <-- Add missing baseline for Cell Y
yline(mean(traceY(y_step_idx)), '--r', 'Steady-State');
area(timeVec(y_base_idx), traceY(y_base_idx), 'FaceAlpha', 0.2, 'FaceColor', 'k', 'EdgeColor', 'none');
area(timeVec(y_step_idx), traceY(y_step_idx), 'FaceAlpha', 0.2, 'FaceColor', 'r', 'EdgeColor', 'none');
title(sprintf('Cell Y – Sweep %d', sweepNum));
xlabel('Time (s)');
ylabel('Current (pA)');
