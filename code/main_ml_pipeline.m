clear; clc; close all;

%% --- 1. Flow Parameters ---
lambda = 0.428;
x0     = 0.0;
n_std  = 1.0;

%% =========================================================
%  PART A — LOAD ORIGINAL RAW DATA
%% =========================================================

fprintf('=================================================\n')
fprintf(' PART A: Loading Original Raw Data\n')
fprintf('=================================================\n')

% Read raw header to get token positions
fid = fopen('isothermal_all_triplets.csv', 'rb');
content = fread(fid, Inf, '*char')';
fclose(fid);

nl = find(content == newline, 1);
header = content(1:nl-1);
if header(end) == char(13); header = header(1:end-1); end

tokens = strsplit(header, ',');
fprintf('Header has %d tokens\n', length(tokens))

% Merge broken column tokens into clean names
merged = {};
for i = 1:length(tokens)
    t = strtrim(tokens{i});
    if contains(lower(t), 'accept')
        merged{end+1} = 'label';
    elseif contains(lower(t), 'reject') || contains(lower(t), 'unlabelled')
        % skip these broken tokens
    elseif isempty(t)
        % skip empty tokens
    else
        merged{end+1} = t;
    end
end

fprintf('Merged to %d column names\n', length(merged))

% Load table (MATLAB reads data correctly as Var1-Var33)
raw = readtable('isothermal_all_triplets.csv');
fprintf('Table has %d columns, %d rows\n', width(raw), height(raw))

% Assign correct column names by position
if width(raw) == length(merged)
    raw.Properties.VariableNames = merged;
    fprintf('Column names assigned successfully\n')
else
    fprintf('Mismatch: %d names vs %d columns\n', length(merged), width(raw))
    fprintf('Column names found:\n')
    for i = 1:length(merged); fprintf('  %d: %s\n', i, merged{i}); end
    error('Fix the column name count above')
end

% Extract key columns
l_norm_raw = raw.l_over_xminusx0;
b_a_raw    = raw.b_over_a;
x_m_raw    = 0.5 * (raw.x_upstream + raw.x_downstream);

fprintf('Mean l/x (raw)        : %.4f\n', mean(l_norm_raw))
fprintf('Std  l/x (raw)        : %.4f\n', std(l_norm_raw))
fprintf('Mean b/a (raw)        : %.4f\n', mean(b_a_raw))

%% =========================================================
%  PART B — APPROACH 1: STATISTICAL THRESHOLD CLEANING
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' PART B: Statistical Threshold Cleaning\n')
fprintf('=================================================\n')

mu    = mean(l_norm_raw);
sig   = std(l_norm_raw);
lower = mu - n_std * sig;
upper = mu + n_std * sig;

fprintf('Mean l/x              : %.4f\n', mu)
fprintf('Std  l/x              : %.4f\n', sig)
fprintf('Keeping range         : %.4f to %.4f\n', lower, upper)

mask_stat      = l_norm_raw >= lower & l_norm_raw <= upper;
raw_stat_clean = raw(mask_stat, :);
l_norm_stat    = raw_stat_clean.l_over_xminusx0;
b_a_stat       = raw_stat_clean.b_over_a;

fprintf('Original rows         : %d\n', height(raw))
fprintf('After stat cleaning   : %d\n', height(raw_stat_clean))
fprintf('Removed               : %d\n', height(raw) - height(raw_stat_clean))
fprintf('Mean l/x (stat clean) : %.4f\n', mean(l_norm_stat))
fprintf('Mean b/a (stat clean) : %.4f\n', mean(b_a_stat))

writetable(raw_stat_clean, 'stat_clean_triplets.csv')
fprintf('Saved stat_clean_triplets.csv\n')

%% =========================================================
%  PART C — APPROACH 2: ML CLEANING (Random Forest)
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' PART C: ML Cleaning (Random Forest)\n')
fprintf('=================================================\n')

ml_clean  = readtable('clean_triplets.csv');
l_norm_ml = ml_clean.l_over_xminusx0;
b_a_ml    = ml_clean.b_over_a;
x_m_ml    = 0.5 * (ml_clean.x_upstream + ml_clean.x_downstream);

fprintf('After ML cleaning     : %d\n', height(ml_clean))
fprintf('Removed               : %d\n', height(raw) - height(ml_clean))
fprintf('Mean l/x (ML clean)   : %.4f\n', mean(l_norm_ml))
fprintf('Mean b/a (ML clean)   : %.4f\n', mean(b_a_ml))

%% =========================================================
%  PART D — COMPARE PDFs
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' PART D: PDF Comparison\n')
fprintf('=================================================\n')

% --- Figure 1: l/(x-x0) comparison ---
figure('Position', [100 100 1500 500]);

subplot(1,3,1)
histogram(l_norm_raw, 80, 'Normalization', 'pdf', ...
          'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(l_norm_raw), 'r-', 'LineWidth', 2, ...
      'Label', sprintf('Mean = %.4f', mean(l_norm_raw)), ...
      'LabelVerticalAlignment', 'bottom')
xlabel('l / (x_m - x_0)', 'FontSize', 12)
ylabel('Probability Density', 'FontSize', 12)
title(sprintf('Original Data\n(n = %d)', height(raw)), 'FontSize', 13)
grid on; box on

subplot(1,3,2)
histogram(l_norm_stat, 80, 'Normalization', 'pdf', ...
          'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(l_norm_stat), 'r-', 'LineWidth', 2, ...
      'Label', sprintf('Mean = %.4f', mean(l_norm_stat)), ...
      'LabelVerticalAlignment', 'bottom')
xlabel('l / (x_m - x_0)', 'FontSize', 12)
ylabel('Probability Density', 'FontSize', 12)
title(sprintf('Statistical Cleaning (±%.0f\\sigma)\n(n = %d)', ...
      n_std, height(raw_stat_clean)), 'FontSize', 13)
grid on; box on

subplot(1,3,3)
histogram(l_norm_ml, 80, 'Normalization', 'pdf', ...
          'FaceColor', [0.1 0.7 0.3], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(l_norm_ml), 'r-', 'LineWidth', 2, ...
      'Label', sprintf('Mean = %.4f', mean(l_norm_ml)), ...
      'LabelVerticalAlignment', 'bottom')
xlabel('l / (x_m - x_0)', 'FontSize', 12)
ylabel('Probability Density', 'FontSize', 12)
title(sprintf('ML Cleaning (Random Forest)\n(n = %d)', height(ml_clean)), 'FontSize', 13)
grid on; box on

saveas(gcf, 'figure1_comparison_lnorm.png')
fprintf('Figure 1 saved: figure1_comparison_lnorm.png\n')

% --- Figure 2: b/a comparison ---
figure('Position', [100 100 1500 500]);

subplot(1,3,1)
histogram(b_a_raw, 80, 'Normalization', 'pdf', ...
          'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(b_a_raw), 'r-', 'LineWidth', 2, ...
      'Label', sprintf('Mean = %.4f', mean(b_a_raw)), ...
      'LabelVerticalAlignment', 'bottom')
xline(1.0, 'g--', 'LineWidth', 2, 'Label', 'b/a = 1')
xlabel('b / a', 'FontSize', 12)
ylabel('Probability Density', 'FontSize', 12)
title(sprintf('Original Data\n(n = %d)', height(raw)), 'FontSize', 13)
xlim([0 8]); grid on; box on

subplot(1,3,2)
histogram(b_a_stat, 80, 'Normalization', 'pdf', ...
          'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(b_a_stat), 'r-', 'LineWidth', 2, ...
      'Label', sprintf('Mean = %.4f', mean(b_a_stat)), ...
      'LabelVerticalAlignment', 'bottom')
xline(1.0, 'g--', 'LineWidth', 2, 'Label', 'b/a = 1')
xlabel('b / a', 'FontSize', 12)
ylabel('Probability Density', 'FontSize', 12)
title(sprintf('Statistical Cleaning (±%.0f\\sigma)\n(n = %d)', ...
      n_std, height(raw_stat_clean)), 'FontSize', 13)
xlim([0 8]); grid on; box on

subplot(1,3,3)
histogram(b_a_ml, 80, 'Normalization', 'pdf', ...
          'FaceColor', [0.1 0.7 0.3], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(b_a_ml), 'r-', 'LineWidth', 2, ...
      'Label', sprintf('Mean = %.4f', mean(b_a_ml)), ...
      'LabelVerticalAlignment', 'bottom')
xline(1.0, 'g--', 'LineWidth', 2, 'Label', 'b/a = 1')
xlabel('b / a', 'FontSize', 12)
ylabel('Probability Density', 'FontSize', 12)
title(sprintf('ML Cleaning (Random Forest)\n(n = %d)', height(ml_clean)), 'FontSize', 13)
xlim([0 8]); grid on; box on

saveas(gcf, 'figure2_comparison_ba.png')
fprintf('Figure 2 saved: figure2_comparison_ba.png\n')

% --- Figure 3: Joint PDFs ---
figure('Position', [100 100 1500 500]);

subplot(1,3,1)
histogram2(l_norm_raw, b_a_raw, 60, ...
           'Normalization', 'count', 'DisplayStyle', 'tile', ...
           'ShowEmptyBins', 'off')
colorbar; colormap('parula')
xlabel('l / (x_m - x_0)', 'FontSize', 12)
ylabel('b / a', 'FontSize', 12)
title(sprintf('Joint PDF — Original\n(n = %d)', height(raw)), 'FontSize', 13)
xlim([0 0.6]); ylim([0 6]); grid on; box on

subplot(1,3,2)
histogram2(l_norm_stat, b_a_stat, 60, ...
           'Normalization', 'count', 'DisplayStyle', 'tile', ...
           'ShowEmptyBins', 'off')
colorbar; colormap('parula')
xlabel('l / (x_m - x_0)', 'FontSize', 12)
ylabel('b / a', 'FontSize', 12)
title(sprintf('Joint PDF — Stat Cleaned\n(n = %d)', ...
      height(raw_stat_clean)), 'FontSize', 13)
xlim([0 0.6]); ylim([0 6]); grid on; box on

subplot(1,3,3)
histogram2(l_norm_ml, b_a_ml, 60, ...
           'Normalization', 'count', 'DisplayStyle', 'tile', ...
           'ShowEmptyBins', 'off')
colorbar; colormap('parula')
xlabel('l / (x_m - x_0)', 'FontSize', 12)
ylabel('b / a', 'FontSize', 12)
title(sprintf('Joint PDF — ML Cleaned\n(n = %d)', height(ml_clean)), 'FontSize', 13)
xlim([0 0.6]); ylim([0 6]); grid on; box on

saveas(gcf, 'figure3_comparison_joint.png')
fprintf('Figure 3 saved: figure3_comparison_joint.png\n')

%% =========================================================
%  PART E — DOWNSTREAM TRENDS
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' PART E: Downstream Trends\n')
fprintf('=================================================\n')

edges    = linspace(min(x_m_ml), max(x_m_ml), 11);
bin_mids = 0.5 * (edges(1:end-1) + edges(2:end));

mean_l_norm_bin   = zeros(1,10);
mean_b_over_a_bin = zeros(1,10);
mean_l_bin        = zeros(1,10);

for i = 1:10
    mask = x_m_ml >= edges(i) & x_m_ml < edges(i+1);
    mean_l_norm_bin(i)   = mean(l_norm_ml(mask));
    mean_b_over_a_bin(i) = mean(b_a_ml(mask));
    mean_l_bin(i)        = mean(ml_clean.l(mask));
end

figure('Position', [100 100 1400 450]);

subplot(1,3,1)
plot(bin_mids, mean_l_norm_bin, 'o-', ...
     'Color', [0.2 0.4 0.8], 'LineWidth', 2, ...
     'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.4 0.8])
hold on
yline(lambda * 0.518, 'r--', 'LineWidth', 1.5, ...
      'Label', sprintf('k\\lambda = %.4f', lambda * 0.518))
xlabel('x_m (downstream position)', 'FontSize', 12)
ylabel('Mean l / (x_m - x_0)', 'FontSize', 12)
title('Normalised Spacing vs Position', 'FontSize', 13)
grid on; box on

subplot(1,3,2)
plot(bin_mids, mean_b_over_a_bin, 'o-', ...
     'Color', [0.9 0.5 0.1], 'LineWidth', 2, ...
     'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.5 0.1])
hold on
yline(1.0, 'g--', 'LineWidth', 2, 'Label', 'b/a = 1 (symmetric)')
yline(mean(b_a_ml), 'r--', 'LineWidth', 1.5, ...
      'Label', sprintf('Mean = %.4f', mean(b_a_ml)))
xlabel('x_m (downstream position)', 'FontSize', 12)
ylabel('Mean b / a', 'FontSize', 12)
title('Asymmetry Ratio vs Position', 'FontSize', 13)
grid on; box on

subplot(1,3,3)
plot(bin_mids, mean_l_bin, 'o-', ...
     'Color', [0.1 0.6 0.1], 'LineWidth', 2, ...
     'MarkerSize', 8, 'MarkerFaceColor', [0.1 0.6 0.1])
xlabel('x_m (downstream position)', 'FontSize', 12)
ylabel('Mean structure size l', 'FontSize', 12)
title('Structure Size Growth vs Position', 'FontSize', 13)
grid on; box on

saveas(gcf, 'figure4_trends.png')
fprintf('Figure 4 saved: figure4_trends.png\n')

%% =========================================================
%  PART F — RANDOM FOREST FEATURE IMPORTANCE
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' PART F: Random Forest Feature Importance\n')
fprintf('=================================================\n')

labelled = readtable('labelled.csv');

fprintf('Label distribution in training data:\n')
tabulate(labelled.label)

feature_names = {'l', 'l_over_xminusx0', 'b_over_a', ...
                 'detJ_upstream',   'traceJ_upstream',   'discJ_upstream', ...
                 'detJ_core',       'traceJ_core',       'discJ_core', ...
                 'detJ_downstream', 'traceJ_downstream', 'discJ_downstream'};

X = [labelled.l, labelled.l_over_xminusx0, labelled.b_over_a, ...
     labelled.detJ_upstream,   labelled.traceJ_upstream,   labelled.discJ_upstream, ...
     labelled.detJ_core,       labelled.traceJ_core,       labelled.discJ_core, ...
     labelled.detJ_downstream, labelled.traceJ_downstream, labelled.discJ_downstream];

Y = categorical(labelled.label);

rng(42)
fprintf('\nTraining Random Forest...\n')
mdl = TreeBagger(100, X, Y, ...
                 'Method',                 'classification', ...
                 'OOBPredictorImportance', 'on', ...
                 'PredictorNames',         feature_names);

oob_error = oobError(mdl);
fprintf('OOB Classification Error : %.4f\n', oob_error(end))
fprintf('OOB Accuracy             : %.2f%%\n', (1 - oob_error(end)) * 100)

importance = mdl.OOBPermutedPredictorDeltaError;
[sorted_imp, idx] = sort(importance, 'ascend');
sorted_names = feature_names(idx);

figure('Position', [100 100 750 520]);
barh(1:12, sorted_imp, 'FaceColor', [0.2 0.4 0.8])
yticks(1:12)
yticklabels(sorted_names)
xlabel('Feature Importance (OOB Permutation Error)', 'FontSize', 12)
title('Feature Importance — What the Model Learned', 'FontSize', 14)
grid on; box on
xlim([0 max(sorted_imp) * 1.2])
saveas(gcf, 'figure5_feature_importance.png')
fprintf('Figure 5 saved: figure5_feature_importance.png\n')

%% =========================================================
%  PART G — MANUAL ML CLEANING
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' PART G: ML Cleaning with Manual Labels\n')
fprintf('=================================================\n')

manual = readtable('manual_labels.csv');
manual = manual(manual.label ~= 2, :);

fprintf('Manually labelled rows : %d\n', height(manual))
fprintf('Label distribution:\n')
tabulate(manual.label)

% Train using 3 features only (manual_labels.csv has no Jacobian columns)
feature_names_manual = {'l', 'l_over_xminusx0', 'b_over_a'};

X_manual = [manual.l, manual.l_over_xminusx0, manual.b_over_a];
Y_manual = categorical(manual.label);

rng(42)
fprintf('\nTraining Random Forest on manual labels...\n')
mdl_manual = TreeBagger(100, X_manual, Y_manual, ...
    'Method',                 'classification', ...
    'OOBPredictorImportance', 'on', ...
    'PredictorNames',         feature_names_manual);

oob_manual = oobError(mdl_manual);
fprintf('OOB Accuracy (manual) : %.2f%%\n', (1 - oob_manual(end)) * 100)

% Apply to full raw dataset using same 3 features
X_full = [raw.l, raw.l_over_xminusx0, raw.b_over_a];

predictions   = predict(mdl_manual, X_full);
predictions   = str2double(cellstr(predictions));
raw_ml_manual = raw(predictions == 1, :);
l_norm_manual = raw_ml_manual.l_over_xminusx0;
b_a_manual    = raw_ml_manual.b_over_a;

fprintf('After manual ML cleaning : %d rows\n', height(raw_ml_manual))
fprintf('Removed                  : %d rows\n', height(raw) - height(raw_ml_manual))
fprintf('Mean l/x                 : %.4f\n', mean(l_norm_manual))
fprintf('Mean b/a                 : %.4f\n', mean(b_a_manual))

writetable(raw_ml_manual, 'manual_ml_clean_triplets.csv')
fprintf('Saved manual_ml_clean_triplets.csv\n')

% --- Figure 6: 4-way comparison ---
figure('Position', [100 100 1800 500]);

subplot(1,4,1)
histogram(l_norm_raw, 80, 'Normalization', 'pdf', ...
    'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(l_norm_raw), 'r-', 'LineWidth', 2, ...
    'Label', sprintf('Mean=%.4f', mean(l_norm_raw)), ...
    'LabelVerticalAlignment', 'bottom')
xlabel('l / (x_m - x_0)', 'FontSize', 11)
ylabel('Probability Density', 'FontSize', 11)
title(sprintf('Original\n(n=%d)', height(raw)), 'FontSize', 12)
grid on; box on

subplot(1,4,2)
histogram(l_norm_stat, 80, 'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(l_norm_stat), 'r-', 'LineWidth', 2, ...
    'Label', sprintf('Mean=%.4f', mean(l_norm_stat)), ...
    'LabelVerticalAlignment', 'bottom')
xlabel('l / (x_m - x_0)', 'FontSize', 11)
ylabel('Probability Density', 'FontSize', 11)
title(sprintf('Stat Cleaned\n(n=%d)', height(raw_stat_clean)), 'FontSize', 12)
grid on; box on

subplot(1,4,3)
histogram(l_norm_ml, 80, 'Normalization', 'pdf', ...
    'FaceColor', [0.1 0.7 0.3], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(l_norm_ml), 'r-', 'LineWidth', 2, ...
    'Label', sprintf('Mean=%.4f', mean(l_norm_ml)), ...
    'LabelVerticalAlignment', 'bottom')
xlabel('l / (x_m - x_0)', 'FontSize', 11)
ylabel('Probability Density', 'FontSize', 11)
title(sprintf('ML Cleaned (Rule Labels)\n(n=%d)', height(ml_clean)), 'FontSize', 12)
grid on; box on

subplot(1,4,4)
histogram(l_norm_manual, 80, 'Normalization', 'pdf', ...
    'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'black', 'LineWidth', 0.3)
hold on
xline(mean(l_norm_manual), 'r-', 'LineWidth', 2, ...
    'Label', sprintf('Mean=%.4f', mean(l_norm_manual)), ...
    'LabelVerticalAlignment', 'bottom')
xlabel('l / (x_m - x_0)', 'FontSize', 11)
ylabel('Probability Density', 'FontSize', 11)
title(sprintf('ML Cleaned (Manual Labels)\n(n=%d)', height(raw_ml_manual)), 'FontSize', 12)
grid on; box on

saveas(gcf, 'figure6_fourway_comparison.png')
fprintf('Figure 6 saved: figure6_fourway_comparison.png\n')

%% =========================================================
%  PART H — FINAL VALIDATION SUMMARY
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' PART H: Final Validation Summary\n')
fprintf('=================================================\n')

k = mean(l_norm_ml) / lambda;

fprintf('Lambda                    : %.3f\n', lambda)
fprintf('x_0                       : %.3f\n', x0)
fprintf('mean l/x raw              : %.4f\n', mean(l_norm_raw))
fprintf('mean l/x stat cleaned     : %.4f\n', mean(l_norm_stat))
fprintf('mean l/x ML cleaned       : %.4f\n', mean(l_norm_ml))
fprintf('k = mean(ML) / lambda     : %.4f\n', k)
fprintf('Theoretical k range       : 0.46 to 0.60\n')

if k >= 0.46 && k <= 0.60
    fprintf('\n*** VALIDATION PASSED *** k = %.4f\n', k)
else
    fprintf('\nWARNING: k = %.4f is outside expected range\n', k)
end

fprintf('\n--- Cleaning Comparison ---\n')
fprintf('Original rows             : %d\n', height(raw))
fprintf('After stat cleaning       : %d (removed %d)\n', ...
        height(raw_stat_clean), height(raw) - height(raw_stat_clean))
fprintf('After ML cleaning         : %d (removed %d)\n', ...
        height(ml_clean), height(raw) - height(ml_clean))

fprintf('\n--- Key Statistics (ML Cleaned) ---\n')
fprintf('Mean b/a                  : %.4f\n', mean(b_a_ml))
fprintf('Std  b/a                  : %.4f\n', std(b_a_ml))
fprintf('Mean l                    : %.6f\n', mean(ml_clean.l))

fprintf('\nAll figures saved. Pipeline complete.\n')
fprintf('=================================================\n')

