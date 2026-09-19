clear; clc; close all;

%% --- Load labelled training data ---
labelled = readtable('labelled.csv');
labelled = labelled(labelled.label ~= 2, :);

fprintf('Training data: %d rows\n', height(labelled))
fprintf('Label distribution:\n')
tabulate(labelled.label)

%% --- Features and labels ---
feature_names = {'l', 'l_over_xminusx0', 'b_over_a', ...
                 'detJ_upstream',   'traceJ_upstream',   'discJ_upstream', ...
                 'detJ_core',       'traceJ_core',       'discJ_core', ...
                 'detJ_downstream', 'traceJ_downstream', 'discJ_downstream'};

X = [labelled.l, labelled.l_over_xminusx0, labelled.b_over_a, ...
     labelled.detJ_upstream,   labelled.traceJ_upstream,   labelled.discJ_upstream, ...
     labelled.detJ_core,       labelled.traceJ_core,       labelled.discJ_core, ...
     labelled.detJ_downstream, labelled.traceJ_downstream, labelled.discJ_downstream];

Y = labelled.label;
Y_cat = categorical(Y);

%% --- Cross validation setup ---
rng(42)
cv = cvpartition(Y, 'KFold', 5);

%% =========================================================
%  MODEL 1 — Random Forest
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' Model 1: Random Forest\n')
fprintf('=================================================\n')

rng(42)
mdl_rf = TreeBagger(100, X, Y_cat, ...
                    'Method', 'classification', ...
                    'OOBPredictorImportance', 'on', ...
                    'PredictorNames', feature_names);

oob_rf = oobError(mdl_rf);
acc_rf = (1 - oob_rf(end)) * 100;
fprintf('OOB Accuracy : %.2f%%\n', acc_rf)

pred_rf  = predict(mdl_rf, X);
pred_rf  = str2double(cellstr(pred_rf));
cm_rf    = confusionmat(Y, pred_rf);
fprintf('Confusion Matrix:\n')
disp(cm_rf)

%% =========================================================
%  MODEL 2 — Support Vector Machine (SVM)
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' Model 2: Support Vector Machine\n')
fprintf('=================================================\n')

rng(42)
mdl_svm  = fitcsvm(X, Y, ...
                   'KernelFunction', 'rbf', ...
                   'Standardize',    true, ...
                   'CrossVal',       'on', ...
                   'CVPartition',    cv);

acc_svm = (1 - kfoldLoss(mdl_svm)) * 100;
fprintf('5-Fold CV Accuracy : %.2f%%\n', acc_svm)

mdl_svm_full = fitcsvm(X, Y, 'KernelFunction', 'rbf', 'Standardize', true);
pred_svm     = predict(mdl_svm_full, X);
cm_svm       = confusionmat(Y, pred_svm);
fprintf('Confusion Matrix:\n')
disp(cm_svm)

%% =========================================================
%  MODEL 3 — k-Nearest Neighbours
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' Model 3: k-Nearest Neighbours\n')
fprintf('=================================================\n')

rng(42)
mdl_knn  = fitcknn(X, Y, ...
                   'NumNeighbors', 5, ...
                   'Standardize',  true, ...
                   'CrossVal',     'on', ...
                   'CVPartition',  cv);

acc_knn = (1 - kfoldLoss(mdl_knn)) * 100;
fprintf('5-Fold CV Accuracy : %.2f%%\n', acc_knn)

mdl_knn_full = fitcknn(X, Y, 'NumNeighbors', 5, 'Standardize', true);
pred_knn     = predict(mdl_knn_full, X);
cm_knn       = confusionmat(Y, pred_knn);
fprintf('Confusion Matrix:\n')
disp(cm_knn)

%% =========================================================
%  MODEL 4 — Decision Tree
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' Model 4: Decision Tree\n')
fprintf('=================================================\n')

rng(42)
mdl_dt  = fitctree(X, Y, ...
                   'CrossVal',    'on', ...
                   'CVPartition', cv);

acc_dt = (1 - kfoldLoss(mdl_dt)) * 100;
fprintf('5-Fold CV Accuracy : %.2f%%\n', acc_dt)

mdl_dt_full = fitctree(X, Y);
pred_dt     = predict(mdl_dt_full, X);
cm_dt       = confusionmat(Y, pred_dt);
fprintf('Confusion Matrix:\n')
disp(cm_dt)

%% =========================================================
%  COMPARISON SUMMARY
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' ML Algorithm Comparison Summary\n')
fprintf('=================================================\n')

models   = {'Random Forest', 'SVM (RBF)', 'k-NN (k=5)', 'Decision Tree'};
accuracies = [acc_rf, acc_svm, acc_knn, acc_dt];

fprintf('%-20s %10s\n', 'Model', 'Accuracy')
fprintf('%s\n', repmat('-', 1, 32))
for i = 1:4
    fprintf('%-20s %9.2f%%\n', models{i}, accuracies(i))
end

[best_acc, best_idx] = max(accuracies);
fprintf('\nBest model: %s (%.2f%%)\n', models{best_idx}, best_acc)

%% =========================================================
%  FIGURE 17 — Accuracy comparison bar chart
%% =========================================================

figure('Position', [100 100 700 500]);
colors_bar = [0.2 0.4 0.8; 0.9 0.3 0.1; 0.1 0.7 0.3; 0.8 0.6 0.1];
b = bar(accuracies, 'FaceColor', 'flat');
for i = 1:4
    b.CData(i,:) = colors_bar(i,:);
end
hold on
yline(90, 'r--', 'LineWidth', 1.5, 'Label', '90% threshold')
xticks(1:4)
xticklabels(models)
ylabel('Accuracy (%)', 'FontSize', 12)
title('ML Algorithm Comparison — Classification Accuracy', 'FontSize', 14)
ylim([80 105])
grid on; box on

% Add accuracy labels on bars
for i = 1:4
    text(i, accuracies(i) + 0.5, sprintf('%.2f%%', accuracies(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold')
end

saveas(gcf, 'figure17_ml_comparison.png')
fprintf('\nFigure 17 saved: figure17_ml_comparison.png\n')

%% =========================================================
%  FIGURE 18 — Apply best model to isothermal dataset
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' Applying best model to isothermal dataset\n')
fprintf('=================================================\n')

% Load raw isothermal data
function raw = load_csv(fname)
    raw = readtable(fname, 'VariableNamingRule', 'preserve');
    if startsWith(raw.Properties.VariableNames{1}, 'Var')
        fid = fopen(fname, 'rb');
        content = fread(fid, Inf, '*char')';
        fclose(fid);
        nl = find(content == newline, 1);
        header = content(1:nl-1);
        if ~isempty(header) && header(end) == char(13)
            header = header(1:end-1);
        end
        tokens = strsplit(header, ',');
        merged = {};
        for i = 1:length(tokens)
            t = char(strtrim(tokens{i}));
            if isempty(t)
            elseif contains(lower(t), 'accept')
                merged{end+1} = 'label';
            elseif contains(lower(t), 'reject') || contains(lower(t), 'unlabelled')
            else
                merged{end+1} = t;
            end
        end
        if width(raw) == length(merged)
            raw.Properties.VariableNames = merged;
        end
    end
end

raw = load_csv('isothermal_all_triplets.csv');

X_full = [raw.l, raw.l_over_xminusx0, raw.b_over_a, ...
          raw.detJ_upstream,   raw.traceJ_upstream,   raw.discJ_upstream, ...
          raw.detJ_core,       raw.traceJ_core,       raw.discJ_core, ...
          raw.detJ_downstream, raw.traceJ_downstream, raw.discJ_downstream];

lambda = 0.428;

% Apply all models and compare k values
pred_rf_full  = str2double(cellstr(predict(mdl_rf, X_full)));
pred_svm_full = predict(mdl_svm_full, X_full);
pred_knn_full = predict(mdl_knn_full, X_full);
pred_dt_full  = predict(mdl_dt_full,  X_full);

k_rf  = mean(raw.l_over_xminusx0(pred_rf_full  == 1)) / lambda;
k_svm = mean(raw.l_over_xminusx0(pred_svm_full == 1)) / lambda;
k_knn = mean(raw.l_over_xminusx0(pred_knn_full == 1)) / lambda;
k_dt  = mean(raw.l_over_xminusx0(pred_dt_full  == 1)) / lambda;

n_rf  = sum(pred_rf_full  == 1);
n_svm = sum(pred_svm_full == 1);
n_knn = sum(pred_knn_full == 1);
n_dt  = sum(pred_dt_full  == 1);

fprintf('\n%-20s %10s %10s\n', 'Model', 'Rows kept', 'k value')
fprintf('%s\n', repmat('-', 1, 42))
fprintf('%-20s %10d %10.4f\n', 'Random Forest',  n_rf,  k_rf)
fprintf('%-20s %10d %10.4f\n', 'SVM',            n_svm, k_svm)
fprintf('%-20s %10d %10.4f\n', 'k-NN',           n_knn, k_knn)
fprintf('%-20s %10d %10.4f\n', 'Decision Tree',  n_dt,  k_dt)

% Figure 18 — k values per model
k_vals_models = [k_rf, k_svm, k_knn, k_dt];
figure('Position', [100 100 700 500]);
b2 = bar(k_vals_models, 'FaceColor', 'flat');
for i = 1:4
    b2.CData(i,:) = colors_bar(i,:);
end
hold on
yline(0.46, 'r--', 'LineWidth', 2, 'Label', 'k=0.46')
yline(0.60, 'r--', 'LineWidth', 2, 'Label', 'k=0.60')
xticks(1:4)
xticklabels(models)
ylabel('k value', 'FontSize', 12)
title('k Value by ML Model — Isothermal Case', 'FontSize', 14)
ylim([0.3 0.7])
grid on; box on
for i = 1:4
    text(i, k_vals_models(i) + 0.005, sprintf('%.4f', k_vals_models(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold')
end
saveas(gcf, 'figure18_k_by_model.png')
fprintf('\nFigure 18 saved: figure18_k_by_model.png\n')

fprintf('\nML comparison complete.\n')
fprintf('=================================================\n')