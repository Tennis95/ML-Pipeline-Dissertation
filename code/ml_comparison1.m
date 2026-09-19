%% ML Comparison - All 8 Cases
% Compares Random Forest, SVM, k-NN, Decision Tree on all reactive datasets

cases = {
    'isothermal_all_triplets.csv',   'Isothermal', 0;
    'hs1f-hs1h-all_triplets.csv', 'phi=1.0',    1.0;
    'hs1f-ls8h_all_triplets.csv', 'phi=8.0',    8.0;
    'hs8h-ls1f_all_triplets.csv', 'phi=0.125',  0.125;
    'hs1f-ls2h_all_triplets.csv', 'phi=2.0',    2.0;
    'hs1f-ls4h_all_triplets.csv', 'phi=4.0',    4.0;
    'hs2h-ls1f_all_triplets.csv', 'phi=0.5',    0.5;
    'hs4h-ls1f_all_triplets.csv', 'phi=0.25',   0.25;
};

nCases = size(cases, 1);
results = struct();

lambda = 0.428;
theory_k_low  = 0.46;
theory_k_high = 0.60;

acc_RF  = zeros(nCases,1);
acc_SVM = zeros(nCases,1);
acc_KNN = zeros(nCases,1);
acc_DT  = zeros(nCases,1);
k_RF    = zeros(nCases,1);
k_SVM   = zeros(nCases,1);
k_KNN   = zeros(nCases,1);
k_DT    = zeros(nCases,1);
k_stat  = zeros(nCases,1);
n_stat  = zeros(nCases,1);
n_RF    = zeros(nCases,1);
n_SVM   = zeros(nCases,1);
n_KNN   = zeros(nCases,1);
n_DT    = zeros(nCases,1);

%% Helper: load CSV with broken header fix
function T = load_csv(fname)
    T = readtable(fname);
    if startsWith(T.Properties.VariableNames{1}, 'Var')
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
                % skip
            else
                merged{end+1} = t;
            end
        end
        if width(T) == length(merged)
            T.Properties.VariableNames = merged;
        end
    end
end

%% Helper: compute k from accepted rows
function k = compute_k(T, lambda)
    vals = T.l_over_xminusx0;
    k = mean(vals) / lambda;
end

%% Main loop
for ci = 1:nCases
    fname = cases{ci,1};
    label = cases{ci,2};
    fprintf('\n=== %s ===\n', label);

    raw = load_csv(fname);

    % Statistical cleaning
    if ~ismember('l', raw.Properties.VariableNames)
        fprintf('  Skipping %s - missing columns\n', fname);
        continue;
    end

    raw.l = raw.x_downstream - raw.x_upstream;
    raw.l_over_xminusx0 = raw.l ./ (0.5*(raw.x_upstream + raw.x_downstream));
    if ~ismember('b_over_a', raw.Properties.VariableNames)
        raw.b_over_a = abs(raw.x_core - raw.x_upstream) ./ abs(raw.x_downstream - raw.x_core);
    end

    mu = mean(raw.l_over_xminusx0);
    sg = std(raw.l_over_xminusx0);
    stat_mask = raw.l_over_xminusx0 > (mu - sg) & raw.l_over_xminusx0 < (mu + sg);
    clean = raw(stat_mask, :);
    n_stat(ci) = height(clean);
    k_stat(ci) = compute_k(clean, lambda);
    fprintf('  Statistical: n=%d, k=%.4f\n', n_stat(ci), k_stat(ci));

    % Auto-label using physical rules
    labels = zeros(height(raw), 1);
    for r = 1:height(raw)
        if raw.l_over_xminusx0(r) >= 0.1 && raw.b_over_a(r) >= 0.3 && raw.b_over_a(r) <= 5.0
            labels(r) = 1;
        end
    end

    if sum(labels==1) < 10 || sum(labels==0) < 10
        fprintf('  Not enough labelled data, skipping ML\n');
        acc_RF(ci) = NaN; acc_SVM(ci) = NaN;
        acc_KNN(ci) = NaN; acc_DT(ci) = NaN;
        k_RF(ci) = NaN; k_SVM(ci) = NaN;
        k_KNN(ci) = NaN; k_DT(ci) = NaN;
        continue;
    end

    features = [raw.l, raw.l_over_xminusx0, raw.b_over_a];
    X = features;
    y = labels;

    cv = cvpartition(y, 'KFold', 5, 'Stratify', true);

    %% Random Forest
    rf_preds = zeros(size(y));
    for f = 1:cv.NumTestSets
        Xtr = X(cv.training(f),:); ytr = y(cv.training(f));
        Xte = X(cv.test(f),:);
        mdl = TreeBagger(50, Xtr, ytr, 'OOBPrediction','off', 'Method','classification');
        p = str2double(predict(mdl, Xte));
        rf_preds(cv.test(f)) = p;
    end
    acc_RF(ci) = mean(rf_preds == y) * 100;
    acc_idx = rf_preds == 1;
    if sum(acc_idx) > 0
        sub = raw(acc_idx,:);
        k_RF(ci) = compute_k(sub, lambda);
        n_RF(ci) = sum(acc_idx);
    else
        k_RF(ci) = NaN; n_RF(ci) = 0;
    end
    fprintf('  RF:  acc=%.1f%%, k=%.4f, n=%d\n', acc_RF(ci), k_RF(ci), n_RF(ci));

    %% SVM
    svm_preds = zeros(size(y));
    for f = 1:cv.NumTestSets
        Xtr = X(cv.training(f),:); ytr = y(cv.training(f));
        Xte = X(cv.test(f),:);
        mdl = fitcsvm(Xtr, ytr, 'KernelFunction','rbf', 'Standardize',true);
        svm_preds(cv.test(f)) = predict(mdl, Xte);
    end
    acc_SVM(ci) = mean(svm_preds == y) * 100;
    acc_idx = svm_preds == 1;
    if sum(acc_idx) > 0
        sub = raw(acc_idx,:);
        k_SVM(ci) = compute_k(sub, lambda);
        n_SVM(ci) = sum(acc_idx);
    else
        k_SVM(ci) = NaN; n_SVM(ci) = 0;
    end
    fprintf('  SVM: acc=%.1f%%, k=%.4f, n=%d\n', acc_SVM(ci), k_SVM(ci), n_SVM(ci));

    %% k-NN
    knn_preds = zeros(size(y));
    for f = 1:cv.NumTestSets
        Xtr = X(cv.training(f),:); ytr = y(cv.training(f));
        Xte = X(cv.test(f),:);
        mdl = fitcknn(Xtr, ytr, 'NumNeighbors',5, 'Standardize',true);
        knn_preds(cv.test(f)) = predict(mdl, Xte);
    end
    acc_KNN(ci) = mean(knn_preds == y) * 100;
    acc_idx = knn_preds == 1;
    if sum(acc_idx) > 0
        sub = raw(acc_idx,:);
        k_KNN(ci) = compute_k(sub, lambda);
        n_KNN(ci) = sum(acc_idx);
    else
        k_KNN(ci) = NaN; n_KNN(ci) = 0;
    end
    fprintf('  KNN: acc=%.1f%%, k=%.4f, n=%d\n', acc_KNN(ci), k_KNN(ci), n_KNN(ci));

    %% Decision Tree
    dt_preds = zeros(size(y));
    for f = 1:cv.NumTestSets
        Xtr = X(cv.training(f),:); ytr = y(cv.training(f));
        Xte = X(cv.test(f),:);
        mdl = fitctree(Xtr, ytr);
        dt_preds(cv.test(f)) = predict(mdl, Xte);
    end
    acc_DT(ci) = mean(dt_preds == y) * 100;
    acc_idx = dt_preds == 1;
    if sum(acc_idx) > 0
        sub = raw(acc_idx,:);
        k_DT(ci) = compute_k(sub, lambda);
        n_DT(ci) = sum(acc_idx);
    else
        k_DT(ci) = NaN; n_DT(ci) = 0;
    end
    fprintf('  DT:  acc=%.1f%%, k=%.4f, n=%d\n', acc_DT(ci), k_DT(ci), n_DT(ci));
end

%% Figure 1: k values per case, all methods
fig1 = figure('Name','k Values - All Cases');
case_labels = cases(:,2);
phi_vals = cell2mat(cases(:,3));
[~, sort_idx] = sort(phi_vals);

x_pos = 1:nCases;
hold on;
plot(x_pos, k_stat(sort_idx), 'k-o', 'LineWidth',1.5, 'DisplayName','Statistical');
plot(x_pos, k_RF(sort_idx),   'b-s', 'LineWidth',1.5, 'DisplayName','Random Forest');
plot(x_pos, k_SVM(sort_idx),  'r-^', 'LineWidth',1.5, 'DisplayName','SVM');
plot(x_pos, k_KNN(sort_idx),  'g-d', 'LineWidth',1.5, 'DisplayName','k-NN');
plot(x_pos, k_DT(sort_idx),   'm-p', 'LineWidth',1.5, 'DisplayName','Decision Tree');
yline(theory_k_low,  '--k', 'k=0.46', 'LabelHorizontalAlignment','left');
yline(theory_k_high, '--k', 'k=0.60', 'LabelHorizontalAlignment','left');
set(gca, 'XTick', x_pos, 'XTickLabel', case_labels(sort_idx), 'XTickLabelRotation', 30);
xlabel('Case'); ylabel('k = mean(l/(x_m-x_0)) / \lambda');
title('k Values Across All Cases and ML Methods');
legend('Location','best');
grid on;
ylim([0.3 0.8]);
saveas(fig1, 'fig19_k_all_cases.png');

%% Figure 2: Accuracy comparison
fig2 = figure('Name','ML Accuracy - All Cases');
bar_data = [acc_RF(sort_idx), acc_SVM(sort_idx), acc_KNN(sort_idx), acc_DT(sort_idx)];
bar(x_pos, bar_data, 'grouped');
set(gca, 'XTick', x_pos, 'XTickLabel', case_labels(sort_idx), 'XTickLabelRotation', 30);
xlabel('Case'); ylabel('Cross-Validation Accuracy (%)');
title('ML Algorithm Accuracy Across All Cases');
legend({'Random Forest','SVM','k-NN','Decision Tree'}, 'Location','best');
grid on;
ylim([0 100]);
saveas(fig2, 'fig20_accuracy_all_cases.png');

%% Figure 3: k vs phi for each method
phi_sorted = phi_vals(sort_idx);
phi_sorted_nonzero = phi_sorted;
phi_sorted_nonzero(phi_sorted == 0) = NaN;

fig3 = figure('Name','k vs Equivalence Ratio');
hold on;
plot(phi_sorted, k_RF(sort_idx),   'b-s', 'LineWidth',1.5, 'DisplayName','Random Forest');
plot(phi_sorted, k_SVM(sort_idx),  'r-^', 'LineWidth',1.5, 'DisplayName','SVM');
plot(phi_sorted, k_KNN(sort_idx),  'g-d', 'LineWidth',1.5, 'DisplayName','k-NN');
plot(phi_sorted, k_DT(sort_idx),   'm-p', 'LineWidth',1.5, 'DisplayName','Decision Tree');
plot(phi_sorted, k_stat(sort_idx), 'k-o', 'LineWidth',1.5, 'DisplayName','Statistical');
yline(theory_k_low,  '--k');
yline(theory_k_high, '--k');
xlabel('Equivalence Ratio \phi'); ylabel('k');
title('k vs Equivalence Ratio \phi — All ML Methods');
legend('Location','best');
grid on;
ylim([0.3 0.8]);
saveas(fig3, 'fig21_k_vs_phi.png');

fprintf('\n=== DONE: figures fig19, fig20, fig21 saved ===\n');