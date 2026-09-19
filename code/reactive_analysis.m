clear; clc; close all;

%% --- Parameters ---
lambda = 0.428;
n_std  = 1.0;

%% --- Define all cases ---
cases = {
    'isothermal_all_triplets.csv',   'Isothermal',  0;
    'hs1f-hs1h-all_triplets.csv', 'phi=1.0',     1.0;
    'hs1f-ls8h_all_triplets.csv', 'phi=8.0',     8.0;
    'hs8h-ls1f_all_triplets.csv', 'phi=0.125',   0.125;
};

n_cases = size(cases, 1);

mean_l_norm = zeros(1, n_cases);
mean_b_a    = zeros(1, n_cases);
k_values    = zeros(1, n_cases);
n_original  = zeros(1, n_cases);
n_cleaned   = zeros(1, n_cases);
case_labels = cases(:, 2);

%% =========================================================
%  HELPER FUNCTION — Load CSV with broken header fix
%% =========================================================

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
        else
            fprintf('WARNING: Column mismatch (%d names, %d cols) for %s\n', ...
                    length(merged), width(raw), fname)
        end
    end
end

%% =========================================================
%  PROCESS EACH CASE
%% =========================================================

for c = 1:n_cases

    fname  = cases{c, 1};
    clabel = cases{c, 2};

    fprintf('\n=================================================\n')
    fprintf(' Processing: %s\n', clabel)
    fprintf('=================================================\n')

    raw = load_csv(fname);
    fprintf('Loaded: %d rows x %d columns\n', height(raw), width(raw))

    l_norm = raw.l_over_xminusx0;
    b_a    = raw.b_over_a;

    mu    = mean(l_norm);
    sig   = std(l_norm);
    lower = mu - n_std * sig;
    upper = mu + n_std * sig;

    mask  = l_norm >= lower & l_norm <= upper;
    clean = raw(mask, :);

    l_norm_clean = clean.l_over_xminusx0;
    b_a_clean    = clean.b_over_a;

    n_original(c)  = height(raw);
    n_cleaned(c)   = height(clean);
    mean_l_norm(c) = mean(l_norm_clean);
    mean_b_a(c)    = mean(b_a_clean);
    k_values(c)    = mean(l_norm_clean) / lambda;

    fprintf('Original rows  : %d\n', n_original(c))
    fprintf('After cleaning : %d\n', n_cleaned(c))
    fprintf('Removed        : %d\n', n_original(c) - n_cleaned(c))
    fprintf('Mean l/x       : %.4f\n', mean_l_norm(c))
    fprintf('Mean b/a       : %.4f\n', mean_b_a(c))
    fprintf('k              : %.4f\n', k_values(c))

    out_name = strrep(fname, '_all_triplets.csv', '_stat_clean.csv');
    if strcmp(fname, 'isothermal_all_triplets.csv')
        out_name = 'isothermal_stat_clean.csv';
    end
    writetable(clean, out_name)
    fprintf('Saved: %s\n', out_name)

end

%% =========================================================
%  SUMMARY TABLE
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' Comparison Summary\n')
fprintf('=================================================\n')
fprintf('\n%-20s %8s %8s %8s %8s %8s\n', ...
        'Case', 'N_orig', 'N_clean', 'mean_l/x', 'mean_b/a', 'k')
fprintf('%s\n', repmat('-', 1, 65))
for c = 1:n_cases
    fprintf('%-20s %8d %8d %8.4f %8.4f %8.4f\n', ...
            case_labels{c}, n_original(c), n_cleaned(c), ...
            mean_l_norm(c), mean_b_a(c), k_values(c))
end

fprintf('\nTheoretical k range: 0.46 to 0.60\n')
for c = 1:n_cases
    if k_values(c) >= 0.46 && k_values(c) <= 0.60
        fprintf('%-20s k=%.4f  PASSED\n', case_labels{c}, k_values(c))
    else
        fprintf('%-20s k=%.4f  OUTSIDE RANGE\n', case_labels{c}, k_values(c))
    end
end

%% =========================================================
%  FIGURE 7 — k values across cases
%% =========================================================

figure('Position', [100 100 800 500]);
bar(1:n_cases, k_values, 'FaceColor', [0.2 0.4 0.8])
hold on
yline(0.46, 'r--', 'LineWidth', 2, 'Label', 'k = 0.46')
yline(0.60, 'r--', 'LineWidth', 2, 'Label', 'k = 0.60')
xticks(1:n_cases)
xticklabels(case_labels)
ylabel('k = mean(l/x) / \lambda', 'FontSize', 12)
title('Structure Spacing Parameter k — All Cases', 'FontSize', 14)
grid on; box on
saveas(gcf, 'figure7_k_comparison.png')
fprintf('\nFigure 7 saved: figure7_k_comparison.png\n')

%% =========================================================
%  FIGURE 8 — mean b/a across cases
%% =========================================================

figure('Position', [100 100 800 500]);
bar(1:n_cases, mean_b_a, 'FaceColor', [0.9 0.5 0.1])
hold on
yline(1.0, 'g--', 'LineWidth', 2, 'Label', 'b/a = 1 (symmetric)')
xticks(1:n_cases)
xticklabels(case_labels)
ylabel('Mean b/a', 'FontSize', 12)
title('Structure Asymmetry b/a — All Cases', 'FontSize', 14)
grid on; box on
saveas(gcf, 'figure8_ba_comparison.png')
fprintf('Figure 8 saved: figure8_ba_comparison.png\n')

%% =========================================================
%  FIGURE 9 — PDFs side by side
%% =========================================================

colors = {[0.5 0.5 0.5], [0.2 0.6 0.9], [0.1 0.7 0.3], [0.9 0.3 0.1]};

figure('Position', [100 100 1800 400]);

for c = 1:n_cases
    fname  = cases{c, 1};
    clabel = cases{c, 2};

    raw2   = load_csv(fname);
    l2     = raw2.l_over_xminusx0;
    mu2    = mean(l2);
    sig2   = std(l2);
    clean2 = raw2(l2 >= mu2 - n_std*sig2 & l2 <= mu2 + n_std*sig2, :);

    subplot(1, n_cases, c)
    histogram(clean2.l_over_xminusx0, 60, 'Normalization', 'pdf', ...
              'FaceColor', colors{c}, 'EdgeColor', 'black', 'LineWidth', 0.3)
    hold on
    xline(mean(clean2.l_over_xminusx0), 'r-', 'LineWidth', 2)
    xlabel('l / (x_m - x_0)', 'FontSize', 11)
    ylabel('PDF', 'FontSize', 11)
    title(sprintf('%s\nn=%d, mean=%.4f', clabel, height(clean2), ...
          mean(clean2.l_over_xminusx0)), 'FontSize', 12)
    grid on; box on
end

saveas(gcf, 'figure9_pdf_all_cases.png')
fprintf('Figure 9 saved: figure9_pdf_all_cases.png\n')

%% =========================================================
%  FIGURE 10 — b/a PDFs side by side
%% =========================================================

figure('Position', [100 100 1800 400]);

for c = 1:n_cases
    fname  = cases{c, 1};
    clabel = cases{c, 2};

    raw2   = load_csv(fname);
    l2     = raw2.l_over_xminusx0;
    mu2    = mean(l2);
    sig2   = std(l2);
    clean2 = raw2(l2 >= mu2 - n_std*sig2 & l2 <= mu2 + n_std*sig2, :);

    subplot(1, n_cases, c)
    histogram(clean2.b_over_a, 60, 'Normalization', 'pdf', ...
              'FaceColor', colors{c}, 'EdgeColor', 'black', 'LineWidth', 0.3)
    hold on
    xline(1.0, 'g--', 'LineWidth', 2)
    xline(mean(clean2.b_over_a), 'r-', 'LineWidth', 2)
    xlabel('b / a', 'FontSize', 11)
    ylabel('PDF', 'FontSize', 11)
    title(sprintf('%s\nn=%d, mean=%.4f', clabel, height(clean2), ...
          mean(clean2.b_over_a)), 'FontSize', 12)
    xlim([0 8]); grid on; box on
end

saveas(gcf, 'figure10_ba_all_cases.png')
fprintf('Figure 10 saved: figure10_ba_all_cases.png\n')

fprintf('\nReactive analysis complete.\n')
fprintf('=================================================\n')