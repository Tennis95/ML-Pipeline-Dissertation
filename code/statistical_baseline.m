clear; clc; close all;

%% --- Parameters ---
lambda = 0.428;
n_std  = 1.0;
x_min  = 0.0;
x_max  = 0.16;

%% --- Define all cases ---
cases = {
    'isothermal_all_triplets.csv',   'Isothermal',  0;
    'hs1f-hs1h-all_triplets.csv', 'phi=1.0',     1.0;
    'hs1f-ls8h_all_triplets.csv', 'phi=8.0',     8.0;
    'hs8h-ls1f_all_triplets.csv', 'phi=0.125',   0.125;
    'hs1f-ls2h_all_triplets.csv', 'phi=2.0',     2.0;
    'hs1f-ls4h_all_triplets.csv', 'phi=4.0',     4.0;
    'hs2h-ls1f_all_triplets.csv', 'phi=0.5',     0.5;
    'hs4h-ls1f_all_triplets.csv', 'phi=0.25',    0.25;
};

n_cases = size(cases, 1);
case_labels = cases(:, 2);

% Storage — full domain
mean_l_norm_full = zeros(1, n_cases);
mean_b_a_full    = zeros(1, n_cases);
k_full           = zeros(1, n_cases);
n_orig           = zeros(1, n_cases);
n_clean_full     = zeros(1, n_cases);

% Storage — spatial filter 0 < x < 0.16m
mean_l_norm_filt = zeros(1, n_cases);
mean_b_a_filt    = zeros(1, n_cases);
k_filt           = zeros(1, n_cases);
n_clean_filt     = zeros(1, n_cases);

%% =========================================================
%  HELPER — Load CSV with broken header fix
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
            fprintf('WARNING: Column mismatch for %s\n', fname)
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

    %% Full domain statistical cleaning
    l_norm = raw.l_over_xminusx0;
    b_a    = raw.b_over_a;
    x_m    = 0.5 * (raw.x_upstream + raw.x_downstream);

    mu    = mean(l_norm);
    sig   = std(l_norm);
    mask_stat = l_norm >= mu - n_std*sig & l_norm <= mu + n_std*sig;
    clean_full = raw(mask_stat, :);

    n_orig(c)          = height(raw);
    n_clean_full(c)    = height(clean_full);
    mean_l_norm_full(c) = mean(clean_full.l_over_xminusx0);
    mean_b_a_full(c)   = mean(clean_full.b_over_a);
    k_full(c)          = mean_l_norm_full(c) / lambda;

    fprintf('--- Full domain ---\n')
    fprintf('Original : %d  |  Cleaned : %d  |  Removed : %d\n', ...
            n_orig(c), n_clean_full(c), n_orig(c)-n_clean_full(c))
    fprintf('Mean l/x : %.4f  |  Mean b/a : %.4f  |  k : %.4f\n', ...
            mean_l_norm_full(c), mean_b_a_full(c), k_full(c))

    writetable(clean_full, strrep(fname, '_all_triplets.csv', '_full_clean.csv'))

    %% Spatial filter: 0 < x_m < 0.16m
    x_m_clean = 0.5 * (clean_full.x_upstream + clean_full.x_downstream);
    mask_spatial = x_m_clean > x_min & x_m_clean < x_max;
    clean_filt = clean_full(mask_spatial, :);

    n_clean_filt(c)    = height(clean_filt);
    mean_l_norm_filt(c) = mean(clean_filt.l_over_xminusx0);
    mean_b_a_filt(c)   = mean(clean_filt.b_over_a);
    k_filt(c)          = mean_l_norm_filt(c) / lambda;

    fprintf('--- Spatial filter (0 < x < 0.16m) ---\n')
    fprintf('Rows in region : %d\n', n_clean_filt(c))
    fprintf('Mean l/x : %.4f  |  Mean b/a : %.4f  |  k : %.4f\n', ...
            mean_l_norm_filt(c), mean_b_a_filt(c), k_filt(c))

    writetable(clean_filt, strrep(fname, '_all_triplets.csv', '_filt_clean.csv'))

end

%% =========================================================
%  SUMMARY TABLE
%% =========================================================

fprintf('\n=================================================\n')
fprintf(' Full Domain Summary\n')
fprintf('=================================================\n')
fprintf('%-15s %8s %8s %8s %8s %8s\n', ...
        'Case','N_orig','N_clean','mean_l/x','mean_b/a','k')
fprintf('%s\n', repmat('-',1,60))
for c = 1:n_cases
    fprintf('%-15s %8d %8d %8.4f %8.4f %8.4f\n', ...
            case_labels{c}, n_orig(c), n_clean_full(c), ...
            mean_l_norm_full(c), mean_b_a_full(c), k_full(c))
end

fprintf('\n=================================================\n')
fprintf(' Spatial Filter Summary (0 < x < 0.16m)\n')
fprintf('=================================================\n')
fprintf('%-15s %8s %8s %8s %8s\n', ...
        'Case','N_filt','mean_l/x','mean_b/a','k')
fprintf('%s\n', repmat('-',1,50))
for c = 1:n_cases
    fprintf('%-15s %8d %8.4f %8.4f %8.4f\n', ...
            case_labels{c}, n_clean_filt(c), ...
            mean_l_norm_filt(c), mean_b_a_filt(c), k_filt(c))
end

%% =========================================================
%  FIGURES — FULL DOMAIN
%% =========================================================

colors = {[0.4 0.4 0.4],[0.2 0.6 0.9],[0.1 0.7 0.3],[0.9 0.3 0.1],...
          [0.6 0.2 0.8],[0.9 0.7 0.1],[0.1 0.5 0.5],[0.8 0.4 0.2]};

% Figure 11 — k values full domain
figure('Position',[100 100 1000 500]);
bar(1:n_cases, k_full, 'FaceColor', [0.2 0.4 0.8])
hold on
yline(0.46,'r--','LineWidth',2,'Label','k=0.46')
yline(0.60,'r--','LineWidth',2,'Label','k=0.60')
xticks(1:n_cases); xticklabels(case_labels)
ylabel('k = mean(l/x) / \lambda','FontSize',12)
title('k Values — All Cases (Full Domain)','FontSize',14)
grid on; box on
saveas(gcf,'figure11_k_all_cases_full.png')
fprintf('\nFigure 11 saved\n')

% Figure 12 — b/a full domain
figure('Position',[100 100 1000 500]);
bar(1:n_cases, mean_b_a_full, 'FaceColor', [0.9 0.5 0.1])
hold on
yline(1.0,'g--','LineWidth',2,'Label','b/a=1')
xticks(1:n_cases); xticklabels(case_labels)
ylabel('Mean b/a','FontSize',12)
title('Mean b/a — All Cases (Full Domain)','FontSize',14)
grid on; box on
saveas(gcf,'figure12_ba_all_cases_full.png')
fprintf('Figure 12 saved\n')

%% =========================================================
%  FIGURES — SPATIAL FILTER 0 < x < 0.16m
%% =========================================================

% Figure 13 — k values filtered
figure('Position',[100 100 1000 500]);
bar(1:n_cases, k_filt, 'FaceColor', [0.2 0.7 0.4])
hold on
yline(0.46,'r--','LineWidth',2,'Label','k=0.46')
yline(0.60,'r--','LineWidth',2,'Label','k=0.60')
xticks(1:n_cases); xticklabels(case_labels)
ylabel('k = mean(l/x) / \lambda','FontSize',12)
title('k Values — All Cases (0 < x < 0.16m)','FontSize',14)
grid on; box on
saveas(gcf,'figure13_k_all_cases_filtered.png')
fprintf('Figure 13 saved\n')

% Figure 14 — b/a filtered
figure('Position',[100 100 1000 500]);
bar(1:n_cases, mean_b_a_filt, 'FaceColor', [0.8 0.3 0.3])
hold on
yline(1.0,'g--','LineWidth',2,'Label','b/a=1')
xticks(1:n_cases); xticklabels(case_labels)
ylabel('Mean b/a','FontSize',12)
title('Mean b/a — All Cases (0 < x < 0.16m)','FontSize',14)
grid on; box on
saveas(gcf,'figure14_ba_all_cases_filtered.png')
fprintf('Figure 14 saved\n')

% Figure 15 — PDF l/x all cases filtered
figure('Position',[100 100 2000 400]);
for c = 1:n_cases
    raw2   = load_csv(cases{c,1});
    l2     = raw2.l_over_xminusx0;
    mu2    = mean(l2); sig2 = std(l2);
    cf     = raw2(l2 >= mu2-n_std*sig2 & l2 <= mu2+n_std*sig2, :);
    xm2    = 0.5*(cf.x_upstream + cf.x_downstream);
    cf     = cf(xm2 > x_min & xm2 < x_max, :);
    subplot(2,4,c)
    histogram(cf.l_over_xminusx0, 50, 'Normalization','pdf',...
              'FaceColor',colors{c},'EdgeColor','black','LineWidth',0.3)
    hold on
    xline(mean(cf.l_over_xminusx0),'r-','LineWidth',2)
    xlabel('l/(x_m-x_0)','FontSize',10)
    ylabel('PDF','FontSize',10)
    title(sprintf('%s\nn=%d, k=%.3f', case_labels{c}, height(cf), k_filt(c)),'FontSize',11)
    grid on; box on
end
sgtitle('PDF of l/(x-x_0) — Filtered Region (0 < x < 0.16m)','FontSize',13)
saveas(gcf,'figure15_pdf_filtered_all.png')
fprintf('Figure 15 saved\n')

% Figure 16 — PDF b/a all cases filtered
figure('Position',[100 100 2000 400]);
for c = 1:n_cases
    raw2   = load_csv(cases{c,1});
    l2     = raw2.l_over_xminusx0;
    mu2    = mean(l2); sig2 = std(l2);
    cf     = raw2(l2 >= mu2-n_std*sig2 & l2 <= mu2+n_std*sig2, :);
    xm2    = 0.5*(cf.x_upstream + cf.x_downstream);
    cf     = cf(xm2 > x_min & xm2 < x_max, :);
    subplot(2,4,c)
    histogram(cf.b_over_a, 50, 'Normalization','pdf',...
              'FaceColor',colors{c},'EdgeColor','black','LineWidth',0.3)
    hold on
    xline(1.0,'g--','LineWidth',2)
    xline(mean(cf.b_over_a),'r-','LineWidth',2)
    xlabel('b/a','FontSize',10)
    ylabel('PDF','FontSize',10)
    title(sprintf('%s\nn=%d, mean=%.3f', case_labels{c}, height(cf), mean(cf.b_over_a)),'FontSize',11)
    xlim([0 8]); grid on; box on
end
sgtitle('PDF of b/a — Filtered Region (0 < x < 0.16m)','FontSize',13)
saveas(gcf,'figure16_ba_filtered_all.png')
fprintf('Figure 16 saved\n')

fprintf('\nAll analysis complete.\n')
fprintf('=================================================\n')