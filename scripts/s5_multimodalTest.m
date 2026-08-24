

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test with emotion regulation data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and prepare data

[emoregDat, ~, imagenames] = load_image_set('emotionreg')

atlas_obj = load_atlas('canlab2018_2mm'); % 489 regions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlation between amygdala and remaining regions

amyROI = fmri_data('..\data\masks\NIDER_AmygdalaROI.nii')
emoregDat.X = extract_roi_averages(emoregDat, amyROI).dat
emoregDat = apply_mask(emoregDat, amyROI, 'invert')

emoregReg = regress(emoregDat)
histogram(emoregReg.t.dat(:,1))

correlations = emoregReg.t.dat(:,1) ./ sqrt(emoregReg.t.dat(:,1).^2 + 33)
mean(correlations)
std(correlations)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %calculate correlation matrix between all regions

compute_corr_stats(emoregDat, atlas_obj)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test with multimodal data

kragelDat = fmri_data('..\data\kragel2018\kragel_2018_nat_neurosci_270_subjects_test_images.nii')
kragelMetaDat = readtable('..\data\kragel2018\kragel_2018_nat_neurosci_270_subjects_test_images_metadata')

kragelMetaDat_unique = unique(kragelMetaDat);

corrOut_kragel = zeros(size(kragelMetaDat_unique, 1), 2)

for i=1:size(kragelMetaDat_unique, 1)
    
    kragelDat_temp = get_wh_image(kragelDat, (kragelMetaDat.Studynumber == i));
    [mTemp sdTemp] = compute_corr_stats(kragelDat_temp, atlas_obj);
    corrOut_kragel(i,1) = mTemp;
    corrOut_kragel(i,2) = sdTemp;
    
end

write(horzcat(kragelMetaDat_unique, array2table(corrOut_kragel, 'VariableNames', {'corrMean', 'corrSD'})), '..\results\tables\kragelAnalyses.csv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resting state networks

%% =========================================================
% Study-by-edge heatmap of pairwise correlations between 7 Buckner networks
% - First row: mean across studies
% - Following rows: Emotion regulation + Kragel studies
% - Repeated subdomains labeled with suffix 1 / 2
%% =========================================================

%% Load network maps
[bucknermaps, mapnames] = load_image_set('bucknerlab');

% Keep first 7 maps
nNetworks = 7;
bucknermaps = bucknermaps.get_wh_image(1:nNetworks);
mapnames = string(mapnames(1:nNetworks));

%% Create edge labels
edgePairs = nchoosek(1:nNetworks, 2);
nEdges = size(edgePairs, 1);

edgeLabels = strings(nEdges,1);
for e = 1:nEdges
    edgeLabels(e) = mapnames(edgePairs(e,1)) + "-" + mapnames(edgePairs(e,2));
end

%% Clean map names a bit if needed
edgeLabels = replace(edgeLabels, "_", " ");

%% =========================================================
% Storage
%% =========================================================
nKragel = size(kragelMetaDat_unique, 1);
nStudiesTotal = nKragel + 1;   % +1 for emotion regulation

corrMat_all = nan(nStudiesTotal, nEdges);
studyLabels = strings(nStudiesTotal,1);
studyLabels(1) = "Emotion regulation";

%% =========================================================
% 1) Emotion regulation dataset
%% =========================================================
M = zeros(size(emoregDat.dat, 2), nNetworks);

for j = 1:nNetworks
    M(:,j) = extract_roi_averages(emoregDat, bucknermaps.get_wh_image(j)).dat;
end

for e = 1:nEdges
    x = M(:, edgePairs(e,1));
    y = M(:, edgePairs(e,2));
    
    valid = ~isnan(x) & ~isnan(y);
    if sum(valid) > 2
        corrMat_all(1,e) = corr(x(valid), y(valid), 'type', 'Pearson');
    end
end

%% =========================================================
% 2) Kragel studies
%% =========================================================
subdomainRaw = string(kragelMetaDat_unique.Subdomain);
subdomainRaw = subdomainRaw(:);

% Recode subdomain labels to cleaner publication-style names
subdomainClean = subdomainRaw;
subdomainClean = replace(subdomainClean, "_", " ");

% Optional manual recoding to match your figure style exactly
subdomainClean = replace(subdomainClean, "Inhibition", "Inhibitory Control");
subdomainClean = replace(subdomainClean, "ResponseSelect", "Response Selection");
subdomainClean = replace(subdomainClean, "Response Selectionion", "Response Selection"); % safeguard
subdomainClean = replace(subdomainClean, "WorkingMem", "Working Memory");
subdomainClean = replace(subdomainClean, "Images", "Aversive Images");
subdomainClean = replace(subdomainClean, "Social", "Social Rejection");
subdomainClean = replace(subdomainClean, "Sounds", "Aversive Sounds");
subdomainClean = replace(subdomainClean, "Mechanical", "Mechanical Pain");
subdomainClean = replace(subdomainClean, "Thermal", "Thermal Pain");
subdomainClean = replace(subdomainClean, "Visceral", "Visceral Pain");

% Add suffixes 1/2 to repeated subdomains
studySubLabels = strings(nKragel,1);
uniqueSubsInOrder = unique(subdomainClean, 'stable');

for s = 1:numel(uniqueSubsInOrder)
    idx = find(subdomainClean == uniqueSubsInOrder(s));
    for k = 1:numel(idx)
        studySubLabels(idx(k)) = uniqueSubsInOrder(s) + " " + string(k);
    end
end

for i = 1:nKragel
    
    kragelDat_temp = get_wh_image(kragelDat, (kragelMetaDat.Studynumber == i));
    
    M = zeros(size(kragelDat_temp.dat, 2), nNetworks);
    for j = 1:nNetworks
        M(:,j) = extract_roi_averages(kragelDat_temp, bucknermaps.get_wh_image(j)).dat;
    end
    
    for e = 1:nEdges
        x = M(:, edgePairs(e,1));
        y = M(:, edgePairs(e,2));
        
        valid = ~isnan(x) & ~isnan(y);
        if sum(valid) > 2
            corrMat_all(i+1,e) = corr(x(valid), y(valid), 'type', 'Pearson');
        end
    end
    
    studyLabels(i+1) = studySubLabels(i);
end

%% =========================================================
% 3) Sort columns by mean correlation across studies
%% =========================================================
meanCorr = mean(corrMat_all, 1, 'omitnan');
[meanCorr_sorted, sortIdx] = sort(meanCorr, 'descend');

corrMat_plot = corrMat_all(:, sortIdx);
edgeLabels_plot = edgeLabels(sortIdx);

%% =========================================================
% 4) Add mean row at top
%% =========================================================
heatmapMat = [meanCorr_sorted; corrMat_plot];
heatmapRowLabels = ["Mean across studies"; studyLabels];

%% =========================================================
% Export source data for heatmap to Excel
%% =========================================================

% Create table with y-axis labels as first column and x-axis labels as columns
sourceTable = array2table(heatmapMat, ...
    'VariableNames', matlab.lang.makeValidName(cellstr(edgeLabels_plot)));

sourceTable.Study = cellstr(heatmapRowLabels);

% Move Study column to the front
sourceTable = movevars(sourceTable, 'Study', 'Before', 1);

% Write to Excel
writetable(sourceTable, 'Figure3B_source_data.xlsx', ...
    'Sheet', 'Heatmap_source_data');

% Also export x- and y-axis labels explicitly on separate sheets
xAxisTable = table((1:numel(edgeLabels_plot))', cellstr(edgeLabels_plot), ...
    'VariableNames', {'Column_index', 'Network_pair'});

yAxisTable = table((1:numel(heatmapRowLabels))', cellstr(heatmapRowLabels), ...
    'VariableNames', {'Row_index', 'Study'});

writetable(xAxisTable, 'Figure3B_source_data.xlsx', ...
    'Sheet', 'X_axis_labels');

writetable(yAxisTable, 'Figure3B_source_data.xlsx', ...
    'Sheet', 'Y_axis_labels');


%% =========================================================
% 5) Plot
%% =========================================================
figure('Color', 'w', 'Position', [100 100 1400 700]);

imagesc(heatmapMat);
% Write mean correlations into the top row
for e = 1:nEdges
    text(e, 1, sprintf('%.2f', meanCorr_sorted(e)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'Color', 'k');
end
colormap(parula);
caxis([-1 1]);
colorbar;

xlabel('Network pair', 'FontWeight', 'bold');
ylabel('Study', 'FontWeight', 'bold');

xticks(1:nEdges);
xticklabels(edgeLabels_plot);
xtickangle(45);

yticks(1:size(heatmapMat,1));
yticklabels(heatmapRowLabels);

title('Pairwise correlations between default mode networks across studies');

% Separator below mean row
hold on;
yline(1.5, 'k-', 'LineWidth', 1.2);

% Improve appearance
set(gca, 'TickLength', [0 0], 'FontSize', 11);
box off;


%% =========================================================
% 5) Plot Vector Nature
%% =========================================================

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2 2 18 7.5]);

ax = axes(fig);

imagesc(ax, heatmapMat);

% Mean correlations in top row
for e = 1:nEdges
    text(ax, e, 1, sprintf('%.2f', meanCorr_sorted(e)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 6, ...
        'FontWeight', 'bold', ...
        'Color', 'k');
end

colormap(ax, parula);
clim(ax, [-1 1]);

cb = colorbar(ax);
cb.FontName = 'Arial';
cb.FontSize = 6;

xlabel(ax, 'Network pair', ...
    'FontName', 'Arial', ...
    'FontSize', 7, ...
    'FontWeight', 'bold');

ylabel(ax, 'Study', ...
    'FontName', 'Arial', ...
    'FontSize', 7, ...
    'FontWeight', 'bold');

xticks(ax, 1:nEdges);
xticklabels(ax, edgeLabels_plot);
xtickangle(ax, 45);

yticks(ax, 1:size(heatmapMat,1));
yticklabels(ax, heatmapRowLabels);

title(ax, ...
    'Pairwise correlations between resting-state network responses across studies', ...
    'FontName', 'Arial', ...
    'FontSize', 8, ...
    'FontWeight', 'bold');

hold(ax, 'on');
yline(ax, 1.5, 'k-', 'LineWidth', 1.2);

% Tick labels
set(ax, ...
    'TickLength', [0 0], ...
    'FontName', 'Arial', ...
    'FontSize', 6, ...
    'Box', 'off');

%% Export as SVG
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'Figure3b.svg', '-dsvg');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global signal effects

ER_globalSignal = nanmean(emoregDat.dat)';

%% =========================================================
% Plot global signal across SUBDOMAINS
% - Emotion regulation shown in black on the left
% - Kragel studies are first cleaned and scaled within study
% - Then the two studies belonging to the same subdomain are combined
% - Points, mean, and 95% CI all share domain color
% - X-axis labels show subdomains
%% =========================================================

%% Storage
nDatasets = size(kragelMetaDat_unique, 1);
globalSignal_all = cell(nDatasets,1);

%% 1) Loop over studies: compute global signal, remove outliers, scale by study SD
for i = 1:nDatasets
    
    kragelDat_temp = get_wh_image(kragelDat, (kragelMetaDat.Studynumber == i));
    
    % Global signal per subject
    y = nanmean(kragelDat_temp.dat)';   
    y = y(:);
    y = y(~isnan(y));
    
    % Remove +/- 1.5 IQR outliers
    q1 = prctile(y, 25);
    q3 = prctile(y, 75);
    iqr_val = q3 - q1;
    
    lower = q1 - 1.5 * iqr_val;
    upper = q3 + 1.5 * iqr_val;
    y = y(y >= lower & y <= upper);
    
    % Scale by study SD (without centering, so zero stays meaningful)
    y_sd = std(y, 'omitnan');
    if numel(y) > 1 && y_sd > 0
        y = y ./ y_sd;
    end
    
    globalSignal_all{i} = y;
end

%% 2) Clean and scale emotion regulation dataset
yER = ER_globalSignal(:);
yER = yER(~isnan(yER));

q1 = prctile(yER, 25);
q3 = prctile(yER, 75);
iqr_val = q3 - q1;

lower = q1 - 1.5 * iqr_val;
upper = q3 + 1.5 * iqr_val;
yER = yER(yER >= lower & yER <= upper);

yER_sd = std(yER, 'omitnan');
if numel(yER) > 1 && yER_sd > 0
    yER = yER ./ yER_sd;
end

%% 3) Build subdomain-level combined datasets
domains = string(kragelMetaDat_unique.Domain);
domains = domains(:);

subdomains = string(kragelMetaDat_unique.Subdomain);
subdomains = subdomains(:);

uniqueSubdomains = unique(subdomains, 'stable');
uniqueSubdomains = uniqueSubdomains(:);
nSubdomains = numel(uniqueSubdomains);

combinedSignal_all = cell(nSubdomains,1);
subdomainDomain = strings(nSubdomains,1);

for s = 1:nSubdomains
    idx = find(subdomains == uniqueSubdomains(s));
    
    % Combine the already-scaled study-level vectors
    y_combined = [];
    for j = 1:numel(idx)
        y_combined = [y_combined; globalSignal_all{idx(j)}];
    end
    
    combinedSignal_all{s} = y_combined;
    
    % Domain for coloring: take the first matching study
    subdomainDomain(s) = domains(idx(1));
end

%% 4) Domain colors
uniqueDomains = unique(domains, 'stable');
uniqueDomains = uniqueDomains(:);
nDomains = numel(uniqueDomains);

domainColors = lines(nDomains);

subdomainColors = zeros(nSubdomains, 3);
for s = 1:nSubdomains
    idx = find(uniqueDomains == subdomainDomain(s), 1);
    subdomainColors(s,:) = domainColors(idx,:);
end

% Clean legend labels: replace underscores with spaces
legendDomainLabels = replace(uniqueDomains, "_", " ");

%% 5) Plot
figure;
hold on;

% Grey dashed reference line at y = 0
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

jitter = 0.18;
pointSize = 25;
meanSize = 120;
capWidth = 0.15;

%% Emotion regulation dataset at x = 1
x0 = 1;
xj = x0 + (rand(size(yER)) - 0.5) * 2 * jitter;

scatter(xj, yER, pointSize, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

n = numel(yER);
m = mean(yER, 'omitnan');
sem = std(yER, 'omitnan') / sqrt(n);
ci = tinv(0.975, n - 1) * sem;

plot([x0 x0], [m - ci, m + ci], 'k', 'LineWidth', 2);
plot([x0 - capWidth, x0 + capWidth], [m - ci, m - ci], 'k', 'LineWidth', 2);
plot([x0 - capWidth, x0 + capWidth], [m + ci, m + ci], 'k', 'LineWidth', 2);
scatter(x0, m, meanSize, 'k', 'filled');

%% Subdomains at x = 2...nSubdomains+1
for s = 1:nSubdomains
    
    y = combinedSignal_all{s};
    x0 = s + 1;
    col = subdomainColors(s,:);
    
    xj = x0 + (rand(size(y)) - 0.5) * 2 * jitter;
    
    scatter(xj, y, pointSize, col, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
    
    n = numel(y);
    m = mean(y, 'omitnan');
    sem = std(y, 'omitnan') / sqrt(n);
    ci = tinv(0.975, n - 1) * sem;
    
    plot([x0 x0], [m - ci, m + ci], '-', 'Color', col, 'LineWidth', 2);
    plot([x0 - capWidth, x0 + capWidth], [m - ci, m - ci], '-', 'Color', col, 'LineWidth', 2);
    plot([x0 - capWidth, x0 + capWidth], [m + ci, m + ci], '-', 'Color', col, 'LineWidth', 2);
    scatter(x0, m, meanSize, col, 'filled');
end

%% 6) Axis labels
xticks(1:(nSubdomains + 1));
xticklabels(["Emotion regulation"; uniqueSubdomains(:)]);
xtickangle(45);

ylabel('Global signal scaled by study SD');
title('Global signal across subdomains');

xlim([0.5, nSubdomains + 1.5]);
box off;

%% 7) Legend
legendHandles = gobjects(nDomains + 1, 1);
legendHandles(1) = scatter(nan, nan, 60, 'k', 'filled');

for d = 1:nDomains
    legendHandles(d + 1) = scatter(nan, nan, 60, domainColors(d,:), 'filled');
end

legendLabels = ["Emotion regulation"; legendDomainLabels(:)];
legend(legendHandles, legendLabels, 'Location', 'eastoutside');