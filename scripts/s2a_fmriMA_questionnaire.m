
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL SECTION

% load MNI template
cd('..\data\masks')
MNIimage = fmri_data('MNI152_T1_2mm_brain_mask.nii')

% load grey matter masks
gray_mask_sparse = fmri_mask_image('gray_matter_mask_sparse.img');
gray_mask = fmri_mask_image('gray_matter_mask.img');
orthviews(gray_mask_sparse)
orthviews(gray_mask)

% load and combine pregistered meta-analytic ROIs
cd('morawetzMA')
morawetzOne = fmri_data('EmotionRegulation_Cluster_1_of_4_cFWE05_001_104.nii')
morawetzTwo = fmri_data('EmotionRegulation_Cluster_2_of_4_cFWE05_001_205.nii')
morawetzThree = fmri_data('EmotionRegulation_Cluster_3_of_4_cFWE05_001_118.nii')
morawetzFour = fmri_data('EmotionRegulation_Cluster_4_of_4_cFWE05_001_114.nii')

netMaskPrereg = image_math(morawetzOne, morawetzTwo, 'plus')

netMaskAll = image_math(netMaskPrereg, morawetzThree, 'plus')
netMaskAll = image_math(netMaskAll, morawetzFour, 'plus')

cd('..\..\..\scripts\')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREREGISTERED QUESTIONNAIRE META-ANALYSIS


%%%%%%%%%%%%%%%%%%%%%%%%
% load data
cd('../data/tMaps_resliced/Questionnaire')
image_path_quest = filenames(fullfile(pwd, char("*.nii")), 'absolute');
fmriData_quest = fmri_data(image_path_quest);


%%%%%%%%%%%%%%%%%%%%%%%%
% Combine studies with >1 effect sizes (currently only Morawetz 2016b

% Subset study names
studyNames_quest = cellstr(fmriData_quest.image_names)
studyNames_quest = cellfun(@(x) regexp(x, '_', 'split'), studyNames_quest, 'UniformOutput', false);
studyNames_quest = cellfun(@(x) x{1}, studyNames_quest, 'UniformOutput', false);

% Identify duplicates
[uniqueValues, indexA, ~] = unique(studyNames_quest);
logicalIndex = true(size(studyNames_quest));
logicalIndex(indexA) = false;
duplicates_quest = studyNames_quest(logicalIndex)

% Replace study-duplicate-images with averaged images
for i = 1:size(duplicates_quest, 1)

    % Create average image
    fmriSubset_temp = get_wh_image(fmriData_quest, ismember(studyNames_quest, duplicates_quest(i)));
    fmriSubset_Mean_temp = mean(fmriSubset_temp);
    fmriSubset_Mean_temp.image_names = append(duplicates_quest{i}, "_quest_resliced.nii");

    % Replace duplicates with average
    fmriData_quest = get_wh_image(fmriData_quest, ~ismember(studyNames_quest, duplicates_quest(i)));
    fmriData_quest = cat(fmriData_quest, fmriSubset_Mean_temp);

    % get updated study names
    studyNames_quest = cellstr(fmriData_quest.image_names)
    studyNames_quest = cellfun(@(x) regexp(x, '_', 'split'), studyNames_quest, 'UniformOutput', false);
    studyNames_quest = cellfun(@(x) x{1}, studyNames_quest, 'UniformOutput', false);

end
fmriData_quest.removed_images = 0 % cat sets this ~= 0. Then, masking doesn't work


%%%%%%%%%%%%%%%%%%%%%%%%
% merge fmri_data with sample size information
cd('../../studyInformation/')
studyInfo = readtable("studyInformation_publication.xlsx");
sampleSizes = studyInfo(:, {'imageNameStem', 'sampleSizeQuest'});
sampleSizes = unique(sampleSizes)
sampleSizes.imageNameStem{1} = 'BenzaitUnpublished'

studyNames_quest = cellstr(fmriData_quest.image_names)
studyNames_quest = cellfun(@(x) regexp(x, '_', 'split'), studyNames_quest, 'UniformOutput', false);
studyNames_quest = cellfun(@(x) x{1}, studyNames_quest, 'UniformOutput', false);
studyNames_quest = cell2table(studyNames_quest, "VariableNames", ["imageNameStem"]);

sampleSize_fMRIOrder = join(studyNames_quest, sampleSizes);
fmriData_quest.covariates = sampleSize_fMRIOrder.sampleSizeQuest;

%setdiff(studyNames_quest.imageNameStem, sampleSizes.imageNameStem)

%%%%%%%%%%%%%%%%%%%%%%%%
% Quality Control

descriptives(fmriData_quest);

% Plots
plot(fmriData_quest)


% registration looks bad for images 2 and 39 (a lot).
orthviews(get_wh_image(fmriData_quest, 1:16))
orthviews(get_wh_image(fmriData_quest, 17:32))
orthviews(get_wh_image(fmriData_quest, 33:39))
sampleSize_fMRIOrder([2,4,18,19,20,21,22,23,24,25, 33],:)

histogram(get_wh_image(fmriData_quest, 1:20), 'byimage', 'color', 'b');
histogram(get_wh_image(fmriData_quest, 21:29), 'byimage', 'color', 'b');
histogram(get_wh_image(fmriData_quest, 29:39), 'byimage', 'color', 'b');

% table of quality measures
qualityTable = array2table(zeros(size(fmriData_quest.image_names, 1), 2)); % Initialize with zeros, change as needed
qualityTable.Properties.VariableNames = ["NumNonNullVoxels", "NumDiffValues"]
namesTable = array2table(string(fmriData_quest.image_names));
namesTable.Properties.VariableNames = "studyNames";

for i = 1:size(qualityTable,1)

    qualityTable.NumNonNullVoxels(i) = sum(fmriData_quest.dat(:,i) ~= 0);
    qualityTable.NumDiffValues(i) = size(unique(fmriData_quest.dat(:,i)), 1);

end
qualityTable_quest = horzcat(namesTable, qualityTable)

format("long")
qualityTable_quest


qualityTable_quest(qualityTable_quest.NumDiffValues < 10000,:) % no maps with implausible values



% 39 (averaged morawetz) has extremely low values, considering the effect
% size. Maybe effect of design. Leaving it in.
% marinmorales (19) also problematic, bc they have way too many significant
% results (> 20,000)



%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Meta-Analysis 

%%%%%%%%%%%%%
% clean dataset
qualityInd = ones(1, size(fmriData_quest.image_names, 1))
[table([1:39]'), sampleSize_fMRIOrder]
qualityInd([2, 39]) = 0
qualityInd = logical(qualityInd)
fmriData_quest_clean = get_wh_image(fmriData_quest, qualityInd)
fmriData_quest_clean.covariates = fmriData_quest_clean.covariates(qualityInd)

%%%%%%%%%%%%%
% meta-analysis on average activity in networks of interest

% first network
netOneAverage = apply_mask(fmriData_quest_clean, morawetzOne);
netOneAverage = mean(netOneAverage.dat);
random_effects_meta_analysis(netOneAverage', fmriData_quest_clean.covariates)

% second network
netTwoAverage = apply_mask(fmriData_quest_clean, morawetzTwo);
netTwoAverage = mean(netTwoAverage.dat);
random_effects_meta_analysis(netTwoAverage', fmriData_quest_clean.covariates)

% write data for Bayesian analyses
study = cellstr(fmriData_quest_clean.image_names);
NOI_outTable = table(study, fmriData_quest_clean.covariates, netOneAverage', netTwoAverage', ...
    'VariableNames', {'study','sampleSize','networkOne','networkTwo'});
%writetable(NOI_outTable, '../../results/bayesfactors/NOI_quest.csv');

% perform whole-brain meta-analysis
metaResults_quest = fmri_meta_analysis(fmriData_quest_clean)

% load grey matter masks
gray_mask_sparse = fmri_mask_image('gray_matter_mask_sparse.img');
gray_mask = fmri_mask_image('gray_matter_mask.img');

% mask results
effectSizeImage_quest_grayMasked = apply_mask(metaResults_quest.effectSizeImage, gray_mask)
tauImage_quest_grayMasked = apply_mask(metaResults_quest.tauImage, gray_mask)

% mask for min number of studies
effectSizeImage_quest_grayMasked_minStud = mask_min_studies(effectSizeImage_quest_grayMasked, 20)
tauImage_quest_grayMasked_minStud = mask_min_studies(tauImage_quest_grayMasked, 20)

max(abs(effectSizeImage_quest_grayMasked_minStud.dat))
max(effectSizeImage_quest_grayMasked_minStud.dat)
min(effectSizeImage_quest_grayMasked_minStud.dat)
mean(effectSizeImage_quest_grayMasked_minStud.dat)


%%%%% effect sizes 
% threshold effect sizes
effectSizeImage_quest_fdr05 = threshold(apply_mask(effectSizeImage_quest_grayMasked_minStud, netMaskPrereg), .05, 'fdr')
effectSizeImage_quest_fdr05 = threshold(apply_mask(effectSizeImage_quest_grayMasked_minStud, netMaskAll), .05, 'fdr')
effectSizeImage_quest_fdr05 = threshold(effectSizeImage_quest_grayMasked_minStud, .05, 'fdr')


%%%%% tau (empty voxels are due to tau being exactly equal to 0)
tauImage_quest_fdr05 = threshold(tauImage_quest_grayMasked_minStud, .05, 'fdr')

writetable(table(effectSizeImage_quest_grayMasked_minStud.dat, tauImage_quest_grayMasked_minStud.dat, 'VariableNames', ["effectSizes_quest", "tau_quest"]), 'quest_effectSizes.csv')

%%%%% export Network data for Bayesian analyses and load for inspection
fmriData_netAll = apply_mask(fmriData_quest_clean, netMaskAll)

datVarNames = arrayfun(@num2str, 1:size(fmriData_netAll.dat',2), 'UniformOutput', false);
datTable = array2table(fmriData_netAll.dat', 'VariableNames', datVarNames);
datTable{:,:}(datTable{:,:} == 0) = NaN;
WB_outTable = [table(study, fmriData_quest_clean.covariates), datTable];
WB_outTable.Properties.VariableNames{2} = 'sampleSize';  
writetable(WB_outTable, '../../results/bayesfactors/WB_quest.csv');


% questReg = region(tauImage_quest_fdr05)
% table(questReg)
% 
% maskCleanSigHet = fmri_data(tauImage_quest_fdr05)
% 
% if tauImage_quest_fdr05.sig == 1
%     maskCleanSigHet.dat = 1;
% else
%     maskCleanSigHet.dat = 0;
% end

% amyHetData = extract_roi_averages(fmriData_quest_clean, maskCleanSigHet)
% amyHetOut = table((amyHetData.dat ./ sqrt(amyHetData.dat.^2 + (fmriData_quest_clean.covariates-2))), fmriData_quest_clean.image_names, fmriData_quest_clean.covariates, 'VariableNames',["amyR","image name stem", "sample size quest"])
% cd('..\..\results\tables')
% write(amyHetOut, 'amyQuest_heterogeneity.csv')

%%%%%% jackknife analysis (takes a while to run!)
jackknife_quest = fmri_meta_analysis_jackknife(fmriData_quest_clean, gray_mask = gray_mask, minStudies=20, thresh_type = 'fdr', thresh = 0.05)
cd('../../results/tables')
writetable(jackknife_quest, 'jackknife_questionnaire.csv')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensitivity analyses


%%%%%%%%%%%%%%%%%%%%%%%%
% power analysis
metaPower(0.06, sampleSizes.sampleSizeQuest, 0.1, 0.05/20000, 100)


metaPower(effectSize, sampleSizes, tau, alpha, numTruePositives)


%%%%%%%%%%%%%%%%%%%%%%%%
% repeat analysis on single studies

% load image names and sample sizes
cd('..\..\data\tMaps_resliced\Questionnaire\')
fmriData_quest_single = fmri_data(image_path_quest);
cd('../../studyInformation/')
sampleSizes_single = readtable("studyInformation_publication.xlsx");
sampleSizes_single = sampleSizes_single(:, {'imageNameStem', 'sampleSizeQuest'});
sampleSizes_single = sortrows(sampleSizes_single, "imageNameStem");
cd('..\tMaps_resliced\Questionnaire\')


% correct sample size for second Min study
fmriData_quest_single.image_names(21,:)
minPosition = 21
T_part1 = sampleSizes_single(1:minPosition-1, :);
T_part2 = sampleSizes_single(minPosition:end, :);
sampleSizes_single = [T_part1; T_part1(end,:); T_part2];

pValues_sep_quest = table(fmriData_quest_single.image_names, sampleSizes_single.sampleSizeQuest, zeros(size(fmriData_quest_single.image_names, 1),1), ...
    'VariableNames', {'studyName', 'sampleSize', 'num_significant_fdr05'})

for i = 1:size(fmriData_quest_single.fullpath, 1)

    stat_temp = statistic_image(fmriData_quest_single.image_names(i,:), 'type', 't', 'dfe', sampleSizes_single.sampleSizeQuest(i)-2);
    stat_temp_gray = apply_mask(stat_temp, gray_mask);
    stat_temp_thresh = threshold(stat_temp_gray, 0.05, 'fdr');
    pValues_sep_quest(i,3) = table(sum(stat_temp_thresh.sig));

end

write(pValues_sep_quest, '..\..\..\results\tables\quest_singleStudies.csv')

% 3 studies have unreasonably high number of significant voxels. probably
% wrong contrasts

% Option: Plot the significant voxels from single studies against network
% of interest



%%%%%%%%%%%%%%%%%%%%%%%%
% extract correlation matrices for 10, 100, 1000 and 10,000 voxels

cd('..\..\results\power')
fmriCorrDat = apply_mask(fmriData_quest_clean, effectSizeImage_quest_grayMasked_minStud)
fmriCorrDat = apply_mask(fmriCorrDat, netMaskPrereg)
CorrDat = fmriCorrDat.dat';

% 10
rng(10)
writematrix(CorrDat(:,randsample(size(CorrDat,2), 10)), 'data10voxels.csv')

% 100
rng(100)
writematrix(CorrDat(:,randsample(size(CorrDat,2), 100)), 'data100voxels.csv')

% 1000
rng(1000)
writematrix(CorrDat(:,randsample(size(CorrDat,2), 1000)), 'data1000voxels.csv')

% 10000
rng(10000)
writematrix(CorrDat(:,randsample(size(CorrDat,2), 10000)), 'data10000voxels.csv')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power analysis final

%%%%%%%%%%%%%%%%%%%
% check empirical distribution of tau in ROI network
tau_prereg = apply_mask(tauImage_quest_grayMasked_minStud, netMaskPrereg)

% has an excess of very low heterogeneity voxels. filter out to be more
% conservative (i.e., choose higher heterogeneity)
histogram(tau_prereg.dat)
tau_prereg.dat = tau_prereg.dat(tau_prereg.dat>0.003)
histogram(tau_prereg.dat)
quantile(tau_prereg.dat, [0.25 0.5 0.75]) % 0.03, 0.4, 0.06


%%%%%%%%%%%%%%%%%%%
% check typical cluster size
% pico perez (2016): 27 clusters with 3098 voxels (n=247). 4 largest
% clusters have 1998 voxels. 200 + 953 + 576 + 269. round for nicer
% numbers, and make largest cluster more conservative: 800, 500, 300, 200 (2000 voxels total)


%%%%%%%%%%%%%%%%%%%
% create cluster masks
netMaskPrereg_reg = region(netMaskPrereg)

% ROI 1: target 800 voxels from region 14
reg_800_source = region2fmri_data(netMaskPrereg_reg(14), netMaskPrereg);
reg_800_small = reg_800_source;
idx = find(reg_800_source.dat ~= 0);
reg_800_small.dat(:) = 0;
reg_800_small.dat(idx(1:800)) = 1;
mask_800 = apply_mask(reg_800_source, reg_800_small);

% ROI 2: target 500 voxels from region 11
reg_500_source = region2fmri_data(netMaskPrereg_reg(11), netMaskPrereg);
reg_500_small = reg_500_source;
idx = find(reg_500_source.dat ~= 0);
reg_500_small.dat(:) = 0;
reg_500_small.dat(idx(1:500)) = 1;
mask_500 = apply_mask(reg_500_source, reg_500_small);

% ROI 3: target 300 voxels from region 12
reg_300_source = region2fmri_data(netMaskPrereg_reg(12), netMaskPrereg);
reg_300_small = reg_300_source;
idx = find(reg_300_source.dat ~= 0);
reg_300_small.dat(:) = 0;
reg_300_small.dat(idx(1:300)) = 1;
mask_300 = apply_mask(reg_300_source, reg_300_small);

% ROI 4: target 200 voxels from region 1
reg_200_source = region2fmri_data(netMaskPrereg_reg(1), netMaskPrereg);
reg_200_small = reg_200_source;
idx = find(reg_200_source.dat ~= 0);
reg_200_small.dat(:) = 0;
reg_200_small.dat(idx(1:200)) = 1;
mask_200 = apply_mask(reg_200_source, reg_200_small);

o2 = canlab_results_fmridisplay([], 'montages', 2);

addblobs(o2, region(mask_800), 'color', [1 0 0]);
addblobs(o2, region(mask_500), 'color', [0 0 1]);
addblobs(o2, region(mask_300), 'color', [0 0.6 0]);
addblobs(o2, region(mask_200), 'color', [0.8 0.4 0]);

%%%%%%%%%%%%%%%%%%%
% calculate cluster correlation
roi_800 = extract_roi_averages(fmriData_quest_clean, mask_800);
roi_500 = extract_roi_averages(fmriData_quest_clean, mask_500);
roi_300 = extract_roi_averages(fmriData_quest_clean, mask_300);
roi_200 = extract_roi_averages(fmriData_quest_clean, mask_200);

roi_mat = [roi_800.dat roi_500.dat roi_300.dat roi_200.dat];
roi_names = {'ROI 800','ROI 500','ROI 300','ROI 200'};

R = corr(roi_mat, 'Rows', 'pairwise');

figure('Color','w','Position',[100 100 700 600]);
imagesc(R);
axis square;
colormap(parula);
caxis([-1 1]);
cb = colorbar;
cb.Label.String = 'Correlation (r)';

set(gca, ...
    'XTick', 1:4, ...
    'YTick', 1:4, ...
    'XTickLabel', roi_names, ...
    'YTickLabel', roi_names, ...
    'TickLength', [0 0], ...
    'FontSize', 11);

xtickangle(45);

for i = 1:size(R,1)
    for j = 1:size(R,2)
        text(j, i, sprintf('%.2f', R(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 11, ...
            'FontWeight', 'bold', ...
            'Color', 'k');
    end
end

title('ROI Correlation Matrix');

saveas(gcf, 'roi_correlation_matrix.png');



% perform simulations (shortcut: MA with all and no clusters, fill in no-cluster values for in between. should double speed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAST SIMULATION FOR 4 CLUSTERS
% - one null whole-brain meta-analysis computed once
% - per iteration, only 1 representative voxel per cluster is meta-analyzed
% - cluster results are copied into all voxels of that cluster
% - scenarios: 1, 2, 3, or 4 true clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% PREPARE DATA

fmriData_sim = apply_mask(fmriData_quest_clean, gray_mask);
n_vec = fmriData_sim.covariates(:);
k = length(n_vec);

%%%%%%%%%%%%%%%%%%%
% CORRELATION MATRIX FOR THE 4 CLUSTERS

roi_800 = extract_roi_averages(fmriData_quest_clean, mask_800);
roi_500 = extract_roi_averages(fmriData_quest_clean, mask_500);
roi_300 = extract_roi_averages(fmriData_quest_clean, mask_300);
roi_200 = extract_roi_averages(fmriData_quest_clean, mask_200);

roi_mat = [roi_800.dat roi_500.dat roi_300.dat roi_200.dat];
R = corr(roi_mat, 'Rows', 'pairwise');
R = (R + R') ./ 2;
R = R + eye(4) * 1e-6;

%%%%%%%%%%%%%%%%%%%
% NULL META-ANALYSIS ON FULL FIXED DATA (DO ONCE)

SI0_tmp = fmri_meta_analysis_perm(fmriData_sim);
SI0 = mask_min_studies(SI0_tmp.effectSizeImage, 20);

%%%%%%%%%%%%%%%%%%%
% INDICES

% indices in fmriData_sim space (for representative voxels / simulated study data)
mask_800_inSim = ~mask_800.removed_voxels;
mask_800_inSim = mask_800_inSim(~fmriData_sim.removed_voxels);
idx_800 = find(mask_800_inSim);

mask_500_inSim = ~mask_500.removed_voxels;
mask_500_inSim = mask_500_inSim(~fmriData_sim.removed_voxels);
idx_500 = find(mask_500_inSim);

mask_300_inSim = ~mask_300.removed_voxels;
mask_300_inSim = mask_300_inSim(~fmriData_sim.removed_voxels);
idx_300 = find(mask_300_inSim);

mask_200_inSim = ~mask_200.removed_voxels;
mask_200_inSim = mask_200_inSim(~fmriData_sim.removed_voxels);
idx_200 = find(mask_200_inSim);

rep_800 = idx_800(1);
rep_500 = idx_500(1);
rep_300 = idx_300(1);
rep_200 = idx_200(1);

% indices in SI0 space (for insertion into statistic_image objects)
mask_800_inSI0 = ~mask_800.removed_voxels;
mask_800_inSI0 = mask_800_inSI0(~SI0.removed_voxels);
idx_800_SI0 = find(mask_800_inSI0);

mask_500_inSI0 = ~mask_500.removed_voxels;
mask_500_inSI0 = mask_500_inSI0(~SI0.removed_voxels);
idx_500_SI0 = find(mask_500_inSI0);

mask_300_inSI0 = ~mask_300.removed_voxels;
mask_300_inSI0 = mask_300_inSI0(~SI0.removed_voxels);
idx_300_SI0 = find(mask_300_inSI0);

mask_200_inSI0 = ~mask_200.removed_voxels;
mask_200_inSI0 = mask_200_inSI0(~SI0.removed_voxels);
idx_200_SI0 = find(mask_200_inSI0);

%%%%%%%%%%%%%%%%%%%
% INDICES IN ROI-THRESHRESHOLDED SPACE (for ROI-specific cluster recovery)

SI0_reg = apply_mask(SI0, netMaskPrereg);

mask_800_inREG = ~mask_800.removed_voxels;
mask_800_inREG = mask_800_inREG(~SI0_reg.removed_voxels);
idx_800_REG = find(mask_800_inREG);

mask_500_inREG = ~mask_500.removed_voxels;
mask_500_inREG = mask_500_inREG(~SI0_reg.removed_voxels);
idx_500_REG = find(mask_500_inREG);

mask_300_inREG = ~mask_300.removed_voxels;
mask_300_inREG = mask_300_inREG(~SI0_reg.removed_voxels);
idx_300_REG = find(mask_300_inREG);

mask_200_inREG = ~mask_200.removed_voxels;
mask_200_inREG = mask_200_inREG(~SI0_reg.removed_voxels);
idx_200_REG = find(mask_200_inREG);

%%%%%%%%%%%%%%%%%%%
% SETTINGS

n_iter = 1000;
mu_r   = 0.10;
tau_r  = 0.08;

pow1_reg = nan(n_iter,1);
pow2_reg = nan(n_iter,1);
pow3_reg = nan(n_iter,1);
pow4_reg = nan(n_iter,1);

pow1_all = nan(n_iter,1);
pow2_all = nan(n_iter,1);
pow3_all = nan(n_iter,1);
pow4_all = nan(n_iter,1);

% whole-brain cluster recovery
rec_800_s1_all = nan(n_iter,1);
rec_800_s2_all = nan(n_iter,1);
rec_800_s3_all = nan(n_iter,1);
rec_800_s4_all = nan(n_iter,1);

rec_500_s2_all = nan(n_iter,1);
rec_500_s3_all = nan(n_iter,1);
rec_500_s4_all = nan(n_iter,1);

rec_300_s3_all = nan(n_iter,1);
rec_300_s4_all = nan(n_iter,1);

rec_200_s4_all = nan(n_iter,1);

% ROI-thresholded cluster recovery
rec_800_s1_reg = nan(n_iter,1);
rec_800_s2_reg = nan(n_iter,1);
rec_800_s3_reg = nan(n_iter,1);
rec_800_s4_reg = nan(n_iter,1);

rec_500_s2_reg = nan(n_iter,1);
rec_500_s3_reg = nan(n_iter,1);
rec_500_s4_reg = nan(n_iter,1);

rec_300_s3_reg = nan(n_iter,1);
rec_300_s4_reg = nan(n_iter,1);

rec_200_s4_reg = nan(n_iter,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

for iter = 1:n_iter
    
    fprintf('Iteration %d / %d\n', iter, n_iter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE 4 CORRELATED CLUSTER VALUES PER STUDY
    
    mu_z = atanh(mu_r);
    obs_t_cluster = nan(k, 4);
    
    for j = 1:k
        
        n_j = n_vec(j);
        
        if n_j <= 3
            error('Study %d has n <= 3, Fisher-z sampling variance is not defined.', j);
        end
        
        total_var_j = tau_r^2 + 1 / (n_j - 3);
        Sigma_j = total_var_j * R;
        
        obs_z_j = mvnrnd(repmat(mu_z, 1, 4), Sigma_j, 1);
        obs_r_j = tanh(obs_z_j);
        obs_r_j = min(max(obs_r_j, -0.999999), 0.999999);
        
        obs_t_j = obs_r_j .* sqrt((n_j - 2) ./ (1 - obs_r_j.^2));
        obs_t_cluster(j, :) = obs_t_j;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % META-ANALYZE ONE REPRESENTATIVE VOXEL PER CLUSTER
    
    dat_800 = fmriData_sim;
    dat_800.dat = fmriData_sim.dat(rep_800, :);
    dat_800.dat(:) = obs_t_cluster(:,1)';
    SI_800_tmp = fmri_meta_analysis_perm(dat_800);
    SI_800 = SI_800_tmp.effectSizeImage;
    
    dat_500 = fmriData_sim;
    dat_500.dat = fmriData_sim.dat(rep_500, :);
    dat_500.dat(:) = obs_t_cluster(:,2)';
    SI_500_tmp = fmri_meta_analysis_perm(dat_500);
    SI_500 = SI_500_tmp.effectSizeImage;
    
    dat_300 = fmriData_sim;
    dat_300.dat = fmriData_sim.dat(rep_300, :);
    dat_300.dat(:) = obs_t_cluster(:,3)';
    SI_300_tmp = fmri_meta_analysis_perm(dat_300);
    SI_300 = SI_300_tmp.effectSizeImage;
    
    dat_200 = fmriData_sim;
    dat_200.dat = fmriData_sim.dat(rep_200, :);
    dat_200.dat(:) = obs_t_cluster(:,4)';
    SI_200_tmp = fmri_meta_analysis_perm(dat_200);
    SI_200 = SI_200_tmp.effectSizeImage;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BUILD SCENARIOS STEPWISE FROM SI0
    
    SI1 = SI0;
    SI1 = insert_cluster_fields(SI1, SI_800, idx_800_SI0);
    
    SI2 = SI1;
    SI2 = insert_cluster_fields(SI2, SI_500, idx_500_SI0);
    
    SI3 = SI2;
    SI3 = insert_cluster_fields(SI3, SI_300, idx_300_SI0);
    
    SI4 = SI3;
    SI4 = insert_cluster_fields(SI4, SI_200, idx_200_SI0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THRESHOLD AND SAVE RESULTS
    
    % scenario 1
    fdr_reg = threshold(apply_mask(SI1, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI1, .05, 'fdr');
    pow1_reg(iter) = sum(fdr_reg.sig) > 0;
    pow1_all(iter) = sum(fdr_all.sig) > 0;
    rec_800_s1_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_800_s1_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    
    % scenario 2
    fdr_reg = threshold(apply_mask(SI2, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI2, .05, 'fdr');
    pow2_reg(iter) = sum(fdr_reg.sig) > 0;
    pow2_all(iter) = sum(fdr_all.sig) > 0;
    rec_800_s2_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_500_s2_all(iter) = sum(fdr_all.sig(idx_500_SI0)) > 0;
    rec_800_s2_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    rec_500_s2_reg(iter) = sum(fdr_reg.sig(idx_500_REG)) > 0;
    
    % scenario 3
    fdr_reg = threshold(apply_mask(SI3, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI3, .05, 'fdr');
    pow3_reg(iter) = sum(fdr_reg.sig) > 0;
    pow3_all(iter) = sum(fdr_all.sig) > 0;
    rec_800_s3_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_500_s3_all(iter) = sum(fdr_all.sig(idx_500_SI0)) > 0;
    rec_300_s3_all(iter) = sum(fdr_all.sig(idx_300_SI0)) > 0;
    rec_800_s3_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    rec_500_s3_reg(iter) = sum(fdr_reg.sig(idx_500_REG)) > 0;
    rec_300_s3_reg(iter) = sum(fdr_reg.sig(idx_300_REG)) > 0;
    
    % scenario 4
    fdr_reg = threshold(apply_mask(SI4, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI4, .05, 'fdr');
    pow4_reg(iter) = sum(fdr_reg.sig) > 0;
    pow4_all(iter) = sum(fdr_all.sig) > 0;
    rec_800_s4_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_500_s4_all(iter) = sum(fdr_all.sig(idx_500_SI0)) > 0;
    rec_300_s4_all(iter) = sum(fdr_all.sig(idx_300_SI0)) > 0;
    rec_200_s4_all(iter) = sum(fdr_all.sig(idx_200_SI0)) > 0;
    rec_800_s4_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    rec_500_s4_reg(iter) = sum(fdr_reg.sig(idx_500_REG)) > 0;
    rec_300_s4_reg(iter) = sum(fdr_reg.sig(idx_300_REG)) > 0;
    rec_200_s4_reg(iter) = sum(fdr_reg.sig(idx_200_REG)) > 0;
    
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean(pow1_reg)
mean(pow2_reg)
mean(pow3_reg)
mean(pow4_reg)

mean(pow1_all)
mean(pow2_all)
mean(pow3_all)
mean(pow4_all)

mean(rec_800_s1_all)
mean(rec_800_s2_all)
mean(rec_800_s3_all)
mean(rec_800_s4_all)

mean(rec_500_s2_all)
mean(rec_500_s3_all)
mean(rec_500_s4_all)

mean(rec_300_s3_all)
mean(rec_300_s4_all)

mean(rec_200_s4_all)

mean(rec_800_s1_reg)
mean(rec_800_s2_reg)
mean(rec_800_s3_reg)
mean(rec_800_s4_reg)

mean(rec_500_s2_reg)
mean(rec_500_s3_reg)
mean(rec_500_s4_reg)

mean(rec_300_s3_reg)
mean(rec_300_s4_reg)

mean(rec_200_s4_reg)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE / APPEND RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = table(pow1_reg, pow2_reg, pow3_reg, pow4_reg, ...
          pow1_all, pow2_all, pow3_all, pow4_all, ...
          rec_800_s1_all, rec_800_s2_all, rec_800_s3_all, rec_800_s4_all, ...
          rec_500_s2_all, rec_500_s3_all, rec_500_s4_all, ...
          rec_300_s3_all, rec_300_s4_all, ...
          rec_200_s4_all, ...
          rec_800_s1_reg, rec_800_s2_reg, rec_800_s3_reg, rec_800_s4_reg, ...
          rec_500_s2_reg, rec_500_s3_reg, rec_500_s4_reg, ...
          rec_300_s3_reg, rec_300_s4_reg, ...
          rec_200_s4_reg);

outfile = sprintf('simulation_results_r%0.2f_t%0.2f.csv', mu_r, tau_r);

if isfile(outfile)
    writetable(T, outfile, 'WriteMode', 'append');
else
    writetable(T, outfile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed effects meta-analysis

cd('../data/tMaps_resliced/revision_ERQcapability')
image_path_quest = filenames(fullfile(pwd, char("*.nii")), 'absolute');
fmriData_ability = fmri_data(image_path_quest);

plot(fmriData_ability)

fmriData_ability.covariates = [163; 105; 176]


% first network
netOneAverage = apply_mask(fmriData_ability, morawetzOne);
netOneAverage = mean(netOneAverage.dat);
fixed_effects_meta_analysis(netOneAverage', fmriData_ability.covariates)

% second network
netTwoAverage = apply_mask(fmriData_ability, morawetzTwo);
netTwoAverage = mean(netTwoAverage.dat);
fixed_effects_meta_analysis(netTwoAverage', fmriData_ability.covariates)


% perform whole-brain meta-analysis
metaResults_quest = fmri_meta_analysis_fixed(fmriData_ability)

% mask results
effectSizeImage_quest_grayMasked = apply_mask(metaResults_quest.effectSizeImage, gray_mask)

% threshold effect sizes
effectSizeImage_quest_fdr05 = threshold(apply_mask(effectSizeImage_quest_grayMasked, netMaskPrereg), .05, 'fdr')
effectSizeImage_quest_fdr05 = threshold(apply_mask(effectSizeImage_quest_grayMasked, netMaskPrereg), .05, 'unc',  'k', 50)
effectSizeImage_quest_fdr05 = threshold(effectSizeImage_quest_grayMasked, .05, 'fdr')
effectSizeImage_quest_fdr05 = threshold(effectSizeImage_quest_grayMasked, .05, 'unc',  'k', 50)

% power analysis
metaPower(0.1336, [163; 105; 176], 0, 0.05, 1)

% export data for bayesian MA
study = cellstr(fmriData_ability.image_names);
NOI_outTable = table(study, fmriData_ability.covariates, netOneAverage', netTwoAverage', ...
    'VariableNames', {'study','sampleSize','networkOne','networkTwo'});
writetable(NOI_outTable, '../../results/bayesfactors/NOI_quest_ability.csv');

fmriData_netAll = apply_mask(fmriData_ability, netMaskAll)
datVarNames = arrayfun(@num2str, 1:size(fmriData_netAll.dat',2), 'UniformOutput', false);
datTable = array2table(fmriData_netAll.dat', 'VariableNames', datVarNames);
datTable{:,:}(datTable{:,:} == 0) = NaN;
WB_outTable = [table(study, fmriData_ability.covariates), datTable];
WB_outTable.Properties.VariableNames{2} = 'sampleSize';  
writetable(WB_outTable, '../../results/bayesfactors/WB_quest_ability.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POWER SIMULATION FOR FIXED-EFFECT FMRI META-ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% PREPARE DATA

fmriData_sim = apply_mask(fmriData_ability, gray_mask);
n_vec = fmriData_sim.covariates(:);
k = length(n_vec);

%%%%%%%%%%%%%%%%%%%
% CORRELATION MATRIX FOR THE 4 CLUSTERS

roi_800 = extract_roi_averages(fmriData_ability, mask_800);
roi_500 = extract_roi_averages(fmriData_ability, mask_500);
roi_300 = extract_roi_averages(fmriData_ability, mask_300);
roi_200 = extract_roi_averages(fmriData_ability, mask_200);

roi_mat = [roi_800.dat roi_500.dat roi_300.dat roi_200.dat];

R = corr(roi_mat, 'Rows', 'pairwise');
R = (R + R') ./ 2;

% Small ridge to avoid numerical problems in mvnrnd
R = R + eye(4) * 1e-6;

%%%%%%%%%%%%%%%%%%%
% NULL META-ANALYSIS ON FULL FIXED DATA, DO ONCE

SI0_tmp = fmri_meta_analysis_fixed(fmriData_sim);
SI0 = SI0_tmp.effectSizeImage;

% Optional, only if this function is still meaningful for your current data.
% With only k = 3 studies, do not use 20 as minimum number of studies.
% For example:
% SI0 = mask_min_studies(SI0_tmp.effectSizeImage, 2);
%
% Or simply keep:
% SI0 = SI0_tmp.effectSizeImage;

%%%%%%%%%%%%%%%%%%%
% INDICES

% indices in fmriData_sim space, for representative voxels / simulated study data
mask_800_inSim = ~mask_800.removed_voxels;
mask_800_inSim = mask_800_inSim(~fmriData_sim.removed_voxels);
idx_800 = find(mask_800_inSim);

mask_500_inSim = ~mask_500.removed_voxels;
mask_500_inSim = mask_500_inSim(~fmriData_sim.removed_voxels);
idx_500 = find(mask_500_inSim);

mask_300_inSim = ~mask_300.removed_voxels;
mask_300_inSim = mask_300_inSim(~fmriData_sim.removed_voxels);
idx_300 = find(mask_300_inSim);

mask_200_inSim = ~mask_200.removed_voxels;
mask_200_inSim = mask_200_inSim(~fmriData_sim.removed_voxels);
idx_200 = find(mask_200_inSim);

rep_800 = idx_800(1);
rep_500 = idx_500(1);
rep_300 = idx_300(1);
rep_200 = idx_200(1);

% indices in SI0 space, for insertion into statistic_image objects
mask_800_inSI0 = ~mask_800.removed_voxels;
mask_800_inSI0 = mask_800_inSI0(~SI0.removed_voxels);
idx_800_SI0 = find(mask_800_inSI0);

mask_500_inSI0 = ~mask_500.removed_voxels;
mask_500_inSI0 = mask_500_inSI0(~SI0.removed_voxels);
idx_500_SI0 = find(mask_500_inSI0);

mask_300_inSI0 = ~mask_300.removed_voxels;
mask_300_inSI0 = mask_300_inSI0(~SI0.removed_voxels);
idx_300_SI0 = find(mask_300_inSI0);

mask_200_inSI0 = ~mask_200.removed_voxels;
mask_200_inSI0 = mask_200_inSI0(~SI0.removed_voxels);
idx_200_SI0 = find(mask_200_inSI0);

%%%%%%%%%%%%%%%%%%%
% INDICES IN ROI-THRESHOLDED SPACE, for ROI-specific cluster recovery

SI0_reg = apply_mask(SI0, netMaskPrereg);

mask_800_inREG = ~mask_800.removed_voxels;
mask_800_inREG = mask_800_inREG(~SI0_reg.removed_voxels);
idx_800_REG = find(mask_800_inREG);

mask_500_inREG = ~mask_500.removed_voxels;
mask_500_inREG = mask_500_inREG(~SI0_reg.removed_voxels);
idx_500_REG = find(mask_500_inREG);

mask_300_inREG = ~mask_300.removed_voxels;
mask_300_inREG = mask_300_inREG(~SI0_reg.removed_voxels);
idx_300_REG = find(mask_300_inREG);

mask_200_inREG = ~mask_200.removed_voxels;
mask_200_inREG = mask_200_inREG(~SI0_reg.removed_voxels);
idx_200_REG = find(mask_200_inREG);

%%%%%%%%%%%%%%%%%%%
% SETTINGS

n_iter = 1000;

% True fixed common effect, on correlation scale
mu_r = 0.15;

% No between-study heterogeneity for the fixed-effect data-generating model.
% This is equivalent to tau_r = 0.
mu_z = atanh(mu_r);

pow1_reg = nan(n_iter,1);
pow2_reg = nan(n_iter,1);
pow3_reg = nan(n_iter,1);
pow4_reg = nan(n_iter,1);

pow1_all = nan(n_iter,1);
pow2_all = nan(n_iter,1);
pow3_all = nan(n_iter,1);
pow4_all = nan(n_iter,1);

% whole-brain cluster recovery
rec_800_s1_all = nan(n_iter,1);
rec_800_s2_all = nan(n_iter,1);
rec_800_s3_all = nan(n_iter,1);
rec_800_s4_all = nan(n_iter,1);

rec_500_s2_all = nan(n_iter,1);
rec_500_s3_all = nan(n_iter,1);
rec_500_s4_all = nan(n_iter,1);

rec_300_s3_all = nan(n_iter,1);
rec_300_s4_all = nan(n_iter,1);

rec_200_s4_all = nan(n_iter,1);

% ROI-thresholded cluster recovery
rec_800_s1_reg = nan(n_iter,1);
rec_800_s2_reg = nan(n_iter,1);
rec_800_s3_reg = nan(n_iter,1);
rec_800_s4_reg = nan(n_iter,1);

rec_500_s2_reg = nan(n_iter,1);
rec_500_s3_reg = nan(n_iter,1);
rec_500_s4_reg = nan(n_iter,1);

rec_300_s3_reg = nan(n_iter,1);
rec_300_s4_reg = nan(n_iter,1);

rec_200_s4_reg = nan(n_iter,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

for iter = 1:n_iter
    
    fprintf('Iteration %d / %d\n', iter, n_iter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE 4 CORRELATED CLUSTER VALUES PER STUDY
    %
    % Fixed-effect data-generating model:
    % all studies share the same true effect mu_z.
    % The only study-level variability is sampling error.
    
    obs_t_cluster = nan(k, 4);
    
    for j = 1:k
        
        n_j = n_vec(j);
        
        if n_j <= 3
            error('Study %d has n <= 3, Fisher-z sampling variance is not defined.', j);
        end
        
        % Fixed-effect sampling variance on Fisher-z scale
        total_var_j = 1 / (n_j - 3);
        
        % Correlated sampling errors across the 4 simulated clusters
        Sigma_j = total_var_j * R;
        
        obs_z_j = mvnrnd(repmat(mu_z, 1, 4), Sigma_j, 1);
        obs_r_j = tanh(obs_z_j);
        obs_r_j = min(max(obs_r_j, -0.999999), 0.999999);
        
        obs_t_j = obs_r_j .* sqrt((n_j - 2) ./ (1 - obs_r_j.^2));
        obs_t_cluster(j, :) = obs_t_j;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % META-ANALYZE ONE REPRESENTATIVE VOXEL PER CLUSTER
    % using fixed-effect meta-analysis
    
    dat_800 = fmriData_sim;
    dat_800.dat = fmriData_sim.dat(rep_800, :);
    dat_800.dat(:) = obs_t_cluster(:,1)';
    SI_800_tmp = fmri_meta_analysis_fixed(dat_800);
    SI_800 = SI_800_tmp.effectSizeImage;
    
    dat_500 = fmriData_sim;
    dat_500.dat = fmriData_sim.dat(rep_500, :);
    dat_500.dat(:) = obs_t_cluster(:,2)';
    SI_500_tmp = fmri_meta_analysis_fixed(dat_500);
    SI_500 = SI_500_tmp.effectSizeImage;
    
    dat_300 = fmriData_sim;
    dat_300.dat = fmriData_sim.dat(rep_300, :);
    dat_300.dat(:) = obs_t_cluster(:,3)';
    SI_300_tmp = fmri_meta_analysis_fixed(dat_300);
    SI_300 = SI_300_tmp.effectSizeImage;
    
    dat_200 = fmriData_sim;
    dat_200.dat = fmriData_sim.dat(rep_200, :);
    dat_200.dat(:) = obs_t_cluster(:,4)';
    SI_200_tmp = fmri_meta_analysis_fixed(dat_200);
    SI_200 = SI_200_tmp.effectSizeImage;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BUILD SCENARIOS STEPWISE FROM SI0
    
    SI1 = SI0;
    SI1 = insert_cluster_fields(SI1, SI_800, idx_800_SI0);
    
    SI2 = SI1;
    SI2 = insert_cluster_fields(SI2, SI_500, idx_500_SI0);
    
    SI3 = SI2;
    SI3 = insert_cluster_fields(SI3, SI_300, idx_300_SI0);
    
    SI4 = SI3;
    SI4 = insert_cluster_fields(SI4, SI_200, idx_200_SI0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THRESHOLD AND SAVE RESULTS
    
    % scenario 1
    fdr_reg = threshold(apply_mask(SI1, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI1, .05, 'fdr');
    
    pow1_reg(iter) = sum(fdr_reg.sig) > 0;
    pow1_all(iter) = sum(fdr_all.sig) > 0;
    
    rec_800_s1_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_800_s1_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    
    % scenario 2
    fdr_reg = threshold(apply_mask(SI2, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI2, .05, 'fdr');
    
    pow2_reg(iter) = sum(fdr_reg.sig) > 0;
    pow2_all(iter) = sum(fdr_all.sig) > 0;
    
    rec_800_s2_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_500_s2_all(iter) = sum(fdr_all.sig(idx_500_SI0)) > 0;
    
    rec_800_s2_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    rec_500_s2_reg(iter) = sum(fdr_reg.sig(idx_500_REG)) > 0;
    
    % scenario 3
    fdr_reg = threshold(apply_mask(SI3, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI3, .05, 'fdr');
    
    pow3_reg(iter) = sum(fdr_reg.sig) > 0;
    pow3_all(iter) = sum(fdr_all.sig) > 0;
    
    rec_800_s3_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_500_s3_all(iter) = sum(fdr_all.sig(idx_500_SI0)) > 0;
    rec_300_s3_all(iter) = sum(fdr_all.sig(idx_300_SI0)) > 0;
    
    rec_800_s3_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    rec_500_s3_reg(iter) = sum(fdr_reg.sig(idx_500_REG)) > 0;
    rec_300_s3_reg(iter) = sum(fdr_reg.sig(idx_300_REG)) > 0;
    
    % scenario 4
    fdr_reg = threshold(apply_mask(SI4, netMaskPrereg), .05, 'fdr');
    fdr_all = threshold(SI4, .05, 'fdr');
    
    pow4_reg(iter) = sum(fdr_reg.sig) > 0;
    pow4_all(iter) = sum(fdr_all.sig) > 0;
    
    rec_800_s4_all(iter) = sum(fdr_all.sig(idx_800_SI0)) > 0;
    rec_500_s4_all(iter) = sum(fdr_all.sig(idx_500_SI0)) > 0;
    rec_300_s4_all(iter) = sum(fdr_all.sig(idx_300_SI0)) > 0;
    rec_200_s4_all(iter) = sum(fdr_all.sig(idx_200_SI0)) > 0;
    
    rec_800_s4_reg(iter) = sum(fdr_reg.sig(idx_800_REG)) > 0;
    rec_500_s4_reg(iter) = sum(fdr_reg.sig(idx_500_REG)) > 0;
    rec_300_s4_reg(iter) = sum(fdr_reg.sig(idx_300_REG)) > 0;
    rec_200_s4_reg(iter) = sum(fdr_reg.sig(idx_200_REG)) > 0;
    
end

toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary_fixed = struct();

summary_fixed.power_reg = [
    mean(pow1_reg)
    mean(pow2_reg)
    mean(pow3_reg)
    mean(pow4_reg)
];

summary_fixed.power_all = [
    mean(pow1_all)
    mean(pow2_all)
    mean(pow3_all)
    mean(pow4_all)
];

summary_fixed.recovery_800_all = [
    mean(rec_800_s1_all)
    mean(rec_800_s2_all)
    mean(rec_800_s3_all)
    mean(rec_800_s4_all)
];

summary_fixed.recovery_500_all = [
    NaN
    mean(rec_500_s2_all)
    mean(rec_500_s3_all)
    mean(rec_500_s4_all)
];

summary_fixed.recovery_300_all = [
    NaN
    NaN
    mean(rec_300_s3_all)
    mean(rec_300_s4_all)
];

summary_fixed.recovery_200_all = [
    NaN
    NaN
    NaN
    mean(rec_200_s4_all)
];

summary_fixed.recovery_800_reg = [
    mean(rec_800_s1_reg)
    mean(rec_800_s2_reg)
    mean(rec_800_s3_reg)
    mean(rec_800_s4_reg)
];

summary_fixed.recovery_500_reg = [
    NaN
    mean(rec_500_s2_reg)
    mean(rec_500_s3_reg)
    mean(rec_500_s4_reg)
];

summary_fixed.recovery_300_reg = [
    NaN
    NaN
    mean(rec_300_s3_reg)
    mean(rec_300_s4_reg)
];

summary_fixed.recovery_200_reg = [
    NaN
    NaN
    NaN
    mean(rec_200_s4_reg)
];

summary_fixed



T = table(pow1_reg, pow2_reg, pow3_reg, pow4_reg, ...
          pow1_all, pow2_all, pow3_all, pow4_all, ...
          rec_800_s1_all, rec_800_s2_all, rec_800_s3_all, rec_800_s4_all, ...
          rec_500_s2_all, rec_500_s3_all, rec_500_s4_all, ...
          rec_300_s3_all, rec_300_s4_all, ...
          rec_200_s4_all, ...
          rec_800_s1_reg, rec_800_s2_reg, rec_800_s3_reg, rec_800_s4_reg, ...
          rec_500_s2_reg, rec_500_s3_reg, rec_500_s4_reg, ...
          rec_300_s3_reg, rec_300_s4_reg, ...
          rec_200_s4_reg);

outfile = sprintf('simulation_results_r15_fixed');

if isfile(outfile)
    writetable(T, outfile, 'WriteMode', 'append');
else
    writetable(T, outfile);
end