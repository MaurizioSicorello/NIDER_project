
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
writetable(NOI_outTable, '../../results/bayesfactors/NOI_quest.csv');

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

