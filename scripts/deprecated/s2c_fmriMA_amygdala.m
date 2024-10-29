% POTENTIAL ANALYSIS: DEMONSTRATE ECOLOGICAL FALLACY BY TESTING WHETHER CONNECTIVITY ANALYSES
% ARE POSITIVE FOR PATHWAYS REPORTED TO BE NEGATIVE

% Load an prepare masks

% MNI image 
cd('..\data\masks')
MNIimage = fmri_data('MNI152_T1_2mm_brain_mask.nii')

% amgdala
amyROI = fmri_data('NIDER_AmygdalaROI.nii')

% load grey matter masks
gray_mask_sparse = fmri_mask_image('gray_matter_mask_sparse.img');
gray_mask = fmri_mask_image('gray_matter_mask.img');

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
% amygdala META-ANALYSIS


%%%%%%%%%%%%%%%%%%%%%%%%
% load data
cd('../data/tMaps_resliced/Amygdala')
image_path_amygdala = filenames(fullfile(pwd, char("*.nii")), 'absolute');
fmriData_amygdala = fmri_data(image_path_amygdala);


%%%%%%%%%%%%%%%%%%%%%%%%
% Combine studies with >1 effect sizes

% Subset study names
studyNames_amygdala = cellstr(fmriData_amygdala.image_names)
studyNames_amygdala = cellfun(@(x) regexp(x, '_', 'split'), studyNames_amygdala, 'UniformOutput', false);
studyNames_amygdala = cellfun(@(x) x{1}, studyNames_amygdala, 'UniformOutput', false);

% Identify duplicates
[uniqueValues, indexA, ~] = unique(studyNames_amygdala);
logicalIndex = true(size(studyNames_amygdala));
logicalIndex(indexA) = false;
duplicates_amygdala = studyNames_amygdala(logicalIndex)

% Replace study-duplicate-images with averaged images
for i = 1:size(duplicates_amygdala, 1)

    % Create average image
    fmriSubset_temp = get_wh_image(fmriData_amygdala, ismember(studyNames_amygdala, duplicates_amygdala(i)));
    fmriSubset_Mean_temp = mean(fmriSubset_temp);
    fmriSubset_Mean_temp.image_names = append(duplicates_amygdala{i}, "_amygdala_resliced.nii");

    % Replace duplicates with average
    fmriData_amygdala = get_wh_image(fmriData_amygdala, ~ismember(studyNames_amygdala, duplicates_amygdala(i)));
    fmriData_amygdala = cat(fmriData_amygdala, fmriSubset_Mean_temp);

    % get updated study names
    studyNames_amygdala = cellstr(fmriData_amygdala.image_names)
    studyNames_amygdala = cellfun(@(x) regexp(x, '_', 'split'), studyNames_amygdala, 'UniformOutput', false);
    studyNames_amygdala = cellfun(@(x) x{1}, studyNames_amygdala, 'UniformOutput', false);

end
fmriData_amygdala.removed_images = 0 % cat sets this ~= 0. Then, masking doesn't work


%%%%%%%%%%%%%%%%%%%%%%%%
% merge fmri_data with sample size information
cd('../../studyInformation/')
sampleSizes = readtable("studyInformation_publication.xlsx");
sampleSizes = sampleSizes(:, {'imageNameStem', 'sampleSizeAmygdala'});
sampleSizes = unique(sampleSizes)
sampleSizes.imageNameStem{1} = 'BenzaitUnpublished'

studyNames_amygdala = cellstr(fmriData_amygdala.image_names)
studyNames_amygdala = cellfun(@(x) regexp(x, '_', 'split'), studyNames_amygdala, 'UniformOutput', false);
studyNames_amygdala = cellfun(@(x) x{1}, studyNames_amygdala, 'UniformOutput', false);
studyNames_amygdala = cell2table(studyNames_amygdala, "VariableNames", ["imageNameStem"]);

sampleSize_fMRIOrder = join(studyNames_amygdala, sampleSizes);
fmriData_amygdala.covariates = sampleSize_fMRIOrder.sampleSizeAmygdala;


%%%%%%%%%%%%%%%%%%%%%%%%
% Quality Control

descriptives(fmriData_amygdala);

% Plots
outlier = plot(fmriData_amygdala)

% check amygdala correlation with itself and correct for sign errors
% (shouldn't correlate to -1, e.g., because many used first eigenvariate,
% but should be negative. some coauthors indicated uncertainty in the direction of the sign)
amySelfCorr = extract_roi_averages(fmriData_amygdala, amyROI).dat;
amyInd = find(amySelfCorr > 0);

for i=1:length(amyInd)

    fmriData_amygdala.dat(:,amyInd(i)) = fmriData_amygdala.dat(:,amyInd(i))*(-1);

end

outlier = plot(fmriData_amygdala)
fmriData_amygdala.image_names(33,:)


% registration looks bad for images 2 (Berboth) and 39
% (morawetz). Two more studies are thresholded. Remove from analysis
% Gaebler/Hofhansel

orthviews(get_wh_image(fmriData_amygdala, 1:20))
orthviews(get_wh_image(fmriData_amygdala, 21:39))
sampleSize_fMRIOrder([8, 14],:)


histogram(get_wh_image(fmriData_amygdala, 1:20), 'byimage', 'color', 'b');
histogram(get_wh_image(fmriData_amygdala, 21:30), 'byimage', 'color', 'b');

% table of quality measures
qualityTable = array2table(zeros(size(fmriData_amygdala.image_names, 1), 2)); % Initialize with zeros, change as needed
qualityTable.Properties.VariableNames = ["NumNonNullVoxels", "NumDiffValues"]
namesTable = array2table(string(fmriData_amygdala.image_names));
namesTable.Properties.VariableNames = "studyNames";

for i = 1:size(qualityTable,1)

    qualityTable.NumNonNullVoxels(i) = sum(fmriData_amygdala.dat(:,i) ~= 0);
    qualityTable.NumDiffValues(i) = size(unique(fmriData_amygdala.dat(:,i)), 1);

end
qualityTable_amygdala = horzcat(namesTable, qualityTable)

format("long")
qualityTable_amygdala


qualityTable_amygdala([1, 34, 29, 30, 20, 13 ],:)



%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Meta-Analysis 

% registration looks bad for images 2 (Berboth), 39
% (morawetz)
% 8 (Gaebler) and 14 (Hofhansel) are empty

% clean dataset
qualityInd = ones(1, size(fmriData_amygdala.image_names, 1))
qualityInd([2, 8, 14, 39]) = 0;
qualityInd = logical(qualityInd)
fmriData_amygdala_clean = get_wh_image(fmriData_amygdala, qualityInd)
fmriData_amygdala_clean.covariates = fmriData_amygdala_clean.covariates(qualityInd)
fmriData_amygdala_clean = apply_mask(fmriData_amygdala_clean, amyROI, 'invert')

% first network
netOneAverage = apply_mask(fmriData_amygdala_clean, morawetzOne);
netOneAverage = mean(netOneAverage.dat);
random_effects_meta_analysis(netOneAverage', fmriData_amygdala_clean.covariates)

% second network
netTwoAverage = apply_mask(fmriData_amygdala_clean, morawetzTwo);
netTwoAverage = mean(netTwoAverage.dat);
random_effects_meta_analysis(netTwoAverage', fmriData_amygdala_clean.covariates)

% perform meta-analysis
metaResults_amygdala = fmri_meta_analysis(fmriData_amygdala_clean)

% grey matter mask
effectSizeImage_amygdala_grayMasked = apply_mask(metaResults_amygdala.effectSizeImage, gray_mask)
tauImage_amygdala_grayMasked = apply_mask(metaResults_amygdala.tauImage, gray_mask)

% mask min studies
effectSizeImage_amygdala_grayMasked_minStud = mask_min_studies(effectSizeImage_amygdala_grayMasked, 20)
tauImage_amygdala_grayMasked_minStud = mask_min_studies(tauImage_amygdala_grayMasked, 20)

% threshold results
effectSizeImage_amygdala_fdr05 = threshold(apply_mask(effectSizeImage_amygdala_grayMasked_minStud, netMaskPrereg), .05, 'fdr')
effectSizeImage_amygdala_fdr05 = threshold(apply_mask(effectSizeImage_amygdala_grayMasked_minStud, netMaskAll), .05, 'fdr')
effectSizeImage_amygdala_fdr05 = threshold(effectSizeImage_amygdala_grayMasked_minStud, .05, 'fdr')

write(effectSizeImage_amygdala_fdr05, 'fname', '..\..\results\images\statsImage_amygdala_thresh.nii', 'thresh')

mean(apply_mask(effectSizeImage_amygdala_grayMasked_minStud, effectSizeImage_amygdala_fdr05).dat)
histogram(apply_mask(effectSizeImage_amygdala_grayMasked_minStud, effectSizeImage_amygdala_fdr05).dat)

%%%%% tau (empty voxels are due to tau being exactly equal to 0)
tauImage_amygdala_fdr05 = threshold(tauImage_amygdala_grayMasked_minStud, .05, 'fdr')

writetable(table(effectSizeImage_amygdala_grayMasked_minStud.dat, tauImage_amygdala_grayMasked_minStud.dat, 'VariableNames', ["effectSizes_quest", "tau_quest"]), 'amygdala_effectSizes.csv')

% jackknife analysis (takes a while to run! results can be read in with the code below)
jackknife_amygdalas = fmri_meta_analysis_jackknife(fmriData_amygdala_clean, gray_mask = gray_mask, minStudies=20, thresh_type = 'fdr', thresh = 0.05)
cd('../../results/tables')
writetable (jackknife_amygdalas, 'jackknife_amygdalas.csv')

jackknife_amygdalas =  readtable('jackknife_amygdalas3.csv')

% histogram(jackknife_amygdalas.numSigEffects - 6, 20)





%%%%%%%%%%%%%%%%%%%%%%%%
% repeat analysis on single studies
cd('..\tMaps_resliced\selfReports')

pValues_sep_amygdala = table(fmriData_amygdala.image_names, fmriData_amygdala.covariates, zeros(size(fmriData_amygdala.image_names, 1),1), ...
    'VariableNames', {'studyName', 'sampleSize', 'num_significant_fdr05'})

for i = 1:size(fmriData_amygdala.fullpath, 1)

    stat_temp = statistic_image(fmriData_amygdala.image_names(i,:), 'type', 't', 'dfe', fmriData_amygdala.covariates(i)-1);
    stat_temp_gray = apply_mask(stat_temp, gray_mask);
    %stat_temp_grayNet = apply_mask(stat_temp_gray, image_math(morawetzOne, morawetzTwo, 'plus'));
    stat_temp_thresh = threshold(stat_temp_gray, 0.05, 'fdr');
    pValues_sep_amygdala(i,3) = table(sum(stat_temp_thresh.sig));

end

find(strcmp('MarinMorales2021_amygdala.nii', cellstr(fmriData_amygdala.image_names)))
find(strcmp('MulejBratec2015_amygdala_resliced.nii', cellstr(fmriData_amygdala.image_names)))
find(strcmp('Steward2021_amygdala_resliced.nii', cellstr(fmriData_amygdala.image_names)))

pValues_sep_amygdala
%steward, mulejbratec and marinmorales have too many significant voxels
%especially considering their sample size

FU_gianaros = statistic_image('Gianaros2020-PIP_amygdala_resliced.nii', 'type', 't', 'dfe', 175);
FU_gianaros_t = threshold(FU_gianaros, 0.05, 'fdr')
FU_gianaros_t_r = region(FU_gianaros_t)
table(FU_gianaros_t_r)
montage(FU_gianaros_t_r)
[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(FU_gianaros_t), 'images_are_replicates', false, 'noverbose');

FU_brehl = statistic_image('Brehl2020_amygdala_resliced.nii', 'type', 't', 'dfe', 175);
FU_brehl_t = threshold(FU_brehl, 0.05, 'fdr')
FU_brehl_t_r = region(FU_brehl_t)
table(FU_brehl_t_r)
[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(FU_brehl_t), 'images_are_replicates', false, 'noverbose');

orthviews(fmri_data('statsImage_rating_thresh_withoutBrehl.nii'))