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
% SELF-RATING META-ANALYSIS


%%%%%%%%%%%%%%%%%%%%%%%%
% load data
cd('../data/tMaps_resliced/SelfReports')
image_path_rating = filenames(fullfile(pwd, char("*.nii")), 'absolute');
fmriData_rating = fmri_data(image_path_rating);


%%%%%%%%%%%%%%%%%%%%%%%%
% Combine studies with >1 effect sizes

% Subset study names
studyNames_rating = cellstr(fmriData_rating.image_names)
studyNames_rating = cellfun(@(x) regexp(x, '_', 'split'), studyNames_rating, 'UniformOutput', false);
studyNames_rating = cellfun(@(x) x{1}, studyNames_rating, 'UniformOutput', false);

% Identify duplicates
[uniqueValues, indexA, ~] = unique(studyNames_rating);
logicalIndex = true(size(studyNames_rating));
logicalIndex(indexA) = false;
duplicates_rating = studyNames_rating(logicalIndex)

% Replace study-duplicate-images with averaged images
for i = 1:size(duplicates_rating, 1)

    % Create average image
    fmriSubset_temp = get_wh_image(fmriData_rating, ismember(studyNames_rating, duplicates_rating(i)));
    fmriSubset_Mean_temp = mean(fmriSubset_temp);
    fmriSubset_Mean_temp.image_names = append(duplicates_rating{i}, "_rating_resliced.nii");

    % Replace duplicates with average
    fmriData_rating = get_wh_image(fmriData_rating, ~ismember(studyNames_rating, duplicates_rating(i)));
    fmriData_rating = cat(fmriData_rating, fmriSubset_Mean_temp);

    % get updated study names
    studyNames_rating = cellstr(fmriData_rating.image_names)
    studyNames_rating = cellfun(@(x) regexp(x, '_', 'split'), studyNames_rating, 'UniformOutput', false);
    studyNames_rating = cellfun(@(x) x{1}, studyNames_rating, 'UniformOutput', false);

end
fmriData_rating.removed_images = 0 % cat sets this ~= 0. Then, masking doesn't work


%%%%%%%%%%%%%%%%%%%%%%%%
% merge fmri_data with sample size information
cd('../../studyInformation/')
sampleSizes = readtable("studyInformation_publication.xlsx");
sampleSizes = sampleSizes(:, {'imageNameStem', 'sampleSizeRating'});
sampleSizes = unique(sampleSizes)

studyNames_rating = cellstr(fmriData_rating.image_names)
studyNames_rating = cellfun(@(x) regexp(x, '_', 'split'), studyNames_rating, 'UniformOutput', false);
studyNames_rating = cellfun(@(x) x{1}, studyNames_rating, 'UniformOutput', false);
studyNames_rating = cell2table(studyNames_rating, "VariableNames", ["imageNameStem"]);

sampleSize_fMRIOrder = join(studyNames_rating, sampleSizes);
fmriData_rating.covariates = sampleSize_fMRIOrder.sampleSizeRating;


%%%%%%%%%%%%%%%%%%%%%%%%
% Quality Control

descriptives(fmriData_rating);

% Plots
plot(fmriData_rating)

% registration looks bad for images 2 and 40 (a lot), 4 (a little)
% 18-25, 33 are missing a bit of the frontal cortex and lower brain
orthviews(get_wh_image(fmriData_rating, 1:20))
orthviews(get_wh_image(fmriData_rating, 21:33))
sampleSize_fMRIOrder([1, 20, 24, 26, 27, 28, 32],:)

histogram(get_wh_image(fmriData_rating, 1:20), 'byimage', 'color', 'b');
histogram(get_wh_image(fmriData_rating, 21:33), 'byimage', 'color', 'b');


% table of quality measures
qualityTable = array2table(zeros(size(fmriData_rating.image_names, 1), 2)); % Initialize with zeros, change as needed
qualityTable.Properties.VariableNames = ["NumNonNullVoxels", "NumDiffValues"]
namesTable = array2table(string(fmriData_rating.image_names));
namesTable.Properties.VariableNames = "studyNames";

for i = 1:size(qualityTable,1)

    qualityTable.NumNonNullVoxels(i) = sum(fmriData_rating.dat(:,i) ~= 0);
    qualityTable.NumDiffValues(i) = size(unique(fmriData_rating.dat(:,i)), 1);

end
qualityTable_rating = horzcat(namesTable, qualityTable)

format("long")
qualityTable_rating


qualityTable_rating(qualityTable_rating.NumDiffValues < 10000,:) % 4 maps have unplausible values

% check whats wrong with images!
qualityTable_rating([1, 24, 26, 27, 28, 32],:)

%%% Berboth (1) and Morawetz (34) normalization did not work well.
% marinmorales (13) added, bc they have way too many study-wise significant results for their sample size
qualityTable_rating([1, 13, 33],:)



%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Meta-Analysis 

% clean dataset
qualityInd = ones(1, size(fmriData_rating.image_names, 1))
qualityInd([1, 33]) = 0 
qualityInd = logical(qualityInd)
fmriData_rating_clean = get_wh_image(fmriData_rating, qualityInd)
fmriData_rating_clean.covariates = fmriData_rating_clean.covariates(qualityInd)

% first network
netOneAverage = apply_mask(fmriData_rating_clean, morawetzOne);
netOneAverage = mean(netOneAverage.dat);
random_effects_meta_analysis(netOneAverage', fmriData_rating_clean.covariates)

% second network
netTwoAverage = apply_mask(fmriData_rating_clean, morawetzTwo);
netTwoAverage = mean(netTwoAverage.dat);
random_effects_meta_analysis(netTwoAverage', fmriData_rating_clean.covariates)

% perform meta-analysis
metaResults_rating = fmri_meta_analysis(fmriData_rating_clean)

% grey matter mask
effectSizeImage_rating_grayMasked = apply_mask(metaResults_rating.effectSizeImage, gray_mask)
tauImage_rating_grayMasked = apply_mask(metaResults_rating.tauImage, gray_mask)

% mask min studies
effectSizeImage_rating_grayMasked_minStud = mask_min_studies(effectSizeImage_rating_grayMasked, 20)
tauImage_rating_grayMasked_minStud = mask_min_studies(tauImage_rating_grayMasked, 20)

% threshold results
effectSizeImage_rating_fdr05 = threshold(apply_mask(effectSizeImage_rating_grayMasked_minStud, netMaskPrereg), .05, 'fdr')
effectSizeImage_rating_fdr05 = threshold(apply_mask(effectSizeImage_rating_grayMasked_minStud, netMaskAll), .05, 'fdr')
effectSizeImage_rating_fdr05 = threshold(effectSizeImage_rating_grayMasked_minStud, .05, 'fdr')

write(effectSizeImage_rating_fdr05, 'fname', '..\..\results\images\statsImage_rating_thresh_withBrehl.nii', 'thresh')
mean(apply_mask(effectSizeImage_rating_grayMasked_minStud, effectSizeImage_rating_fdr05).dat)

r = autolabel_regions_using_atlas(region(effectSizeImage_rating_fdr05));
montage(r, 'regioncenters', 'colormap');
drawnow, snapnow
montage(r)
[~, ~, regTabOut] = table(r)

regTabOut.Sign = regTabOut.maxZ < 0
regTabOut = sortrows(regTabOut, {'Sign', 'modal_label_descriptions', 'Volume'})
regTabOut.Atlas_regions_covered = [];
regTabOut.region_index = [];
regTabOut.Sign = [];
cd('..\..\results\tables\')
writetable(regTabOut, 'regionTable_rating_withBrehl.xlsx');

rating_effectSize = table(effectSizeImage_rating_fdr05.dat, effectSizeImage_rating_fdr05.sig, 'VariableNames', ["effectSizes_rating_withBrehl", "sig_rating_withBrehl"])

[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(effectSizeImage_rating_grayMasked_minStud), 'images_are_replicates', false, 'noverbose');

tauImage_rating_fdr05 = threshold(apply_mask(tauImage_rating_grayMasked_minStud, netMaskAll), .05, 'fdr')
tauImage_rating_fdr05 = threshold(tauImage_rating_grayMasked_minStud, .05, 'fdr')

% jackknife analysis (takes a while to run! results can be read in with the code below)
jackknife_ratings = fmri_meta_analysis_jackknife(fmriData_rating_clean, gray_mask = gray_mask, minStudies=20, thresh_type = 'fdr', thresh = 0.05)
cd('../../results/tables')
writetable(jackknife_ratings, 'jackknife_ratings.csv')






%%%%%%%%%%%%%%%%%%%%%%%%
% repeat analysis on single studies
cd('..\..\data\tMaps_resliced\selfReports')

pValues_sep_rating = table(fmriData_rating.image_names, fmriData_rating.covariates, zeros(size(fmriData_rating.image_names, 1),1), ...
    'VariableNames', {'studyName', 'sampleSize', 'num_significant_fdr05'})

for i = 1:size(fmriData_rating.fullpath, 1)

    stat_temp = statistic_image(fmriData_rating.image_names(i,:), 'type', 't', 'dfe', fmriData_rating.covariates(i)-1);
    stat_temp_gray = apply_mask(stat_temp, gray_mask);
    %stat_temp_grayNet = apply_mask(stat_temp_gray, image_math(morawetzOne, morawetzTwo, 'plus'));
    stat_temp_thresh = threshold(stat_temp_gray, 0.05, 'fdr');
    pValues_sep_rating(i,3) = table(sum(stat_temp_thresh.sig));

end

pValues_sep_rating
write(pValues_sep_rating, '..\..\..\results\tables\rating_singleStudies.csv')


FU_gianaros = statistic_image('Gianaros2020-PIP_rating_resliced.nii', 'type', 't', 'dfe', 175);
FU_gianaros_t = threshold(FU_gianaros, 0.05, 'fdr')
FU_gianaros_t_r = region(FU_gianaros_t)
table(FU_gianaros_t_r)
montage(FU_gianaros_t_r)
[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(FU_gianaros_t), 'images_are_replicates', false, 'noverbose');

FU_brehl = statistic_image('Brehl2020_rating_resliced.nii', 'type', 't', 'dfe', 175);
FU_brehl_t = threshold(FU_brehl, 0.05, 'fdr')
FU_brehl_t_r = region(FU_brehl_t)
table(FU_brehl_t_r)
[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(FU_brehl_t), 'images_are_replicates', false, 'noverbose');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% repeat MA without brehl

% clean dataset
qualityInd = ones(1, size(fmriData_rating.image_names, 1))
qualityInd([1, 2, 33]) = 0 % CURRENTLY WITHOUT BREHL (2)
qualityInd = logical(qualityInd)
fmriData_rating_clean = get_wh_image(fmriData_rating, qualityInd)
fmriData_rating_clean.covariates = fmriData_rating_clean.covariates(qualityInd)

% first network
netOneAverage = apply_mask(fmriData_rating_clean, morawetzOne);
netOneAverage = mean(netOneAverage.dat);
random_effects_meta_analysis(netOneAverage', fmriData_rating_clean.covariates)

% second network
netTwoAverage = apply_mask(fmriData_rating_clean, morawetzTwo);
netTwoAverage = mean(netTwoAverage.dat);
random_effects_meta_analysis(netTwoAverage', fmriData_rating_clean.covariates)

% perform meta-analysis
metaResults_rating = fmri_meta_analysis(fmriData_rating_clean)

% grey matter mask
effectSizeImage_rating_grayMasked = apply_mask(metaResults_rating.effectSizeImage, gray_mask)
tauImage_rating_grayMasked = apply_mask(metaResults_rating.tauImage, gray_mask)

% mask min studies
effectSizeImage_rating_grayMasked_minStud = mask_min_studies(effectSizeImage_rating_grayMasked, 20)
tauImage_rating_grayMasked_minStud = mask_min_studies(tauImage_rating_grayMasked, 20)

% threshold results
effectSizeImage_rating_fdr05 = threshold(apply_mask(effectSizeImage_rating_grayMasked_minStud, netMaskPrereg), .05, 'fdr')
effectSizeImage_rating_fdr05 = threshold(apply_mask(effectSizeImage_rating_grayMasked_minStud, netMaskAll), .05, 'fdr')
effectSizeImage_rating_fdr05 = threshold(effectSizeImage_rating_grayMasked_minStud, .05, 'fdr')

rating_effectSize_withoutBrehl = table(effectSizeImage_rating_fdr05.dat, effectSizeImage_rating_fdr05.sig, 'VariableNames', ["effectSizes_rating_withoutBrehl", "sig_rating_withoutBrehl"])

writetable(rating_effectSize, 'rating_effectSizes_withBrehl.csv')
writetable(rating_effectSize_withoutBrehl, 'rating_effectSizes_withoutBrehl.csv')


write(effectSizeImage_rating_fdr05, 'fname', '..\..\results\images\statsImage_rating_thresh_withoutBrehl.nii', 'thresh')
mean(apply_mask(effectSizeImage_rating_grayMasked_minStud, effectSizeImage_rating_fdr05).dat)

r = autolabel_regions_using_atlas(region(effectSizeImage_rating_fdr05));
montage(r, 'regioncenters', 'colormap');
drawnow, snapnow
montage(r)
[~, ~, regTabOut] = table(r)

regTabOut.Sign = regTabOut.maxZ < 0
regTabOut = sortrows(regTabOut, {'Sign', 'modal_label_descriptions', 'Volume'})
regTabOut.Atlas_regions_covered = [];
regTabOut.region_index = [];
regTabOut.Sign = [];
cd('..\..\results\tables\')
writetable(regTabOut, 'regionTable_rating_withoutBrehl.xlsx');

[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(effectSizeImage_rating_fdr05), 'images_are_replicates', false, 'noverbose');
[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(ttest(fmriData_rating_clean)), 'images_are_replicates', false, 'noverbose');

tauImage_rating_fdr05 = threshold(apply_mask(tauImage_rating_grayMasked_minStud, netMaskAll), .05, 'fdr')
tauImage_rating_fdr05 = threshold(tauImage_rating_grayMasked_minStud, .05, 'fdr')

