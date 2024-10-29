% Gist: Meta-Analysis works. There is a rounding error on the 7th decimal
% place for t-values when assigning them to the fmri_data object.
% Also, be careful when using canlab functionality:
% The effect size image output is a p-type statistic_image object,
% but the field for 'dat' contains correlations. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correspondence with meta-analysis in R

% load example data
testDat = readtable('mcdanielsTestData.csv')

% convert r to t
testDat.ti = testDat.ri  .* sqrt(testDat.ni - 2) ./ sqrt(1 - testDat.ri .^2)
testDat_short = testDat(11:15,:)

% run core meta-analysis function on first 5 rows for decimal precision check (identical!)
random_effects_meta_analysis(testDat_short.ti, testDat_short.ni)
%random_effects_meta_analysis_psychomet(testDat_short.ti, testDat_short.ni, relX = 0.8, relY = 0.8, restrFactorX=1.5, restrFactorY=1.2)

% load fmri data
cd('../../data/tMaps_resliced/Questionnaire')
image_path_test = filenames(fullfile(pwd, char("*.nii")), 'absolute');
fmriData_test = fmri_data(image_path_test(1:5, 1));

% mask fmri data
cd('..\..\masks')
amyMask = fmri_data('NIDER_AmygdalaROI.nii')
fmriData_test_amyMasked = apply_mask(fmriData_test, amyMask)
cd('..\..\testAnalyses\testMetaAnalysisCode')

% combine test data with fmri data and add sample size info
repEffects = transpose(repmat(testDat_short.ti, 1, 100));

% Step 2: Determine the number of zeros to insert
numZeros = 20;  % Example number of zeros to insert
% Step 3: Randomly select positions
numElements = length(repEffects);
rng(10)
randomPositions = randperm(numElements, numZeros);
% Step 4: Insert the zeros
resultArray = repEffects;
resultArray(randomPositions) = 0;
repEffects = resultArray

fmriData_test_amyMasked.dat(201:300,:) = repEffects;
fmriData_test_amyMasked.covariates = testDat_short.ni;

% perform random_effects_meta_analysis from fmri_data object
random_effects_meta_analysis(testDat_short.ti, testDat_short.ni)
random_effects_meta_analysis(fmriData_test_amyMasked.dat(201,:)', fmriData_test_amyMasked.covariates)
testDat_short.ti - fmriData_test_amyMasked.dat(201,:)' % deviates after 7 decimals
% reason: assignment to dataframe appears to round differently
testDat_short.ti - repEffects(1,:)'
testDat_short.ti - fmriData_test_amyMasked.dat(201,:)'

% perform fmri_meta_analysis
random_effects_meta_analysis(fmriData_test_amyMasked.dat(202,:)', fmriData_test_amyMasked.covariates)
random_effects_meta_analysis(fmriData_test_amyMasked.dat(203,:)', fmriData_test_amyMasked.covariates)
McDonald_MA_struct = fmri_meta_analysis(fmriData_test_amyMasked);
McDonald_MA_struct.effectSizeImage.dat(202,:)

testEffect = McDonald_MA_struct.effectSizeImage;
testVariance = McDonald_MA_struct.tauImage;

vertcat(testEffect.dat(203), testEffect.p(203), testEffect.N(203), testEffect.ste(203))
vertcat(testVariance.dat(203), testVariance.p(203), testEffect.N(203))
orthviews(testEffect)

% save and reload
write(testEffect, 'fname', 'statisticImageTest.nii')
testEffect_reload = statistic_image('statisticImageTest.nii')
orthviews(testEffect_reload)

% save and reload after thresholding 
testEffect_thresh = threshold(testEffect, 0.05, 'fdr')
orthviews(testEffect_thresh)
write(testEffect_thresh, 'thresh', 'fname', 'statisticImageTest_thresh.nii')
testEffect_reload_thresh = statistic_image('statisticImageTest_thresh.nii')
orthviews(testEffect_reload_thresh)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meta-Analysis of pain images

%%%%%%%%%%%%%%%%%%%%%%%%
% load data
cd('painMaps_neurovault')
image_path_test = filenames(fullfile(pwd, char("*.nii.gz")), 'absolute');

cd('..\..\..\data\masks')
mni_templ = fmri_data('MNI152_T1_2mm_brain_mask.nii');
cd('..\..\testAnalyses\testMetaAnalysisCode\painMaps_neurovault\')

image_names = filenames(fullfile(pwd, char("*.nii")), 'absolute');

% for i=1:length(image_names)
% 
%     [~,name,~] = fileparts(image_names{i});
%     name
%     scn_map_image(image_names{i}, mni_templ, 'write', strcat(name, '_resliced.nii'));
% 
% end

image_names_resliced = filenames(fullfile(pwd, '*_resliced.nii'), 'absolute')
painData = fmri_data(image_names_resliced)


%%%%%%%%%%%%%%%%%%%%%%%%
% merge fmri_data with sample size information
cd('../')
sampleSizes = readtable("studyInfo_neurovault.xlsx");
sampleSizes = sampleSizes(:, {'imageNameStem', 'sample_size'});
sampleSizes = unique(sampleSizes)

studyNames = cellstr(painData.image_names)
studyNames = cell2table(studyNames, "VariableNames", ["imageNameStem"]);

sampleSize_fMRIOrder = join(studyNames, sampleSizes);
painData.covariates = sampleSize_fMRIOrder.sample_size;


%%%%%%%%%%%%%%%%%%%%%%%%
% Quality Control

descriptives(painData);

% Plots
plot(painData)

orthviews(painData)

histogram(painData)
histogram(painData, 'byimage', 'color', 'b');

% table of quality measures
qualityTable = array2table(zeros(size(painData.image_names, 1), 2)); % Initialize with zeros, change as needed
qualityTable.Properties.VariableNames = ["NumNonNullVoxels", "NumDiffValues"]
namesTable = array2table(string(painData.image_names));
namesTable.Properties.VariableNames = "studyNames";

for i = 1:size(qualityTable,1)

    qualityTable.NumNonNullVoxels(i) = sum(painData.dat(:,i) ~= 0);
    qualityTable.NumDiffValues(i) = size(unique(painData.dat(:,i)), 1);

end
qualityTable_pain = horzcat(namesTable, qualityTable)

format("long")
qualityTable_pain


%%%%%%%%%%%%%%%%%%%%%%%%
% conduct MA on six images of correct data type

% subset images
painData_clean = get_wh_image(painData, logical([0,1,1,1,1,1,1]));
painData_clean.image_names
orthviews(painData_clean)

% perform MA
painMA_results = fmri_meta_analysis(painData_clean)
painMA_results2 = fmri_meta_analysis_psychomet(painData_clean)
save("painMA.mat", painMA_results)
load("painMA.mat")


effectSize_image= painMA_results.effectSizeImage
tauImage = painMA_results.tauImage

% super important to mask the tau image for the minimum number of studies
% bc otherwise can return inf when min studies = 1 vor a voxel.
% can fix this in future versions
tauImage = mask_min_studies(tauImage, 3)

% plot results
histogram(effectSize_image.dat)
histogram(effectSize_image.p)
histogram(effectSize_image.N)
histogram(tauImage.dat)
histogram(tauImage.p)
orthviews(effectSize_image)

% subset voxels covered by at least 3 studies
maskMinN = amyMask
maskMinN.dat = zeros(size(maskMinN.dat,1), 1)
maskMinN.dat(effectSize_image.N >= 5) = 1
effectSize_image_minFive = apply_mask(effectSize_image, maskMinN)
orthviews(effectSize_image_minFive)

% threshold image
painResults_fdr05 = threshold(effectSize_image_minFive, 0.05, 'fdr')
orthviews(painResults_fdr05)

% test results with neurosynth decoding
[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data(painResults_fdr05), 'images_are_replicates', false, 'noverbose');

% plot against meta-analytic pain image from canlab (Zunhammer study)
painMAcanlab = fmri_data('full_pain_g_p_FDR.nii.gz')


o2 = canlab_results_fmridisplay([], 'multirow', 2);

o2 = montage(region(painResults_fdr05), o2, 'wh_montages', 1:2, 'colormap');
o2 = title_montage(o2, 2, 'Matlab random effects Meta-Analysis (q < .05 corrected)');

o2 = montage(region(painMAcanlab), o2, 'wh_montages', 3:4, 'colormap');
o2 = title_montage(o2, 4, 'Zunhammer Meta-Analysis (q < .05 corrected)');


[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(painMAcanlab, 'images_are_replicates', false, 'noverbose');




