%% Minimal demo: load example data and assign sample sizes
% This demo uses pre-resliced example images included in the repository.

%% 1. Load example images

image_names = filenames( ...
    fullfile('painMaps_neurovault', '*_resliced.nii'), ...
    'absolute');

painData = fmri_data(image_names);


%% 2. Add study sample sizes
% Sample sizes are required by the meta-analysis functions.

sampleSizes = readtable('studyInfo_neurovault.xlsx');
sampleSizes = sampleSizes(:, {'imageNameStem', 'sample_size'});
sampleSizes = unique(sampleSizes);

studyNames = cell2table( ...
    cellstr(painData.image_names), ...
    'VariableNames', {'imageNameStem'});

sampleSize_fMRIOrder = join(studyNames, sampleSizes);

painData.covariates = sampleSize_fMRIOrder.sample_size;


%% 3. Select six example studies

painData_clean = get_wh_image( ...
    painData, logical([0,1,1,1,1,1,1]));


%% 4. Run voxel-wise random-effects meta-analysis

painMA_results = fmri_meta_analysis(painData_clean);

effectSize_image = painMA_results.effectSizeImage;
tauImage         = painMA_results.tauImage;


%% 5. Restrict results to adequately covered voxels

tauImage = mask_min_studies(tauImage, 3);

maskMinN = effectSize_image;
maskMinN.dat = double(effectSize_image.N >= 5);

effectSize_image_minFive = apply_mask( ...
    effectSize_image, maskMinN);


%% 6. FDR correction and visualization

painResults_fdr05 = threshold( ...
    effectSize_image_minFive, 0.05, 'fdr');

orthviews(painResults_fdr05);

disp('Demo completed successfully.')