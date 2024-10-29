function results = fmri_meta_analysis_jackknife(InputData, optional_args)

    arguments
       InputData fmri_data
       optional_args.thresh_type string = "fdr"
       optional_args.thresh (1,1) double = 0.05 
       optional_args.minStudies (1,1) double = 20
       optional_args.gray_mask (1,1) fmri_data = []
       optional_args.masking (1,1) fmri_data = []
    end

    % Extract values from optional_args
    thresh_type = optional_args.thresh_type;
    thresh = optional_args.thresh;
    minStudies = optional_args.minStudies;
    gray_mask = optional_args.gray_mask;
    masking = optional_args.masking;
    
    % check dimensions for correctness
    if (size(InputData.dat,2) + size(InputData.covariates,1) + size(InputData.image_names,1)) ~= size(InputData.dat,2)*3
        error('Number of images, image-names and/or sample sizes must match');
    end

    % empty table for results
    resultsTable = array2table(nan(size(InputData.image_names,1), 6), 'VariableNames', ["imageNames", "numSigEffects", "numSigTau", "average_effect", "lower_limit", "upper_limit"]);
    resultsTable.imageNames = cellstr(InputData.image_names);

    % loop through studies, leaving each study out once
    for i = 1:size(InputData.dat,2)
        
        % print iteration
        fprintf('\n\n Study: %d out of %d \n\n', i, size(InputData.dat,2));
        
        % subset relevant studies
        studyInd = ones(1, size(InputData.image_names, 1));
        studyInd(i) = 0;
        studyInd = logical(studyInd);
        fmriData_temp = get_wh_image(InputData, studyInd);
        fmriData_temp.covariates = fmriData_temp.covariates(studyInd);

        % perform MA
        metaResults_temp = fmri_meta_analysis(fmriData_temp);
        effectSizeImage_temp = metaResults_temp.effectSizeImage;
        tauImage_temp = metaResults_temp.tauImage;

        % gray-matter masking
        if ~isempty(gray_mask)
            effectSizeImage_temp = apply_mask(metaResults_temp.effectSizeImage, gray_mask);
            tauImage_temp = apply_mask(metaResults_temp.tauImage, gray_mask);
        end

        % min-study masking
        effectSizeImage_temp = mask_min_studies(effectSizeImage_temp, minStudies);
        tauImage_temp = mask_min_studies(tauImage_temp, minStudies);

        % optional masking
        if ~isempty(masking)
            effectSizeImage_temp = apply_mask(effectSizeImage_temp, masking);
            tauImage_temp = apply_mask(tauImage_temp, masking);
        end

        % statistical thresholding
        effectSizeImage_temp_thresh = threshold(effectSizeImage_temp, thresh, thresh_type);
        tauImage_temp_thresh = threshold(tauImage_temp, thresh, thresh_type);

        resultsTable(i,2) = table(sum(effectSizeImage_temp_thresh.p <= effectSizeImage_temp_thresh.threshold));
        resultsTable(i,3) = table(sum(tauImage_temp_thresh.p <= tauImage_temp_thresh.threshold));
        
        effectSize_masked = apply_mask(effectSizeImage_temp, effectSizeImage_temp_thresh);
        resultsTable(i,4) = table(mean(effectSize_masked.dat));
        
        if isempty(min(effectSize_masked.dat))
            lower_li = NaN;
            upper_li = NaN;
        else
            lower_li = table(min(effectSize_masked.dat));
            upper_li = table(max(effectSize_masked.dat));

        resultsTable(i,5) = lower_li;
        resultsTable(i,6) = upper_li;


    end

    results = resultsTable;

end


