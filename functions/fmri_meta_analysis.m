% Takes an fmri_data object as input, which must contain the image-wise
% sample sizes in the "covariates" substructure.
% Optional argument sets the minimum number of studies, which have valid
% values for a voxel for the meta-analysis to be conducted on that voxel.
% Defaults to minStudies = 1.
% Returns two statistical images. First image contains correlations in the
% dat substructure as well as p-values and voxel-wise sample sizes.
% Second image contains estimates of tau in the dat substructure as well
% as p-values for excess variance

% Issue: the core function "random_effects_meta_analysis" checks the dimensions
% of (a) number of images and (b) number of sample size. But it shows
% no error when wrapped by the "fmri_meta_analysis" function.
% Probably due to subsetting [line 33]


function outStructure = fmri_meta_analysis(InputData)
    
    % create statistic image to fill with results
    statImageEffect = ttest(InputData);
    statImageEffect.type = 'p';
    statImageVariance = statImageEffect;

    sprintf("Processing...")
    fprintf('Iteration: ');

    for i = 1:size(InputData.dat, 1)
            
            fprintf('%d', i);
            fprintf(repmat('\b', 1, numel(num2str(i))));
            
            % subset values unequal 0 (indicates voxel not included in image)
            rowInd = ~(isnan(InputData.dat(i,:)) | InputData.dat(i,:) == 0);
            
            % perform MA and save results
            voxelWiseResults = random_effects_meta_analysis(InputData.dat(i,rowInd)', InputData.covariates(rowInd));

            statImageEffect.dat(i, 1) = voxelWiseResults.pooled_corr;
            statImageEffect.p(i, 1) = voxelWiseResults.p_value;
            statImageEffect.ste(i,1) = voxelWiseResults.pooled_se;
            statImageEffect.N(i,1) = voxelWiseResults.df+1;   

            statImageVariance.dat(i, 1) = voxelWiseResults.tau;
            statImageVariance.p(i, 1) = voxelWiseResults.Q_pvalue;
            statImageVariance.N(i,1) = voxelWiseResults.df+1;
            
            % fmri_data objects cause problems when values are exactly.
            % This is a problem for tau, as it is set to 0 when negative.
            if voxelWiseResults.tau == 0
                statImageVariance.dat(i, 1) = 0.0000000000001;
            end

            
    end

    outStructure = struct();
    outStructure.effectSizeImage = statImageEffect;
    outStructure.tauImage = statImageVariance;

end

