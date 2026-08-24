function outStructure = fmri_meta_analysis_fixed(InputData)
    
    % Create statistic image to fill with fixed-effect meta-analysis results
    statImageEffect = ttest(InputData);
    statImageEffect.type = 'p';

    sprintf("Processing...")
    fprintf('Iteration: ');

    for i = 1:size(InputData.dat, 1)
            
        fprintf('%d', i);
        fprintf(repmat('\b', 1, numel(num2str(i))));
        
        % Subset values unequal 0, indicating voxel not included in image
        rowInd = ~(isnan(InputData.dat(i,:)) | InputData.dat(i,:) == 0);
        
        % Perform fixed-effect meta-analysis and save results
        voxelWiseResults = fixed_effects_meta_analysis( ...
            InputData.dat(i,rowInd)', ...
            InputData.covariates(rowInd) ...
        );

        statImageEffect.dat(i, 1) = voxelWiseResults.pooled_corr;
        statImageEffect.p(i, 1) = voxelWiseResults.p_value;
        statImageEffect.ste(i, 1) = voxelWiseResults.pooled_se;
        statImageEffect.N(i, 1) = voxelWiseResults.df + 1;
            
    end

    outStructure = struct();
    outStructure.effectSizeImage = statImageEffect;

end