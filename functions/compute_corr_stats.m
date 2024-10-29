function [mean_val, std_val] = compute_corr_stats(fMRIdata, atlas_obj)

    % Initialize roiTab with zeros
    roiTab = zeros(size(fMRIdata.dat, 2), size(atlas_obj.labels, 2));

    % Loop through each atlas label
    for i = 1:size(atlas_obj.labels, 2)
        try
            % Extract ROI averages for the current atlas region
            roiTab(:,i) = extract_roi_averages(fMRIdata, select_atlas_subset(atlas_obj, i)).dat;
        catch
            % If there's an error, skip to the next iteration
            warning('Error occurred with atlas region %d. Skipping to next.', i);
            continue;
        end
    end

    % Remove columns with all zeros
    non_zero_columns = any(roiTab ~= 0, 1);
    roiTab = roiTab(:, non_zero_columns);

    % Compute the correlation matrix
    corrMat = corr(roiTab);

    % Extract the lower triangular part of the correlation matrix
    lower_tri = tril(corrMat, -1);
    lower_tri_array = lower_tri(lower_tri ~= 0);

    % Calculate the mean and standard deviation of the lower triangular elements
    mean_val = mean(lower_tri_array);
    std_val = std(lower_tri_array);

end