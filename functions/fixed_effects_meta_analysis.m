function result = fixed_effects_meta_analysis(t_values, sample_sizes)
    % Check if the input sizes match
    if size(t_values,1) ~= size(sample_sizes,1)
        error('Number of t-values and sample sizes must match');
    end
    
    % Remove missing values
    missingInd = isnan(t_values) | isinf(t_values);
    t_values = t_values(~missingInd);
    sample_sizes = sample_sizes(~missingInd);
    
    % Degrees of freedom for 2nd level
    dfbetween = length(t_values) - 1;

    % Degrees of freedom for each study
    dfwithin = sample_sizes - 2;

    % Convert t-values to correlations
    correlations = t_values ./ sqrt(t_values.^2 + dfwithin);

    % Fisher's Z transformation
    zs = atanh(correlations);

    % Variances of Fisher's Z
    variances = 1 ./ (sample_sizes - 3);

    % Fixed-effect inverse-variance weights
    weights = 1 ./ variances;

    % Weighted mean of Fisher's Z, fixed-effect model
    pooled_z_fixed = sum(weights .* zs) / sum(weights);

    % Standard error, z- and p-value
    SE = sqrt(1 / sum(weights));
    Z_statistic = pooled_z_fixed / SE;
    p_value_oneTailed = 1 - normcdf(abs(Z_statistic));
    p_value_twoTailed = 2 * p_value_oneTailed;

    % Fisher-backwards transformed estimate with CIs
    pooled_corr_fixed = tanh(pooled_z_fixed);
    LLCI = tanh(pooled_z_fixed - 1.96 * SE); 
    ULCI = tanh(pooled_z_fixed + 1.96 * SE);

    % Q-statistic for heterogeneity
    Q_statistic = sum(weights .* (zs - pooled_z_fixed).^2);
    Q_pvalue = 1 - chi2cdf(Q_statistic, dfbetween);

    % Store results in a structure
    result.pooled_corr = pooled_corr_fixed;
    result.pooled_z = pooled_z_fixed;
    result.pooled_se = SE;
    result.LLCI = LLCI;
    result.ULCI = ULCI;
    result.df = dfbetween;
    result.z_value = Z_statistic;
    result.p_value = p_value_twoTailed;
    result.Q_statistic = Q_statistic;
    result.Q_pvalue = Q_pvalue;
end