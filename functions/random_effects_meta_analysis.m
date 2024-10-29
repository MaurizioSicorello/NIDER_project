
function result = random_effects_meta_analysis(t_values, sample_sizes)
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
    % Variances
    variances = 1 ./ (sample_sizes - 3);
    % Weight calculation
    weights = 1 ./ variances;
    % Weighted mean of Fisher's Z (fixed effect)
    pooled_z_fixed = sum(weights .* zs) / sum(weights);
    % Q-statistic
    Q_statistic = sum(weights .* (zs - pooled_z_fixed).^2);
    Q_pvalue = 1 - chi2cdf(Q_statistic, dfbetween);
    % Parameter C
    C = sum(weights) - (sum(weights.^2) / sum(weights));
    % Tau-squared
    tau_square = (Q_statistic - dfbetween) / C;
    if tau_square < 0
        tau_square = 0;
    end
    % Random effect weights
    weights_star = 1 ./ (variances + tau_square);
    % Weighted mean of Fisher's Z (random effect)
    pooled_z_random = sum(weights_star .* zs) / sum(weights_star);
    % Standard error, z- and p-value
    SE = sqrt(1 / sum(weights_star));
    Z_statistic = pooled_z_random / SE;
    p_value_oneTailed = 1 - normcdf(abs(Z_statistic));
    p_value_twoTailed = 2 * p_value_oneTailed;
    % Fisher-backwards transformed estimate with CIs
    pooled_z_random_out = tanh(pooled_z_random);
    LLCI = tanh(pooled_z_random - 1.96 * SE); 
    ULCI = tanh(pooled_z_random + 1.96 * SE);
  
    % Store results in a structure
    result.pooled_corr = pooled_z_random_out;
    result.pooled_se = SE;
    result.df = dfbetween;
    result.p_value = p_value_twoTailed;
    result.tau = sqrt(tau_square);
    result.Q_statistic = Q_statistic;
    result.Q_pvalue = Q_pvalue;
end