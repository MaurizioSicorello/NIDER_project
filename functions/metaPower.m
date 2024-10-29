function results = metaPower(effectSize, sampleSizes, tau, alpha, numTruePositives)
    
    % based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5590730/

    %%%%%% convert effect size
    % fisher transformation
    effectSize_z = atanh(effectSize);
    
    %%%%%% calculate standard error
    % Variances
    variances = 1 ./ (sampleSizes - 3);
    % fixed effect weights
    weights = 1 ./ variances
    % Random effect weights
    weights_star = 1 ./ (variances + tau^2);
    % Standard error
    SE_random = sqrt(1 / sum(weights_star));
    SE_fixed = sqrt(1 / sum(weights));
    
    %%%%%% non-centrality parameter
    z_random = effectSize_z/SE_random;
    z_fixed = effectSize_z/SE_fixed;
    
    %%%%%% power
    % critical z-value
    z_crit = norminv(1-alpha/2);
    randomPower = 1 + normcdf(-z_crit + z_random) - normcdf(z_crit + z_random);
    fixedPower = 1 + normcdf(-z_crit + z_fixed) - normcdf(z_crit + z_fixed);
    
    %%%%%% probability of at least 1 true positive test
    randomPower = 1 - (1-randomPower)^numTruePositives;
    fixedPower = 1 - (1-fixedPower)^numTruePositives;

    %%%%%% return results
    results.random = randomPower;
    results.fixed = fixedPower;

end