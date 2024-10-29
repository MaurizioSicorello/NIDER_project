% approach from: https://onlinelibrary.wiley.com/doi/10.1111/peps.12122

function results = correctedCorrelations(r, optional_args)
    
    % handle arguments
    arguments
       r double
       optional_args.relX (1,1) double = 1 
       optional_args.relY (1,1) double = 1
       optional_args.restrFactorX (1,1) double = 1
       optional_args.restrFactorY (1,1) double = 1
    end

    % Extract values from optional_args
    relX = optional_args.relX;
    relY = optional_args.relY;
    restrFactorX = optional_args.restrFactorX;
    restrFactorY = optional_args.restrFactorY;

    % correct r for measurement error
    r_disatten = r ./ sqrt(relX * relY);

    % correct r for restricted variance
    uX = 1/restrFactorX;
    uY = 1/restrFactorY;

    uT = sqrt((relX*uX^2)/(1+relX*uX^2-uX^2));
    uP = sqrt((relY*uY^2)/(1+relY*uY^2-uY^2));

    r_corrected = r_disatten*uT*uP + sqrt((1-uT^2)*(1-uP^2));

    results = r_corrected;

end
    