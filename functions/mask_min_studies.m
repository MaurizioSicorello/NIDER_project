function maskedImage = mask_min_studies(statisticImage, minStudies)

    maskMinN = fmri_data(statisticImage);
    maskMinN.dat = zeros(size(maskMinN.dat,1), 1);
    maskMinN.dat(statisticImage.N >= minStudies) = 1;
    maskedImage = apply_mask(statisticImage, maskMinN);

end



