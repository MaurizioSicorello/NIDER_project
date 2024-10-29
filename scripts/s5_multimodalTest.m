



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test with emotion regulation data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and prepare data

[emoregDat, ~, imagenames] = load_image_set('emotionreg')

atlas_obj = load_atlas('canlab2018_2mm'); % 489 regions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlation between amygdala and remaining regions

amyROI = fmri_data('..\data\masks\NIDER_AmygdalaROI.nii')
emoregDat.X = extract_roi_averages(emoregDat, amyROI).dat
emoregDat = apply_mask(emoregDat, amyROI, 'invert')

emoregReg = regress(emoregDat)
histogram(emoregReg.t.dat(:,1))

correlations = emoregReg.t.dat(:,1) ./ sqrt(emoregReg.t.dat(:,1).^2 + 33)
mean(correlations)
std(correlations)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %calculate correlation matrix between all regions

compute_corr_stats(emoregDat, atlas_obj)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test with multimodal data

kragelDat = fmri_data('..\data\kragel2018\kragel_2018_nat_neurosci_270_subjects_test_images.nii')
kragelMetaDat = readtable('..\data\kragel2018\kragel_2018_nat_neurosci_270_subjects_test_images_metadata')

kragelMetaDat_unique = unique(kragelMetaDat);

corrOut_kragel = zeros(size(kragelMetaDat_unique, 1), 2)

for i=1:size(kragelMetaDat_unique, 1)
    
    kragelDat_temp = get_wh_image(kragelDat, (kragelMetaDat.Studynumber == i));
    [mTemp sdTemp] = compute_corr_stats(kragelDat_temp, atlas_obj);
    corrOut_kragel(i,1) = mTemp;
    corrOut_kragel(i,2) = sdTemp;
    
end

write(horzcat(kragelMetaDat_unique, array2table(corrOut_kragel, 'VariableNames', {'corrMean', 'corrSD'})), '..\results\tables\kragelAnalyses.csv')

