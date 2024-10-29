

% load MNI template
cd('../data/masks')
mni_templ = fmri_data('MNI152_T1_2mm_brain_mask.nii');


%load images
cd('../tMaps_raw')
image_names = filenames(fullfile(pwd, '*.nii'), 'absolute')

for i=1:length(image_names)
    
    [~,name,~] = fileparts(image_names{i});
    name
    scn_map_image(image_names{i}, mni_templ, 'write', strcat(name, '_resliced.nii'));
    
end

image_names_resliced = filenames(fullfile(pwd, '*_resliced.nii'), 'absolute')


% screen resliced images
cd('../tMaps_resliced/Amygdala')
image_names = filenames(fullfile(pwd, '*.nii'), 'absolute')
data_obj = fmri_data(image_names)
plot(data_obj)
descriptives(data_obj)

data_obj = fmri_data(image_names{1})
hist(data_obj.dat)

cd('../Questionnaire')
image_names = filenames(fullfile(pwd, '*.nii'), 'absolute')
data_obj = fmri_data(image_names)
plot(data_obj)
descriptives(data_obj)


cd('../SelfReports')
image_names = filenames(fullfile(pwd, '*.nii'), 'absolute')
data_obj = fmri_data(image_names)
plot(data_obj)
descriptives(data_obj)




% reslicing led to strange errors for some images. Worked with the
% following alternative procedure
netaQuestResliced = resample_space(fmri_data('NetaUnpublished_quest_ERQ.nii'), mni_templ)
orthviews(netaQuestResliced)
write(netaQuestResliced, 'fname', 'NetaUnpublished_quest_ERQ_resliced.nii', 'overwrite');

powersQuestResliced = resample_space(fmri_data('PowersUnpublished_quest_ERQ.nii'), mni_templ)
orthviews(powersQuestResliced)
write(powersQuestResliced, 'fname', 'PowersUnpublished_quest_ERQ_resliced.nii', 'overwrite');

steinfurthQuestResliced = resample_space(fmri_data('Steinfurth2013_quest_ERQ.nii'), mni_templ)
orthviews(steinfurthQuestResliced)
write(steinfurthQuestResliced, 'fname', 'steinfurthQuestResliced.nii', 'overwrite');

netaRatingResliced = resample_space(fmri_data('NetaUnpublished_rating.nii'), mni_templ)
orthviews(netaRatingResliced)
write(netaRatingResliced, 'fname', 'NetaUnpublished_rating_resliced.nii', 'overwrite');

powersRatingResliced = resample_space(fmri_data('PowersUnpublished_rating.nii'), mni_templ)
orthviews(powersRatingResliced)
write(powersRatingResliced, 'fname', 'PowersUnpublished_rating_resliced.nii', 'overwrite');

steinfurthRatingResliced = resample_space(fmri_data('Steinfurth2013_rating.nii'), mni_templ)
orthviews(steinfurthRatingResliced)
write(steinfurthRatingResliced, 'fname', 'Steinfurth2013_rating_resliced.nii', 'overwrite');

netaamyResliced = resample_space(fmri_data('NetaUnpublished_amy.nii'), mni_templ)
orthviews(netaamyResliced)
write(netaamyResliced, 'fname', 'NetaUnpublished_amy_resliced.nii', 'overwrite');

powersamyResliced = resample_space(fmri_data('PowersUnpublished_amy.nii'), mni_templ)
orthviews(powersamyResliced)
write(powersamyResliced, 'fname', 'PowersUnpublished_amy_resliced.nii', 'overwrite');

steinfurthamyResliced = resample_space(fmri_data('Steinfurth2013_amy.nii'), mni_templ)
orthviews(steinfurthamyResliced)
write(steinfurthamyResliced, 'fname', 'Steinfurth2013_amy_resliced.nii', 'overwrite');



netaamyResliced = resample_space(fmri_data('Steward2021_quest_ERQ_2.nii'), mni_templ)
orthviews(netaamyResliced)
write(netaamyResliced, 'fname', 'Steward2021_quest_ERQ_resliced.nii', 'overwrite');

netaamyResliced = resample_space(fmri_data('Steward2021_amy_2.nii'), mni_templ)
orthviews(netaamyResliced)
netaamyResliced.dat = netaamyResliced.dat*(-1); % correction coding error (amygdala was positively correlated with itself in the reverse contrast)
orthviews(netaamyResliced)
write(netaamyResliced, 'fname', 'Steward2021_amy_resliced.nii', 'overwrite');

netaamyResliced = resample_space(fmri_data('Steward2021_rating_2.nii'), mni_templ)
orthviews(netaamyResliced)
write(netaamyResliced, 'fname', 'Steward2021_rating_resliced.nii', 'overwrite');



netaamyResliced = resample_space(fmri_data('MinUnpublished_quest_ERQ-SE.nii'), mni_templ)
orthviews(netaamyResliced)
write(netaamyResliced, 'fname', 'MinUnpublished_quest_ERQ-SE_resliced.nii', 'overwrite');






%data_obj = fmri_data('C:\Users\mauri\OneDrive\Arbeit\Projekte\NIDER\NIDER_project\data\tMaps_raw\tstat1.nii')
%testname = 'C:\Users\mauri\OneDrive\Arbeit\Projekte\NIDER\NIDER_project\data\tMaps_raw\tstat1.nii';
%scn_map_image(testname, mni_templ, 'write', strcat(testname, '_resliced.nii'));
