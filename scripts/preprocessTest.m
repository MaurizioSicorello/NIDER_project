% load MNI template
cd('../data/masks')
mni_templ = fmri_data('MNI152_T1_2mm_brain_mask.nii');


%load images
cd('../tMaps_raw')

image_names = filenames(fullfile(pwd, '*.nii'), 'absolute')

for i=1:length(image_names)
    
    [~,name,~] = fileparts(image_names{i});
    scn_map_image(image_names{i}, mni_templ, 'write', strcat(name, '_resliced.nii'));

end

image_names_resliced = filenames(fullfile(pwd, '*_resliced.nii'), 'absolute')
%%%% write code to automatically move images!!



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
