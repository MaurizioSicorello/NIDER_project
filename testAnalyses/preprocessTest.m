% load MNI template
mni_templ = fmri_data('MNI152_T1_2mm_brain_mask.nii');

cd('./testDataFear')

% old code. unpacks gz files
% image_names = filenames(fullfile(pwd, '*.gz'), 'absolute')
% data_obj = fmri_data(image_names);

image_names = filenames(fullfile(pwd, '*.nii'), 'absolute')

for i=1:length(image_names)
    
    [~,name,~] = fileparts(image_names{i});
    scn_map_image(image_names{i}, mni_templ, 'write', strcat(name, '_resliced.nii'));

end

image_names_resliced = filenames(fullfile(pwd, '*_resliced.nii'), 'absolute')


plot(data_obj)
descriptives(data_obj)