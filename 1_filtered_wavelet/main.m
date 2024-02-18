current_file = mfilename('fullpath');
[current_path, ~, ~] = fileparts(current_file);
[current_path, ~, ~] = fileparts(current_path);

gsp_path = fullfile(current_path,'1_filtered_wavelet','gspbox');
utils_path = fullfile(current_path,'1_filtered_wavelet','utils');
addpath(genpath(gsp_path),genpath(utils_path))

gsp_start
gsp_install

cd ..
current_path = Filter_identify_general(4,2,0,3);
% current_path = Filter_identify_general(40,0,0,3);

% [to_DFYANG_dirname] = HCP_Bold_collect(current_path,1); % test: 0; retest: 1.
