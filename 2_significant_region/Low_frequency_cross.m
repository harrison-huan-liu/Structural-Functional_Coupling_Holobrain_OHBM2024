current_file = mfilename('fullpath');
[current_path, ~, ~] = fileparts(current_file);
[current_path, ~, ~] = fileparts(current_path);

scan_char = {'test','retest'};
for scan = 0:1
    for task_index = 1:7
        fprintf('Scan: %s ;Task: %d !\n', scan_char{scan+1}, task_index)
        count_i = 1;
        task_filename = ['HCP_task_data_', scan_char{scan+1}, '_labeled_', num2str(task_index), '.mat'];
        task_path = fullfile(current_path, 'HCP_YA', task_filename);
        HCP_task_data = load(task_path);
        Wavelet_filePath = fullfile(current_path, 'output', 'subdict', 'wavelet_360_single_original_pn.mat');
        load(Wavelet_filePath)
        sample_num = size(HCP_task_data.HCP_task_data,1);
        subtask_label = unique(HCP_task_data.HCP_task_data(:,2));
        for eigen_index = 1:10
            index_i = 1;
            upd = textprogressbar(10*360*sample_num/size(subtask_label, 1), 'barlength', 40, ...
                            'updatestep', 10, ...
                            'startmsg', 'Waiting... ',...
                            'endmsg', ' Yay!', ...
                            'showbar', true, ...
                            'showremtime', true, ...
                            'showactualnum', true, ...
                            'barsymbol', '+', ...
                            'emptybarsymbol', '-');
            for node=1:360
                for i=1:size(subtask_label, 1):sample_num
                    upd(count_i)
                    for j=1:size(subdict,2)
                        bold_vector = subdict{1,eigen_index}(node,:)*HCP_task_data.HCP_task_data{i,1}';
                        interference = subdict{1,j}(node,:)*HCP_task_data.HCP_task_data{i,1}';
                        Boxplot(index_i,j) = bold_vector*interference';
                    end
                    index_i = index_i + 1;
                    count_i = count_i + 1;
                end
            end

            wavelet_num = size(subdict,2);
            savefilename = ['Boxplot_', scan_char{scan+1}, '_', num2str(task_index), '_', num2str(wavelet_num), '_', num2str(eigen_index), '.mat'];
            savepath = fullfile(current_path, 'output', 'cfc_matrix', savefilename);
            save(savepath, 'Boxplot', '-v7.3')
        end
    end
end