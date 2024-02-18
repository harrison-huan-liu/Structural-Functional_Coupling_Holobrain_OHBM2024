function [] = rearrange_cfc(task_index,scan,CFC_form)
current_file = mfilename('fullpath');
[current_path, ~, ~] = fileparts(current_file);
[current_path, ~, ~] = fileparts(current_path);

scan_char = {'test','retest'};
Bold_fileName = ['HCP_task_data_', scan_char{scan+1}, '_labeled_', num2str(task_index), '.mat'];
Bold_filePath = fullfile(current_path, 'HCP_YA', Bold_fileName);
Wavelet_filePath = fullfile(current_path, 'output', 'subdict', 'wavelet_360_single_original_pn.mat');
load(Bold_filePath)
load(Wavelet_filePath)

% node_select = [60,120];
node_select = 1:360;
CFC_node = cell(size(node_select,2),1);

subtask_label = unique(HCP_task_data(:,2));
power = cell(1,size(HCP_task_data,1)/size(subtask_label, 1));

wavelet_num = size(subdict,2);
for i = 1:size(node_select,2)
    CFC_node{i,1} = zeros(wavelet_num,wavelet_num);
end

for node_i = 1:size(node_select,2)
    for j = 1:(size(HCP_task_data,1)/size(subtask_label, 1))
        power{1,j} = zeros(wavelet_num,size(HCP_task_data{1,1},1));
    end
    fprintf("The node: %d!\n", node_select(node_i))
    upd = textprogressbar(wavelet_num+(size(HCP_task_data,1)/size(subtask_label, 1)), 'barlength', 10, ...
                         'updatestep', 10, ...
                         'startmsg', 'Waiting... ',...
                         'endmsg', ' Yay!', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '+', ...
                         'emptybarsymbol', '-');
    ave_cfc = zeros(wavelet_num,wavelet_num);
    for wavelet_i = 1:wavelet_num
        count_i = 0;
        wavelet = subdict{1,wavelet_i};
        for sample_i = 1:size(subtask_label, 1):size(HCP_task_data,1)
            count_i = count_i + 1;
            temp_power = wavelet * (HCP_task_data{sample_i, 1}');
            power{1,count_i}(wavelet_i,:) = temp_power(node_select(node_i),:);
        end
        upd(wavelet_i)
    end
    for i = 1:(size(HCP_task_data,1)/size(subtask_label, 1))
        if CFC_form
            temp_corr = power{1,i}*power{1,i}'; % corr(power{1,i}', 'type', 'pearson');
        else
            temp_corr = corr(power{1,i}', 'type', 'pearson');
        end
        ave_cfc = ave_cfc + temp_corr;
        upd(wavelet_num+i)
    end
    CFC_node{node_i} = ave_cfc/(size(HCP_task_data,1)/size(subtask_label, 1));
end

cfc_savename = ['CFC_node_', scan_char{scan+1}, '_', num2str(task_index), '_single_', num2str(CFC_form), '_', num2str(wavelet_num), '.mat'];
cfc_savepath = fullfile(current_path, 'output', 'cfc_matrix', cfc_savename);
save(cfc_savepath, 'CFC_node', '-v7.3')
end


