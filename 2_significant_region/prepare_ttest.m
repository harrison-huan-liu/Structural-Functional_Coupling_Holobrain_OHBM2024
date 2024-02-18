function [transform_ttest, ave_cfc, count_all, PID, Task_name] = prepare_ttest(current_path, task_i, labeled, with_wavelet, scan, wavelet_general, filter_bandnumber)
% test: 0; retest: 1. (scan)
diag_only = 1;
calculate_cfc = 1;
with_wavelet_char = {'bold','wavelet_default','wavelet_general'};
scan_char = {'test','retest'};

if with_wavelet==0
    wavelet_general=0;
end

if wavelet_general==0
    filter_bandnumber=3;
end

fprintf('prepare task: %d\n', task_i)
sample_i = 1;
BOLD_signal_filename = ['HCP_task_data_', scan_char{scan+1}, '_labeled_', num2str(task_i), '.mat'];
BOLD_signal_filepath = fullfile(current_path, 'HCP_YA', BOLD_signal_filename);

BOLD_signal = load(BOLD_signal_filepath);
subtask_label = unique(BOLD_signal.HCP_task_data(:,2));
ave_cfc = cell(size(subtask_label, 1),1);
all_cfc = cell(2,2);
count_all = zeros(1, size(subtask_label, 1));

% transform_ttest = cell(1080, 1080);
% for crossmat_i = 1:1080
%     for crossmat_j = 1:crossmat_i
%         transform_ttest{crossmat_i,crossmat_j} = zeros(size(BOLD_signal.HCP_task_data, 1), 2);
%     end
% end

PID = cell(2, 1); % size(BOLD_signal.HCP_task_data, 1)
Task_name = cell(2, 1); % size(BOLD_signal.HCP_task_data, 1)

subdict = {};
for m = 1:size(subtask_label, 1)
    if with_wavelet
        if wavelet_general % 1: identify the number of wavelet; 0: 3 bands default
            filename = ['wavelet_', num2str(filter_bandnumber), '.mat'];
            filepath = fullfile(current_path, 'output', 'subdict', filename);
            a = load(filepath);
            structname = fieldnames(a);
            subdict = getfield(a, structname{1});
            if diag_only
                for w_i = 1:size(subdict,2)
                    for w_j = 1:size(subdict,2)
                        ave_cfc{m,w_i+(w_j-1)*filter_bandnumber} = zeros(360,1);
                    end
                end
            else
                ave_cfc{m,1} = zeros(360*size(subdict,2),360*size(subdict,2));
            end
        else
            for k = 1:3
                filename = ['wavelet', num2str(k), '.mat'];
                filepath = fullfile(current_path, 'output', 'subdict', filename);
                a = load(filepath);
                structname = fieldnames(a);
                wavelet = getfield(a, structname{1});
                subdict = [subdict,wavelet];
            end
            ave_cfc{m,1} = zeros(1080,1080);
        end
    else
        ave_cfc{m,1} = zeros(360,360);
    end
end

% transform_ttest = zeros([size(ave_cfc{1,1}),size(BOLD_signal.HCP_task_data, 1)]);
cfc = zeros('like', ave_cfc{1,1});
diag_cfc = cell(1, size(subdict,2));

if filter_bandnumber>10
    select_sample_number = 1;
    update_step = 1;
else
    select_sample_number = size(BOLD_signal.HCP_task_data, 1);
    update_step = 10;
end

upd = textprogressbar(select_sample_number, 'barlength', 20, ...
                                                        'updatestep', update_step, ...
                                                        'startmsg', 'Waiting... ',...
                                                        'endmsg', ' Calculate end! ', ...
                                                        'showbar', true, ...
                                                        'showremtime', true, ...
                                                        'showactualnum', true, ...
                                                        'barsymbol', '+', ...
                                                        'emptybarsymbol', '-');

transform_ttest = [];
count_i = 1;
for i = 1:select_sample_number
    if ~isempty(BOLD_signal.HCP_task_data{i,1})
        if calculate_cfc
            power = [];
            if with_wavelet
                if wavelet_general % 1: identify the number of wavelet; 0: 3 bands default
                    for j = 1:size(subdict,2)
                        if diag_only
                            for k = 1:size(subdict,2)
                                temp_power_j = subdict{1,j}*(BOLD_signal.HCP_task_data{i,1}');
                                temp_power_k = subdict{1,k}*(BOLD_signal.HCP_task_data{i,1}');
                                temp_diag_cfc = corr(temp_power_j', temp_power_k', 'type', 'pearson');
                                diag_cfc{1, j+(k-1)*filter_bandnumber} = diag(temp_diag_cfc);
                            end
                        else
                            power = [power; subdict{1,j}*BOLD_signal.HCP_task_data{i,1}'];
                        end
                    end
                else
                    for k = 1:3
                        power = [power; subdict{1,k}*BOLD_signal.HCP_task_data{i,1}'];
                    end
                end
            else
                power = BOLD_signal.HCP_task_data{i,1}';
            end
            if ~diag_only
                cfc = corr(power', 'type', 'pearson');
            end
            % cfc = power*power';

            % if strcmp(BOLD_signal.HCP_task_data{i,2}, subtask_label{1,1})
            %     transform_ttest(:,:,sample_i) = cfc;
            %     sample_i=sample_i+1;
            % end
            
            % transform_ttest(:,:,i) = cfc;
            % for crossmat_i = 1:1080
            %     for crossmat_j = 1:crossmat_i
            %         transform_ttest{crossmat_i,crossmat_j}(i,task_i) = cfc(crossmat_i,crossmat_j);
            %     end
            % end

            for j = 1:size(subtask_label, 1)
                if strcmp(BOLD_signal.HCP_task_data{i,2}, subtask_label{j,1})
                    if diag_only
                        count_all(1,j) = count_all(1,j)+1;
                        for w_i = 1:size(subdict,2)*size(subdict,2)
                            % all_cfc{count_all(1,j),j} = cfc;
                            ave_cfc{j,w_i} = ave_cfc{j,w_i} + diag_cfc{1,w_i};
                        end
                    else
                        count_all(1,j) = count_all(1,j)+1;
                        % all_cfc{count_all(1,j),j} = cfc;
                        ave_cfc{j,1} = ave_cfc{j,1} + cfc;
                    end
                end
            end
        end

        if labeled
            PID{count_i,1} = BOLD_signal.HCP_task_data{i,3};
            Task_name{count_i,1} = BOLD_signal.HCP_task_data{i,2};
        end
        count_i = count_i + 1;
    end
    upd(i);
end

if calculate_cfc
    for len = 1:size(subtask_label, 1)
        if diag_only
            for w_i = 1:size(subdict,2)*size(subdict,2)
                ave_cfc{len,w_i} = ave_cfc{len,w_i} / count_all(1,len);
            end
        else
            ave_cfc{len,1} = ave_cfc{len,1} / count_all(1,len);
        end
    end
    % all_cfc_select = all_cfc(1:10,1:size(subtask_label, 1));

    % all_cfc_savepath = [current_path, '/data/to_DFYANG/all_cfc_', with_wavelet_char{with_wavelet+wavelet_general+1}, '_', scan_char{scan+1}, '_', num2str(task_i), '.mat'];
    ave_cfc_savename = ['ave_cfc_', with_wavelet_char{with_wavelet+wavelet_general+1}, '_', scan_char{scan+1}, '_', num2str(task_i), '_', num2str(filter_bandnumber), '_', num2str(select_sample_number), '.mat'];
    ave_cfc_savepath = fullfile(current_path, 'output', 'cfc_matrix', ave_cfc_savename);
    % save(all_cfc_savepath, 'all_cfc_select')
    save(ave_cfc_savepath, 'ave_cfc', '-v7.3')
end
end