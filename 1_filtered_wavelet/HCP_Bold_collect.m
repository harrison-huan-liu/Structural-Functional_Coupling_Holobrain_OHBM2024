function [to_DFYANG_dirname] = HCP_Bold_collect(current_path,scan)
    HCP_datafolder = fullfile(current_path, 'HCP_YA');
    to_DFYANG_filename = dir(HCP_datafolder);
    scan_char = {'test','retest'};
        
    to_DFYANG_dirname = strings(7,1);
    strategy_i = 1;
    for i = 3:size(to_DFYANG_filename, 1)
        if to_DFYANG_filename(i).isdir == 1
            to_DFYANG_dirname(strategy_i, 1) = to_DFYANG_filename(i).name;
            strategy_i = strategy_i + 1;
        end
    end

    HCP_task_datafolder = strings(size(to_DFYANG_dirname, 1),1);
    for i = 1:size(to_DFYANG_dirname, 1)
        HCP_task_data = cell(2, 2); % cell(size(HCP_task_filename, 1), 2)
        save_sample_hcp_task_i = 1;
        fprintf('Task_folder: %d %s\n', i, to_DFYANG_dirname(i,1))
        HCP_task_datafolder(i, 1) = fullfile(HCP_datafolder, to_DFYANG_dirname(i,1));
        HCP_task_file = dir(fullfile(HCP_task_datafolder(i, 1),'*.csv'));
        HCP_task_labelfile = dir(fullfile(HCP_task_datafolder(i, 1),'*.xlsx'));
        HCP_task_labelfilepath = fullfile(HCP_task_datafolder(i, 1), HCP_task_labelfile(1).name);
        HCP_task_filename = {HCP_task_file.name}';
        TASK_LABEL = [];
        TASK_LABEL = readtable(HCP_task_labelfilepath, 'ReadVariableNames', true, 'VariableNamingRule', 'preserve'); % , 'VariableNamesRange', '1:1'

        defaultsamplePid = '000000';
        upd = textprogressbar(size(HCP_task_filename,1) , defaultsamplePid, 'barlength', 20, ...
                                                        'updatestep', 10, ...
                                                        'startmsg', 'Waiting... ',...
                                                        'endmsg', ' Yay!', ...
                                                        'showbar', true, ...
                                                        'showremtime', true, ...
                                                        'showactualnum', true, ...
                                                        'barsymbol', '+', ...
                                                        'emptybarsymbol', '-');

        for sample_hcp_task_i = 1:size(HCP_task_filename,1)
            temp_hcp_task_filename = fullfile(HCP_task_datafolder(i, 1), HCP_task_filename{sample_hcp_task_i});
            HCP_task_bold = [];
            HCP_task_bold = readcell(temp_hcp_task_filename); % [~,~,HCP_task_bold] = csvread(temp_hcp_task_filename);
            temp_hcp_bold = [];
            for time_i = 1:size(HCP_task_bold,1)
                % temp = strsplit(HCP_task_bold{time_i});
                for size_i = 1:size(HCP_task_bold,2)
                    temp_hcp_bold(time_i,size_i) = HCP_task_bold{time_i,size_i};% str2num(temp{1,size_i});
                end
            end
            TASK_information = [];
            TASK_information = strsplit(HCP_task_filename{sample_hcp_task_i}, '_');
            if TASK_LABEL{1,3}{1} == '0'
                task_begin = 2;
            else
                task_begin = 1;
            end
            for task_label_i = task_begin:size(TASK_LABEL,1) % readtable can't make the first column as variablename
                TASK_LABEL_ID = [];
                TASK_LABEL_ID = strsplit(TASK_LABEL{task_label_i, 2}{1}, '_');
                if scan == 0
                    scan_order = 'LR';
                else
                    scan_order = 'RL';
                end
                if strcmp(TASK_information{1,1}, TASK_LABEL_ID{1,1}) && strcmp(TASK_LABEL_ID{1,3}, scan_order) && strcmp(TASK_information{1,3}, scan_order) && strcmp(TASK_information{1,2}, TASK_LABEL_ID{1,2})% 'LR': test; 'RL': retest
                    LABEL_NAME = [];
                    LABEL_NAME = unique(table2cell(TASK_LABEL(task_begin,3:end))); % unique(TASK_LABEL{task_begin,3:end})
                    for label_name_i = 1:size(LABEL_NAME, 2)
                        if ~strcmp(LABEL_NAME{1,label_name_i}, 'rest')
                            select_timepoint = 1;
                            save_hcp_bold = [];
                            for timepoint = 1:size(temp_hcp_bold, 1)
                                if strcmp(TASK_LABEL{task_label_i, timepoint + 2}{1}, LABEL_NAME{1,label_name_i})
                                    save_hcp_bold(select_timepoint,:) = temp_hcp_bold(timepoint,:);
                                    select_timepoint = select_timepoint + 1;
                                end
                            end
                            % fprintf('%s \n', LABEL_NAME{1,label_name_i})
                            HCP_task_data{save_sample_hcp_task_i,1} = save_hcp_bold;
                            HCP_task_data{save_sample_hcp_task_i,2} = LABEL_NAME{1,label_name_i};
                            HCP_task_data{save_sample_hcp_task_i,3} = TASK_LABEL_ID{1,1};
                            save_sample_hcp_task_i = save_sample_hcp_task_i + 1;
                        end
                    end
                end
            end
            
            upd(sample_hcp_task_i,TASK_information{1,1});
        end
        fprintf('split_sample_number: %d \n', save_sample_hcp_task_i)
        save_filename = ['HCP_task_data_', scan_char{scan+1}, '_labeled_', num2str(i),'.mat'];
        save_filepath = fullfile(HCP_datafolder, save_filename);
        save(save_filepath,'HCP_task_data')
    end
end
