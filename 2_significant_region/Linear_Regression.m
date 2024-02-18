clc
clear
close all
Prepare_Linear_Regression(7,7,1,1,0,1,360)
Prepare_Linear_Regression(7,7,1,1,1,1,360)

function Prepare_Linear_Regression(task_start,task_end,GLM_difference,with_wavelet, scan, wavelet_general, filter_bandnumber)
    labeled = 1;
    save_sampe_information = 0;
    current_file = mfilename('fullpath');
    [current_path, ~, ~] = fileparts(current_file);
    [current_path, ~, ~] = fileparts(current_path);
    scan_char = {'test','retest'};

    label_information_path = fullfile(current_path, 'output', 'unrestricted_HCP_YA.csv');
    label_information = readcell(label_information_path); % , 1, 'A1:VJ1207'

    for task_i = task_start:task_end % task_start = 2; task_end = 7;
        [~,~,~,PID_i,Task_name_i] = prepare_ttest(current_path,task_i,labeled,with_wavelet, scan, wavelet_general, filter_bandnumber); % transform_ttest_i
        for task_j = 1:task_i-1
            fprintf('Task_data: %d and %d\n', task_i, task_j)
            [~,~,~,PID_j,Task_name_j] = prepare_ttest(current_path,task_j,labeled,with_wavelet, scan, wavelet_general, filter_bandnumber); % transform_ttest_j
            % linearregression_pvalue = zeros(size(transform_ttest_i,1),size(transform_ttest_i,1));% {task_i,task_j}
            % linearregression_estimate = zeros(size(transform_ttest_i,1),size(transform_ttest_i,1));
            
            % temp_i_path = [current_path, '/data/to_DFYANG/transform_ttest_', num2str(task_i), '.mat'];
            % temp_j_path = [current_path, '/data/to_DFYANG/transform_ttest_', num2str(task_j), '.mat'];
            % temp_i = load(temp_i_path)
            % temp_j = load(temp_j_path)
            % transform_ttest_i = temp_i.transform_ttest
            % transform_ttest_j = temp_j.transform_ttest


            %% store the age, gender, and label information
            if save_sampe_information == 1
                if GLM_difference == 0 % not completed
                    [transform_ttest_i, transform_ttest_j, PID_com] = Select_same_PID(PID_i, PID_j, transform_ttest_i, transform_ttest_j);
                    % PID_com = PID_i;
                    [Age_com, Gender_com] = Align_Confunder(label_information, PID_com);
                else
                    LABEL_i = zeros(size(PID_i,1), 1);
                    LABEL_j = ones(size(PID_j,1), 1);
                    LABEL = [LABEL_i;LABEL_j];
                    PID_com = [PID_i;PID_j];
                    [Age_i, Gender_i] = Align_Confunder(label_information, PID_i);
                    [Age_j, Gender_j] = Align_Confunder(label_information, PID_j);
                    Age_com = [Age_i;Age_j];
                    Gender_com = [Gender_i;Gender_j];
                    
                    Age_com_savename = ['Age_com_', scan_char{scan+1}, '_', num2str(task_i), '_', num2str(task_j), '.mat'];
                    Age_com_savepath = fullfile(current_path, 'output', 'Age_com', Age_com_savename);
                    Gender_com_savename = ['Gender_com_', scan_char{scan+1}, '_', num2str(task_i), '_', num2str(task_j), '.mat'];
                    Gender_com_savepath = fullfile(current_path, 'output', 'Gender_com', Gender_com_savename);
                    LABEL_savename = ['LABEL_', scan_char{scan+1}, '_', num2str(task_i), '_', num2str(task_j), '.mat'];
                    LABEL_savepath = fullfile(current_path, 'output', 'LABEL',LABEL_savename);
                    save(Age_com_savepath, 'Age_com')
                    save(Gender_com_savepath, 'Gender_com')
                    save(LABEL_savepath, 'LABEL')
                end
            end


            % for crossmat_i = 1:size(transform_ttest_i,1)
            %     if mod(crossmat_i, 50)==0
            %         fprintf('start linear regression: %d\n', crossmat_i)
            %     end
            %     for crossmat_j = 1:crossmat_i
            %         task1_point = squeeze(transform_ttest_i(crossmat_i,crossmat_j,:));
            %         task2_point = squeeze(transform_ttest_j(crossmat_i,crossmat_j,:));
                    
            %         if all(eig(task1_point)>0) && all(eig(task2_point)>0)
            %             if GLM_difference == 0
            %                 com_length = min(size(transform_ttest_i,3),size(transform_ttest_j,3));
                            
            %                 X = [Age_com(1:com_length,1),Gender_com(1:com_length,1),task1_point];
                            
            %                 % mdl = fitlm(X,task2_point);
            %                 % linearregression_pvalue(crossmat_i,crossmat_j) = mdl.Coefficients{4,4};
            %                 % linearregression_estimate(crossmat_i,crossmat_j) = mdl.Coefficients{4,1};

            %                 [~,~,~,~,mdl] = regress(task2_point, X);
            %                 linearregression_pvalue(crossmat_i,crossmat_j) = mdl{3,3};

            %                 % if mdl.Coefficients{4,4}<10^-20

            %                 % [~,ttest_value(crossmat_i,crossmat_j)] = ttest(task1_point, task2_point);
            %                 % if ttest_value(crossmat_i,crossmat_j) < 0.000001
            %                 % [ttest_value(crossmat_i,crossmat_j),~,~] = permutationTest(task1_point, task2_point, 1000000);
            %                 % fprintf('p-value: %d\n', ttest_value(crossmat_i,crossmat_j))
            %                 % {task_i,task_j}
            %             else
            %                 CFC_point = [task1_point;task2_point];
                        
            %                 X = [ones(size(Age_com,1),1),Age_com,Gender_com,CFC_point];

            %                 % mdl = fitlm(X,LABEL);
            %                 % linearregression_pvalue(crossmat_i,crossmat_j) = mdl.Coefficients{4,4};
            %                 % linearregression_estimate(crossmat_i,crossmat_j) = mdl.Coefficients{4,1};

            %                 % [~,~,~,~,mdl] = regress(LABEL, X);
            %                 % linearregression_pvalue(crossmat_i,crossmat_j) = mdl{3,3};
            %             end
            %         end
            %     end
            % end
            % % mdl.Coefficients
            % if GLM_difference == 0
            %     linearregression_pvalue_savepath = [current_path, '/data/to_DFYANG/linearregression_pvalue_', num2str(task_i), '_', num2str(task_j), '.mat'];
            %     linearregression_estimate_savepath = [current_path, '/data/to_DFYANG/linearregression_estimate_', num2str(task_i), '_', num2str(task_j), '.mat'];
            % else if with_wavelet == 0
            %     linearregression_pvalue_savepath = [current_path, '/data/to_DFYANG/linearregression_pvalue_bold_', num2str(task_i), '_', num2str(task_j), '.mat'];
            %     linearregression_estimate_savepath = [current_path, '/data/to_DFYANG/linearregression_estimate_bold_', num2str(task_i), '_', num2str(task_j), '.mat'];
            % else
            %     linearregression_pvalue_savepath = [current_path, '/data/to_DFYANG/linearregression_pvalue_ondifference_', num2str(task_i), '_', num2str(task_j), '.mat'];
            %     linearregression_estimate_savepath = [current_path, '/data/to_DFYANG/linearregression_estimate_ondifference_', num2str(task_i), '_', num2str(task_j), '.mat'];
            % end
            % save(linearregression_pvalue_savepath, 'linearregression_pvalue')
            % save(linearregression_estimate_savepath, 'linearregression_estimate')
        end
    end
end