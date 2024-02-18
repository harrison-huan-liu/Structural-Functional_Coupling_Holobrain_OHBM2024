function Save_task_data()
    current_file = mfilename('fullpath');
    [current_path, ~, ~] = fileparts(current_file);
    [current_path, ~, ~] = fileparts(current_path);
    labeled = 1;
    for task = 1:7
        fprintf('Task_data: %d \n', task)
        [transform_ttest,~,~,PID,Task_name] = prepare_ttest(current_path,task,labeled);
        transform_ttest_savepath = [current_path, '/data/to_DFYANG/transform_ttest_', num2str(task), '.mat'];
        save(transform_ttest_savepath, 'transform_ttest', '-v7.3')
        PID_savepath = [current_path, '/data/to_DFYANG/PID_', num2str(task), '.mat'];
        save(PID_savepath, 'PID')
        Task_name_savepath = [current_path, '/data/to_DFYANG/Task_name_', num2str(task), '.mat'];
        save(Task_name_savepath, 'Task_name')
    end
end