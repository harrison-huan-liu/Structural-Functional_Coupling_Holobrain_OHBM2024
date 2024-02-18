function task_based_T_test()
current_file = mfilename('fullpath');
[current_path, ~, ~] = fileparts(current_file);
for task_i = 1:7
    allcfcpath = [current_path, '\data\to_DFYANG\all_cfc_', num2str(task_i), '.mat'];
    load(allcfcpath)
    % labelnumberpath = [current_path, '\data\to_DFYANG\label_number_', num2str(task_i), '.mat'];
    % load(labelnumberpath)
    all_cfc{1,1}/label_num(1,1);
end
[p(m+1,n), ~, ~]=permutationTest(sample1(:,j)', sample2(:,j)', 1000);
end