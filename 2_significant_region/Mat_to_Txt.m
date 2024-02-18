function Mat_to_Txt(task_start,task_end)
    current_file = mfilename('fullpath');
    [current_path, ~, ~] = fileparts(current_file);
    [current_path, ~, ~] = fileparts(current_path);

    % atlaspath = fullfile(current_path, '\data\HCP_MMP_360.xlsx');
    % atlas = readcell(atlaspath);
    % cortex={'Primary_Visual','MT+_Complex_and_Neighboring_Visual_Areas','Dorsal_Stream_Visual','Early_Visual','Ventral_Stream_Visual','Somatosensory_and_Motor','Premotor','Posterior_Cingulate','Early_Auditory','Temporo-Parieto_Occipital_Junction','Dorsolateral_Prefrontal','Temporo-Parieto-Occipital_Junction','Superior_Parietal','Paracentral_Lobular_and_Mid_Cingulate','Anterior_Cingulate_and_Medial_Prefrontal','Orbital_and_Polar_Frontal','Inferior_Frontal','Posterior_Opercular','Insular_and_Frontal_Opercular','Auditory_Association','Inferior_Parietal','Medial_Temporal','Lateral_Temporal'}';

    for task_i = task_start:task_end % task_start = 2; task_end = 7;
        for task_j = 1:task_i-1
            filename = [current_path, '/data/to_DFYANG/linearregression_pvalue_', num2str(task_i), '_', num2str(task_j), '.mat'];
            % filename = [current_path, '/data/to_DFYANG/ttest_value_', num2str(task_i), '_', num2str(task_j), '.mat'];
            load(filename);

            % transform_cfc = cell(2,1);
            % transform_cfc{1} = zeros(1080,1080);
            % transform_cfc{2} = zeros(1080,1080);
            % order = 1;
            % kkkk = zeros(23,1);
            % for i = 1:23
            %     for j = 2:361
            %         if strcmp(cortex{i},atlas{j,7})
            %             for k = 1:2
            %                 transform_cfc{k}(:,order) = cfc{k}(:,j-1);
            %                 transform_cfc{k}(order,:) = cfc{k}(j-1,:);
            %                 transform_cfc{k}(:,360+order) = cfc{k}(:,360+j-1);
            %                 transform_cfc{k}(360+order,:) = cfc{k}(360+j-1,:);
            %                 transform_cfc{k}(:,720+order) = cfc{k}(:,720+j-1);
            %                 transform_cfc{k}(720+order,:) = cfc{k}(720+j-1,:);
            %             end
            %             order = order + 1;
            %         end
            %     end
            %     kkkk(i) = order-1;
            % end

            for k=1:3
                for kk=1:3
                    temp = linearregression_pvalue(360*k-359:360*k,360*kk-359:360*kk);
                    % temp = ttest_value(360*k-359:360*k,360*kk-359:360*kk);
                    for row = 1:size(temp,1)
                        for col = 1:row
                            if temp(row,col)<10e-5
                                temp(row,col)=1;
                                temp(col,row)=1;
                            else
                                temp(row,col)=0;
                                temp(col,row)=0;
                            end
                        end
                    end
                    aaa = sum(temp,'all');
                    fprintf('%d\n', aaa)
                    % savepath = [current_path, '/data/to_DFYANG/ttest_value_', num2str(task_i), '_', num2str(task_j), '_', num2str(k), '_', num2str(kk), '.txt'];
                    savepath = [current_path, '/data/to_DFYANG/linearregression_pvalue_', num2str(task_i), '_', num2str(task_j), '_', num2str(k), '_', num2str(kk), '.txt'];
                    % save savepath -ascii temp
                    save(savepath, 'temp', '-ascii')
                end
            end
        end
    end

end