function Generate_CFC(upsampling, task_start, task_end, with_wavelet)
% upsignalFile = dir(fullfile('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\HCP_dynamics\upsampled_signal\','*.mat'));
current_file = mfilename('fullpath');
[current_path, ~, ~] = fileparts(current_file);
[current_path, ~, ~] = fileparts(current_path);
ave_cfc = cell(3,2);
all_power = cell(2,1);
all_power{1} = [];
all_power{2} = [];
labeled = 0;

if upsampling
    upsignaldatafilepath = fullfile(current_path, '/multi-boundary_ODE/upsampled_signal/','*.mat');
    upsignalFile = dir(upsignaldatafilepath);
    upsignalName = {upsignalFile.name}';
    for k = 1:size(upsignalName, 1)
    %     filename = ['D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\HCP_dynamics\upsampled_signal\',upsignalName{k}];
        filename = [current_path, '/multi-boundary_ODE/upsampled_signal/',upsignalName{k}];
        load(filename)
        ave_cfc{floor((k-1)/2)+1,mod(k-1,2)+1} = zeros(360,360);
        filename = [current_path, '/data/to_DFYANG/wavelet', num2str(floor((k-1)/2)+1),'.mat'];
        a = load(filename);
        structname = fieldnames(a);
        subdict = getfield(a,structname{1});
        power = subdict*upsampled_signal;
        samp = mod(k-1,2)+1;
        if samp==1
            all_power{1}=[all_power{1};power];
        else
            all_power{2}=[all_power{2};power];
        end
        ave_cfc{floor((k-1)/2)+1,mod(k-1,2)+1} = diag(dot(power', power'));
    %     figure
    %     imagesc(ave_cfc{floor((k-1)/2)+1,mod(k-1,2)+1})
    end
else
    % ttest_value = cell(2,2);
    for task_i = task_start:task_end % task_start = 2; task_end = 7;
        for task_j = 1:task_i-1
            fprintf('Task_data: %d and %d\n', task_i, task_j)
            [transform_ttest_i,~,~,PID_i,~] = prepare_ttest(current_path,task_i,labeled,with_wavelet);
            [transform_ttest_j,~,~,PID_j,~] = prepare_ttest(current_path,task_j,labeled,with_wavelet);
            ttest_value = zeros(size(transform_ttest_i,1),size(transform_ttest_i,1));% {task_i,task_j}

            % [transform_ttest_i, transform_ttest_j, PID_com] = Select_same_PID(PID_i, PID_j, transform_ttest_i, transform_ttest_j)

            for crossmat_i = 1:size(transform_ttest_i,1)
                if mod(crossmat_i, 50)==0
                    fprintf('start ttest: %d\n', crossmat_i)
                end
                for crossmat_j = 1:crossmat_i
                    com_length = min(size(transform_ttest_i,3),size(transform_ttest_j,3));
                    task1_point = squeeze(transform_ttest_i(crossmat_i,crossmat_j,1:com_length));
                    task2_point = squeeze(transform_ttest_j(crossmat_i,crossmat_j,1:com_length));
                    [~,ttest_value(crossmat_i,crossmat_j)] = ttest(task1_point, task2_point);
                    % if ttest_value(crossmat_i,crossmat_j) < 0.000001
                    % [ttest_value(crossmat_i,crossmat_j),~,~] = permutationTest(task1_point, task2_point, 1000000);
                    % fprintf('p-value: %d\n', ttest_value(crossmat_i,crossmat_j))
                    % {task_i,task_j}
                end
            end

            % savepath = [current_path, '\data\to_DFYANG\all_cfc_', num2str(task_i), '.mat'];
            % save(savepath, 'all_cfc', '-v7.3')
            % savepath = [current_path, '\data\to_DFYANG\ave_cfc_', num2str(task_i), '.mat'];
            % save(savepath, 'ave_cfc')
            % savenumberpath = [current_path, '\data\to_DFYANG\label_number_', num2str(task_i), '.mat'];
            % save(savenumberpath, 'label_number')
            ttest_savepath = [current_path, '/data/to_DFYANG/ttest_value_', num2str(task_i), '_', num2str(task_j), '.mat'];
            save(ttest_savepath, 'ttest_value')
        end
    end
    % for j_1 = 2:7
    %     for j_2 = 1:j_1-1
    %         fprintf('ttest_task: %d and %d\n', j_1, j_2)
    %         ttest_value{j_1,j_2} = zeros(1080,1080);
    %         for crossmat_i = 1:1080
    %             for crossmat_j = 1:crossmat_i
    %                 [ttest_value{j_1,j_2}(crossmat_i,crossmat_j),~,~] = permutationTest(transform_ttest{crossmat_i,crossmat_j}(:,j1), transform_ttest{crossmat_i,crossmat_j}(:,j2), 1000);
    %             end
    %         end
    %     end
    % end
end


% cfc = cell(2,1);
% for i=1:2
%     cfc{i}=cov(all_power{i}');
%     figure
%     imagesc(cfc{i})
% end

% atlaspath = fullfile(current_path, '\data\HCP_MMP_360.xlsx');
% [~,~,atlas] = xlsread(atlaspath,1,'A1:M361');
% cortex={'Primary_Visual','MT+_Complex_and_Neighboring_Visual_Areas','Dorsal_Stream_Visual','Early_Visual','Ventral_Stream_Visual','Somatosensory_and_Motor','Premotor','Posterior_Cingulate','Early_Auditory','Temporo-Parieto_Occipital_Junction','Dorsolateral_Prefrontal','Temporo-Parieto-Occipital_Junction','Superior_Parietal','Paracentral_Lobular_and_Mid_Cingulate','Anterior_Cingulate_and_Medial_Prefrontal','Orbital_and_Polar_Frontal','Inferior_Frontal','Posterior_Opercular','Insular_and_Frontal_Opercular','Auditory_Association','Inferior_Parietal','Medial_Temporal','Lateral_Temporal'}';
% cfc = cell(2,1);
% for label_i = 1:2
%     cfc{label_i,1} = ave_cfc{label_i,1}/label_number(label_i,1);
% end
% 
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
% 
% plot_cfc = zeros(1080,1080);
% for i = 1:2
%     figure
%     for j = 1:1080
%         for k = 1:1080
%             if transform_cfc{i}(j,k)<0
%                 plot_cfc(j,k) = -sqrt(-transform_cfc{i}(j,k));
%             else
%                 plot_cfc(j,k) = sqrt(transform_cfc{i}(j,k));
%             end
%         end
%     end
%     imagesc(plot_cfc)
%     axis square  % 让图片变成方形，增加美观度
%     colormap(jet);
% end
% 
% figure
% for k=1:23
%     for m=1:3
%         x=(360*(m-1)+kkkk(k))*ones(1080,1);
%         y=1:1080;
%         plot(x,y)
%         hold on
%         plot(y,x)
%         hold on
%     end
% end
% 
% 
% ave_bold_cfc = zeros(270,270);
% ave_bold_cfc_cn = zeros(270,270);
% ave_bold_cfc_ad = zeros(270,270);
% ave_fMRI_all_cfc = zeros(270,270);
% ave_fMRI_cfc = cell(3,1);
% for fMRI_scale = 1:3
%     ave_fMRI_cfc{fMRI_scale} = zeros(90,90);
% end
% 
% windows = 30;
% step = 15;
% c = load('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Localized_Spectral_Graph_Filter_Frames_Code\aal116_fMRI.mat');
% sample_num = 0;
% cn_bold_num = 0;
% ad_bold_num = 0;
% for sample_bold_i = 1:135
% %     str_fMRI = ['D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\ADNI_txt\',aal_subcode];
% %     temp_fMRI = load(str_fMRI);
%     sample_bold = c.aal116_fMRI(sample_bold_i).BOLD;
%     fMRI_freq_1 = [];
%     fMRI_freq_2 = [];
%     fMRI_freq_3 = [];
%     fMRI_freq = [];
%     for node_i=1:90
%         fMRI_freq_1(:,node_i) = band_filter_signal(sample_bold(:,node_i),0.5,0.03,0.04);
%         fMRI_freq_2(:,node_i) = band_filter_single(sample_bold(:,node_i),0.5,0.03,0.06);
%         fMRI_freq_3(:,node_i) = band_filter_single(sample_bold(:,node_i),0.5,0.06,0.09);
%     end
%     fMRI_freq = [fMRI_freq_1,fMRI_freq_2,fMRI_freq_3];
%     for sliding_time=1:floor(size(fMRI_freq,1)/windows)
%         temp_bold_cfc = cov(fMRI_freq(step*(sliding_time-1)+1:step*(sliding_time-1)+windows,:));
%         ave_bold_cfc = ave_bold_cfc + temp_bold_cfc;
%     end
% 
%     fMRI_all_power = [];
%     for fMRI_scale = 1:3
%         norm_wavelet = normalize(real(subdict_all{sample_bold_i,fMRI_scale}));
%         fMRI_power = norm_wavelet*sample_bold(:,1:90)';
%         fMRI_all_power = [fMRI_all_power;fMRI_power];
%     end
%     for sliding_time=1:floor(size(fMRI_all_power,2)/windows)
%         for fMRI_scale = 1:3
%             fMRI_cov = cov(fMRI_power(:,step*(sliding_time-1)+1:step*(sliding_time-1)+windows)');
%             ave_fMRI_cfc{fMRI_scale} = ave_fMRI_cfc{fMRI_scale} + fMRI_cov;
%         end
%         temp_fMRI_all_cfc = cov(fMRI_all_power(:,step*(sliding_time-1)+1:step*(sliding_time-1)+windows)');
%         ave_fMRI_all_cfc = ave_fMRI_all_cfc + temp_fMRI_all_cfc;
%         sample_num = sample_num + 1;
%     end
%     if strcmp(c.aal116_fMRI(sample_bold_i).Label,'CN')
%         ave_bold_cfc_cn = ave_bold_cfc_cn + temp_bold_cfc;
%         cn_bold_num = cn_bold_num + 1;
%     elseif strcmp(c.aal116_fMRI(sample_bold_i).Label,'LMCI')||strcmp(c.aal116_fMRI(sample_bold_i).Label,'AD')
%         ave_bold_cfc_ad = ave_bold_cfc_ad + temp_bold_cfc;
%         ad_bold_num = ad_bold_num + 1;
%     end
% end
% 
% % % show the covariance matrix of signals in different frequency
% % for fMRI_scale = 1:3
% %     for i=1:90
% %         ave_fMRI_cfc{fMRI_scale}(i,i)=0;
% %     end
% %     figure
% %     imagesc(ave_fMRI_cfc{fMRI_scale}/sample_num)
% %     set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% % end
% 
% % % normalize by the maximum of each frequency
% % for i=1:3
% %     for j=1:3
% %         ave_bold_cfc((i-1)*90+1:i*90,(j-1)*90+1:j*90) = ave_bold_cfc((i-1)*90+1:i*90,(j-1)*90+1:j*90)/max(max(ave_bold_cfc((i-1)*90+1:i*90,(j-1)*90+1:j*90)));
% %     end
% % end
% % for i=1:270
% % %     ave_fMRI_all_cfc(:,i) = ave_fMRI_all_cfc(:,i)/ave_fMRI_all_cfc(i,i);
% %     ave_fMRI_all_cfc(i,i)=0;
% %     ave_bold_cfc(i,i)=0;
% % end
% 
% ave_fMRI_all_cfc = ave_fMRI_all_cfc/sample_num;
% ave_bold_cfc = ave_bold_cfc/sample_num;
% ave_bold_cfc_cn = ave_bold_cfc_cn/cn_bold_num;
% ave_bold_cfc_ad = ave_bold_cfc_ad/ad_bold_num;
% ave_bold_cfc(ave_bold_cfc>1)=1;
% ave_bold_cfc(ave_bold_cfc<-1)=-1;
% figure
% imagesc(ave_bold_cfc)
% set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% figure
% imagesc(ave_bold_cfc_cn)
% set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% figure
% imagesc(ave_bold_cfc_ad)
% set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% figure
% imagesc(ave_fMRI_all_cfc)
% set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% 
% 
% windows_up = 3000;
% step_up = 1500;
% count_scale_num = zeros(3,1);
% temp_power_freq = cell(135,3);
% for i = 1:size(upsignalName,1)
%     split_upsignalName = strsplit(upsignalName{i,1},'_');
%     split_scale = split_upsignalName{1,2};
%     split_sample = split_upsignalName{1,3};
%     split_sample = strsplit(split_sample,'.');
%     split_sample = split_sample{1,1};
%     for j = 1:3
%         if strcmp(split_scale,num2str(j))
%             temp_filename = ['D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Identified pipeline server\upsampled signal\',upsignalName{i,1}];
%             b = load(temp_filename);
%             temp_signal = real(b.upsampled_signal);
%             norm_wavelet = normalize(real(subdict_all{str2num(split_sample),j}));
%             temp_power = norm_wavelet*temp_signal;
% %             temp_power_freq_1 = [];
% %             temp_power_freq_2 = [];
% %             temp_power_freq_3 = [];
% %             temp_power_freq = [];
%             if j==1
%                 for node_i = 1:90
%                     temp_power_freq{str2num(split_sample),j}(:,node_i) = band_filter_single(temp_power(node_i,:),50,0.1,0.3);
%                 end
%             elseif j==2
%                 for node_i = 1:90
%                     temp_power_freq{str2num(split_sample),j}(:,node_i) = band_filter_single(temp_power(node_i,:),50,0.3,0.6);
%                 end
%             else
%                 for node_i = 1:90
%                     temp_power_freq{str2num(split_sample),j}(:,node_i) = band_filter_single(temp_power(node_i,:),50,0.6,1);
%                 end
%             end
% %             temp_power_freq=[temp_power_freq_1,temp_power_freq_2,temp_power_freq_3];
% 
% %             % identify frequency
% %             color=[1 0 0;0 1 0;0 0 1;0.5 1 1;1 1 0.5;1 0.5 1; 0 0 0.5; 0.5 0 0;0 0.5 0;1 0.5 0.5; 0.5 1 0.5;0.5 0.5 1;1 1 0;0 1 1;1 0 1];
% %             Fs = 50;            % Sampling frequency                    
% %             L = size(b.upsampled_signal(1,:),2);             % Length of signal
% %             
% %             S = temp_power; % Graph_temp(j,1).Harmonics';
% %             Y=fft(S(1,1:L));
% %             P2=abs(Y/L);
% %             P1=P2(1:L/2+1);
% %             P1(2:end-1)=2*P1(2:end-1);
% % 
% %             f = Fs*(0:(L/2))/L;
% %             plot(f,P1,'color',color(j,:))
% %             hold on
% %             title("Amplitude Spectrum of BOLD signals")
% %             xlabel("f (Hz)")
% %             ylabel("|P1(f)|")
% 
% %             % calculate the covariance matrix of signals in different frequency
% %             for sliding_time=1:floor(size(temp_power,2)/windows_up)
% %                 temp_cfc = cov(temp_power(:,step_up*(sliding_time-1)+1:step_up*(sliding_time-1)+windows_up)');
% %                 ave_cfc{j} = ave_cfc{j} + temp_cfc;
% %                 count_scale_num(j) = count_scale_num(j) + 1;
% %             end
%         end
%     end
% end
% 
% % calculate cfc
% count_scale_freq_num = 0;
% % temp_power_all_freq_cfc = zeros(270,270);
% ave_power_all_freq_cfc = zeros(270,270);
% ave_power_all_freq_cfc_cn = zeros(270,270);
% ave_power_all_freq_cfc_ad = zeros(270,270);
% cn_all_freq_num = 0;
% ad_all_freq_num = 0;
% for i=1:135
%     if ~isempty(temp_power_freq{i,1})&&~isempty(temp_power_freq{i,2})&&~isempty(temp_power_freq{i,3})
%         temp_power_all_freq=[temp_power_freq{i,1},temp_power_freq{i,2},temp_power_freq{i,3}];
%         for sliding_time=1:floor(size(temp_power_all_freq,1)/windows_up)
%             temp_power_all_freq_cfc = cov(temp_power_all_freq(step_up*(sliding_time-1)+1:step_up*(sliding_time-1)+windows_up,:));
%             ave_power_all_freq_cfc = ave_power_all_freq_cfc + temp_power_all_freq_cfc;
%             count_scale_freq_num = count_scale_freq_num + 1;
%         end
%         if strcmp(c.aal116_fMRI(i).Label,'CN')
%             ave_power_all_freq_cfc_cn = ave_power_all_freq_cfc_cn + temp_power_all_freq_cfc;
%             cn_all_freq_num = cn_all_freq_num + 1;
%         elseif strcmp(c.aal116_fMRI(i).Label,'LMCI')||strcmp(c.aal116_fMRI(i).Label,'AD')
%             ave_power_all_freq_cfc_ad = ave_power_all_freq_cfc_ad + temp_power_all_freq_cfc;
%             ad_all_freq_num = ad_all_freq_num + 1;
%         end
%     end
% end
% 
% % % show the covariance matrix of signals in different frequency
% % for k = 1:3
% %     ave_cfc{k} = ave_cfc{k}/count_scale_num(k);
% %     figure
% %     for l=1:90
% %         ave_cfc{k}(l,l)=0;
% %     end
% %     imagesc(ave_cfc{k})
% % end
% 
% % % normalize by the maximum of each frequency
% % for i=1:3
% %     for j=1:3
% %         ave_power_all_freq_cfc((i-1)*90+1:i*90,(j-1)*90+1:j*90) = ave_power_all_freq_cfc((i-1)*90+1:i*90,(j-1)*90+1:j*90)/max(max(ave_power_all_freq_cfc((i-1)*90+1:i*90,(j-1)*90+1:j*90)));
% %     end
% % end
% % for i=1:90*3
% % %     ave_power_all_freq_cfc(:,i) = ave_power_all_freq_cfc(:,i)/ave_power_all_freq_cfc(i,i);
% %     ave_power_all_freq_cfc(i,i)=0;
% % end
% 
% ave_power_all_freq_cfc = ave_power_all_freq_cfc/count_scale_freq_num;
% ave_power_all_freq_cfc_cn = ave_power_all_freq_cfc_cn/cn_all_freq_num;
% ave_power_all_freq_cfc_ad = ave_power_all_freq_cfc_ad/ad_all_freq_num;
% ave_power_all_freq_cfc(ave_power_all_freq_cfc>1)=1;
% ave_power_all_freq_cfc(ave_power_all_freq_cfc<-1)=-1;
% figure
% imagesc(ave_power_all_freq_cfc)
% set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% figure
% imagesc(ave_power_all_freq_cfc_cn)
% set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% figure
% imagesc(ave_power_all_freq_cfc_ad)
% set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
% 
% save('ave_fMRI_all_cfc_2.mat','ave_fMRI_all_cfc','-v6')
% save('ave_bold_cfc_2.mat','ave_bold_cfc','-v6')
% save('ave_power_all_freq_cfc_2.mat','ave_power_all_freq_cfc','-v6')
% 
% save('ave_bold_cfc_cn.mat','ave_bold_cfc_cn','-v6')
% save('ave_power_all_freq_cfc_cn.mat','ave_power_all_freq_cfc_cn','-v6')
% save('ave_bold_cfc_ad.mat','ave_bold_cfc_ad','-v6')
% save('ave_power_all_freq_cfc_ad.mat','ave_power_all_freq_cfc_ad','-v6')
% 
% %% 20230512 exploring
% upsignalFile = dir(fullfile('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Identified pipeline server\upsampled signal\','*.mat'));
% upsignalName = {upsignalFile.name}';
% 
% c = load('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Localized_Spectral_Graph_Filter_Frames_Code\aal116_fMRI.mat');
% for sample_bold_i = 1:135
%     sample_bold = c.aal116_fMRI(sample_bold_i).BOLD;
%     fMRI_cov = cov(sample_bold(:,1:90));
% end
% 
% filename = ['D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Localized_Spectral_Graph_Filter_Frames_Code\subdict.mat'];
% a = load(filename);
% structname = fieldnames(a);
% subdict = getfield(a,structname{1});
% 
% ave_cfc = cell(90,1);
% for k = 1:90
%     ave_cfc{k} = zeros(90,90);
% end
% 
% all_cfc = zeros(8100,8100);
% for i = 1:size(upsignalName,1)
%     all_power = [];
%     split_upsignalName = strsplit(upsignalName{i,1},'_');
%     split_scale = split_upsignalName{1,2};
%     split_sample = split_upsignalName{1,3};
%     split_sample = strsplit(split_sample,'.');
%     split_sample = split_sample{1,1};
%     if strcmp(split_scale,num2str(3))
%         for j=1:90
%             temp_filename = ['D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Identified pipeline server\upsampled signal\',upsignalName{i,1}];
%             b = load(temp_filename);
%             temp_signal = real(b.upsampled_signal);
%             temp_power = real(subdict{str2num(split_sample),j})*temp_signal;
% %             for 
% %             end
% 
% %             % identify frequency
% %             color=[1 0 0;0 1 0;0 0 1;0.5 1 1;1 1 0.5;1 0.5 1; 0 0 0.5; 0.5 0 0;0 0.5 0;1 0.5 0.5; 0.5 1 0.5;0.5 0.5 1;1 1 0;0 1 1;1 0 1];
% %             Fs = 50;            % Sampling frequency                    
% %             L = size(b.upsampled_signal(1,:),2);             % Length of signal
% % 
% %             S = temp_power; % Graph_temp(j,1).Harmonics';
% %             Y=fft(S(1,1:L));
% %             P2=abs(Y/L);
% %             P1=P2(1:L/2+1);
% %             P1(2:end-1)=2*P1(2:end-1);
% % 
% %             f = Fs*(0:(L/2))/L;
% %             plot(f,P1,'color',color(mod(j,15)+1,:))
% %             hold on
% %             title("Amplitude Spectrum of BOLD signals")
% %             xlabel("f (Hz)")
% %             ylabel("|P1(f)|")
% 
%             % calculate cfc
%             temp_cfc = cov(temp_power');
%             ave_cfc{j} = ave_cfc{j} + temp_cfc;
%             all_power = [all_power;temp_power(:,1:10000)];
%         end
%     end
%     temp_all_cfc = cov(all_power');
%     if sum(sum(1*isnan(temp_all_cfc)))==0
%         all_cfc = all_cfc + temp_all_cfc;
%     end
% end
% 
% cfc = zeros(8100,8100);
% for i=1:8100
%     for j=1:8100
%         b_i = mod(i-1,90)+1;
%         b_ii = (i-b_i)/90;
%         b_j = mod(j-1,90)+1;
%         b_jj = (j-b_j)/90;
%         cfc(b_ii+90*b_i,b_jj+90*b_j)=temp_all_cfc(i,j);
%     end
% end
% figure
% imagesc(histeq(cfc))
% figure
% imagesc(histeq(temp_all_cfc))
% 
% save('all_power.mat','all_power','-v6')
% save('all_cfc.mat','all_cfc','-v6')
% 
% for k = 10:10:90
% %     ave_cfc{k} = ave_cfc{k}/count_scale_num(k);
%     figure
%     for l=1:90
%         ave_cfc{k}(l,l)=0;
%     end
%     imagesc(ave_cfc{k})
% end
% 
% 
% b = load('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Localized_Spectral_Graph_Filter_Frames_Code\aal116_fMRI.mat');
% 
% figure
% Fs = 50;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 100*size(b.aal116_fMRI(1).BOLD,1);             % Length of signal
% t = (0:L-1)*T;        % Time vector
% 
% color=[1 0 0;0 1 0;0 0 1;0.5 1 1;1 1 0.5;1 0.5 1; 0 0 0.5; 0.5 0 0;0 0.5 0;1 0.5 0.5; 0.5 1 0.5;0.5 0.5 1;1 1 0;0 1 1;1 0 1];
% 
% for i = 1:5
%     for j = 1:3
%         subplot(5,3,j+3*(i-1))
%         scatter(t(1:100:19600),b.aal116_fMRI(1).BOLD(1:196,i),5,color(1,:)) %/max(b.aal116_fMRI(1).BOLD)
%         hold on
%         plot(t(1:19600),upsignal{1,j}(i,1:19600),'color',color(2,:)) %/max(upsignal{1,j}(i,1:19600))
%         hold on
%         filename = ['D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\Localized_Spectral_Graph_Filter_Frames_Code\subdict',num2str(j),'.mat'];
%         a = load(filename);
%         structname = fieldnames(a);
%         subdict_all(:,j) = getfield(a,structname{1});
%         temp{j}=subdict_all{1,j}*upsignal{1,j};
%         S=temp{j}(i,:); %/max(temp(i,:))
%         plot(t(1:19600),S(1:19600),'color',color(3,:))
%         set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
%     end
% end
% 
% ave_cfc = cell(3,1);
% for i=1:1
%     ave_cfc{i} = zeros(90,90);
%     for sliding_window=1:99
%         cfc{sliding_window,i}=cov(temp{i}(:,100*(sliding_window-1)+1:100*(sliding_window-1)+200)');
%         ave_cfc{i} = ave_cfc{i} + cfc{sliding_window,i};
%     end
%     ave_cfc{i} = ave_cfc{i}/100;
% end
% 
% save('ave_cfc.mat','ave_cfc','-v6')
% 
% figure
% subplot(2,1,2)
% for j = 1:3
%     S = temp{1,4-j}; % Graph_temp(j,1).Harmonics';
%     Y=fft(S(1,1:19600));
%     P2=abs(Y/L);
%     P1=P2(1:L/2+1);
%     P1(2:end-1)=2*P1(2:end-1);
% 
%     f = Fs*(0:(L/2))/L;
%     plot(f,P1,'color',color(4-j,:))
%     hold on
%     title("Amplitude Spectrum of BOLD signals")
%     xlabel("f (Hz)")
%     ylabel("|P1(f)|")
% end
% legend(leg_str)
% 
% save('upsignal.mat','upsignal','-v6')
% figure
% plot(upsignal{1, 1}(:,1))
% % for k=1:90
% %     temp(:,k)=band_filter_signal(upsignal{1}(:,k),50);
% % end
% 
% 
% filter_temp=cell(135,3);
% cfc=cell(135,3);
% temp_cfc=cell(135,1);
% ave_cfc = zeros(270,270);
% for sample_i=1:3
%     temp_temp=[];
%     for i=1:3
% %         for k=1:90
% %             filter_temp{sample_i,i}(k,:)=band_filter_single(temp{sample_i,i}(k,:),50,0.5*i-0.1,0.5*i+0.1);
% %         end
%         filter_temp=temp;
%         cfc{sample_i,i}=cov(filter_temp{sample_i,i}');
%         temp_temp=[temp_temp,filter_temp{sample_i,i}'];
%     end
%     temp_cfc{sample_i} = cov(temp_temp);
%     temp_cfc{sample_i}(temp_cfc{sample_i}<0)=0;
%     ave_cfc = ave_cfc + temp_cfc{sample_i};
% end
% ave_cfc=ave_cfc/3;
% 
% % figure
% % imagesc(cfc{1,1})
% % figure
% % imagesc(cfc{1,2})
% % figure
% % imagesc(cfc{1,3})
% % figure
% % imagesc(temp_cfc{1})
% figure
% imagesc(ave_cfc)
% set(gca, 'Fontname', 'Times New Roman','FontSize',20,'FontWeight','bold');
% colormap(jet)