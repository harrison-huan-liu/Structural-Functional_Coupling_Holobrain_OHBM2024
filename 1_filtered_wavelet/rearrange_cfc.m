clc
clear
close all

% post_process_CFC()
% calculate_CFC()

current_file = mfilename('fullpath');
[current_path, ~, ~] = fileparts(current_file);
[current_path, ~, ~] = fileparts(current_path);

HCP_task_data_path = fullfile(current_path, 'HCP_YA', 'HCP_task_data_test_labeled_4.mat');
load(HCP_task_data_path);
% FC_matrix_4 = calculate_ave_FC(HCP_task_data_test_labeled_4.HCP_task_data);
% HCP_task_data_test_labeled_1 = load('C:\Holobrain\HCP_YA\HCP_task_data_retest_labeled_1.mat');
% FC_matrix_1 = calculate_ave_FC(HCP_task_data_test_labeled_1.HCP_task_data);
% 
% FC_matrix_d = FC_matrix_1 - FC_matrix_4;
% figure
% imagesc(FC_matrix_d)
% clim([-max(max(FC_matrix_d)) max(max(FC_matrix_d))])
% colormap("jet")
% colorbar

HCP_MMP_360_path = fullfile(current_path, 'HCP_YA', 'HCP_MMP_360.xlsx');
HCP_MMP_360 = readcell(HCP_MMP_360_path);
Motor_cortex = ["Somatosensory_and_Motor", "Premotor"];
Motor_node = zeros(360,1);
for node_i = 2:361
    if strcmp(HCP_MMP_360{node_i,7},Motor_cortex{1})
        Motor_node(node_i-1,1) = 1;
    end
end
Motor_node = ones(360,1);


subtask_label = unique(HCP_task_data(:,2));

Network_path = fullfile(current_path, 'HCP_YA', 'MMP_DiffusionConnectivity_HCP_avg56.mat');
load(Network_path)
% not_include_index = [2,3,5];
% not_include_index = [2,3,4,5,6,7,8,10,11,13];
not_include_index = [];
wavelet_num = 360-size(not_include_index,2);
filtered_signal = cell(2,2);
filtered_signal_node = cell(2,2);

SC_avg56(SC_avg56<0.5)=0;
path = zeros(360,360);
diag_line = diag(ones(360,1));
SC_avg56 = SC_avg56 + diag_line;
ave_cfc = cell(2,2);
for rad_i = 3:10
    fprintf('Rad_i: %d!\n', rad_i)
    path = SC_avg56^rad_i;
    
    path(path>0)=1;
    wavelet_path = fullfile(current_path, 'output', 'subdict', 'wavelet_360_single_original_pn.mat');
    load(wavelet_path)
    
    for i = 1:360
        subdict{1,i} = subdict{1,i}.*path;
    end
    
    pic_i = 1;
    for figure_i=1:360
        % if rem(pic_i-1,36)==0
        %     figure
        % end
        % set(gcf, 'unit', 'centimeters', 'position', [10 5 45 30]);
        if ~ismember(figure_i,not_include_index)
            % subplot(6,6,rem(pic_i-1,36)+1)
            for sample_num = 1:size(subtask_label,1):size(HCP_task_data,1)
                filtered_signal{floor(sample_num/size(subtask_label,1))+1,pic_i} = subdict{1,figure_i}*HCP_task_data{sample_num,1}';
            end
            % imagesc(subdict{1,figure_i}) % filtered_signal{1,pic_i}
            % clim([-max(max(subdict{1,figure_i})) max(max(subdict{1,figure_i}))])
            % colormap("jet")
            % colorbar
            pic_i = pic_i + 1;
        end
    end
    
    for sample_i = 1:size(filtered_signal,1)
        for node_i = 1:360
            filtered_signal_node{sample_i,node_i} = zeros(wavelet_num,size(HCP_task_data{1,1},1)); % 100
            for fre_i = 1:size(filtered_signal,2)
                filtered_signal_node{sample_i,node_i}(fre_i,:) = filtered_signal{sample_i,fre_i}(node_i,1:size(HCP_task_data{1,1},1)); % 1:100
            end
        end
    end
    
    CFC_pattern = cell(size(filtered_signal,1),size(filtered_signal_node,2));
    
    for plt_i = 1:size(filtered_signal_node,2)
        ave_cfc{plt_i,rad_i} = zeros(wavelet_num,wavelet_num);
        % if rem(plt_i-1,36)==0
        %     figure
        % end
        % set(gcf, 'unit', 'centimeters', 'position', [10 5 45 30]);
        % subplot(6,6,rem(plt_i-1,36)+1)
        for sample_j = 1:size(filtered_signal,1)
            % CFC_pattern{sample_j,plt_i} = corr(filtered_signal_node{sample_j,plt_i}', 'type', 'Pearson');
            CFC_pattern{sample_j,plt_i} = filtered_signal_node{sample_j,plt_i}*filtered_signal_node{sample_j,plt_i}';
            ave_cfc{plt_i,rad_i} = ave_cfc{plt_i,rad_i} + CFC_pattern{sample_j,plt_i};
        end
        ave_cfc{plt_i,rad_i} = ave_cfc{plt_i,rad_i}/size(filtered_signal,1);
        % show_mat = ave_cfc{plt_i};
        % % show_mat(abs(show_mat)<0.3)=0;
        % imagesc(show_mat)
        % temp_matrix = ave_cfc{plt_i};
        % temp_max = temp_matrix - diag(diag(temp_matrix));
        % % for repeat = 1:10
        % %     temp_max(find(temp_max == max(temp_max))) = mean(temp_max);
        % % end
        % clim([-max(max(temp_max)) max(max(temp_max))])
        % colormap("jet")
        % colorbar
    end
    
    node_ave_cfc = zeros(wavelet_num,wavelet_num);
    figure
    for i = 1:360
        if Motor_node(i) == 1
            node_ave_cfc = node_ave_cfc + ave_cfc{i,rad_i};
        end
    end
    node_ave_cfc = node_ave_cfc/360;
    % Nr = normalize(nthroot(node_ave_cfc,3),"range");
    % imagesc(Nr*2-1) % 
    imagesc(nthroot(node_ave_cfc,3)) % 
    max_value = max(max(node_ave_cfc));
    set(gca, 'Fontname', 'Times New Roman','FontSize',18,'FontWeight','bold');
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    c = colorbar;
    % c.Ticks = [-1 -0.5 0 0.5 1];
    % c.TickLabels = {'-1','-0.125','0','0.125','1'};
    c.Limits = [-max_value,max_value];
    colormap("jet")
    filename = ['CFC_', num2str(rad_i), '.jpg'];
    filepath = fullfile(current_path, 'output', filename);
    saveas(gcf,filepath)
end



function [FC_matrix] = calculate_ave_FC(HCP_task_data)
subtask_label = unique(HCP_task_data(:,2));
FC_matrix = zeros(360,360);
for i = 1:size(subtask_label, 1):size(HCP_task_data,1)
    temp_FC_matrix = HCP_task_data{i, 1}'*HCP_task_data{i, 1};
    FC_matrix = FC_matrix + temp_FC_matrix/size(HCP_task_data{i, 1},1); % 
end
FC_matrix = FC_matrix/(size(HCP_task_data,1)/size(subtask_label, 1));

figure
imagesc(FC_matrix)
clim([-max(max(FC_matrix)) max(max(FC_matrix))])
colormap("jet")
colorbar
end

function [] = calculate_CFC()
wavelet = load('wavelet_360.mat');

% node_select = [60,120];
node_select = 1:360;
CFC_node = cell(size(node_select,2),1);

subtask_label = unique(HCP_task_data(:,2));
power = cell(1,size(HCP_task_data,1)/size(subtask_label, 1));

for i = 1:size(node_select,2)
    CFC_node{i,1} = zeros(360,360);
end

for node_i = 1:size(node_select,2)
    for j = 1:(size(HCP_task_data,1)/size(subtask_label, 1))
        power{1,j} = zeros(360,size(HCP_task_data{1,1},1));
    end
    fprintf("The node: %d!\n", node_select(node_i))
    upd = textprogressbar(360+(size(HCP_task_data,1)/size(subtask_label, 1)), 'barlength', 20, ...
                         'updatestep', 10, ...
                         'startmsg', 'Waiting... ',...
                         'endmsg', ' Yay!', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '+', ...
                         'emptybarsymbol', '-');
    ave_cfc = zeros(360,360);
    for wavelet_i = 1:360
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
        temp_corr = power{1,i}*power{1,i}'; % corr(power{1,i}', 'type', 'pearson');
        ave_cfc = ave_cfc + temp_corr;
        upd(360+i)
    end
    CFC_node{node_i} = ave_cfc/(size(HCP_task_data,1)/size(subtask_label, 1));
end

figure
set(gcf, 'unit', 'centimeters', 'position', [10 5 40 15]);
for node_i = 1:size(node_select,2)
    subplot(1,2,node_i)
    imagesc(CFC_node{node_i})
    colorbar
end
end


function [] = post_process_CFC()
% CFC_node_retest_1 = load("CFC_node_retest_1.mat");
% CFC_node_retest_4 = load("CFC_node_retest_4.mat");
% CFC_node_test_1 = load("CFC_node_test_1.mat");
% CFC_node_test_4 = load("CFC_node_test_4.mat");

CFC_node_retest_1 = load("CFC_node_retest_1_single.mat");
CFC_node_retest_4 = load("CFC_node_retest_4_single.mat");
CFC_node_test_1 = load("CFC_node_test_1_single.mat");
CFC_node_test_4 = load("CFC_node_test_4_single.mat");

% plot_max_pooling_CFC(CFC_node_retest_1.CFC_node);
% plot_max_pooling_CFC(CFC_node_retest_4.CFC_node);
% plot_max_pooling_CFC(CFC_node_test_1.CFC_node);
% plot_max_pooling_CFC(CFC_node_test_4.CFC_node);

plot_ave_CFC(CFC_node_retest_1.CFC_node);
plot_ave_CFC(CFC_node_retest_4.CFC_node);
plot_ave_CFC(CFC_node_test_1.CFC_node);
plot_ave_CFC(CFC_node_test_4.CFC_node);

plot_node_CFC(CFC_node_retest_1.CFC_node);
plot_node_CFC(CFC_node_retest_4.CFC_node);
plot_node_CFC(CFC_node_test_1.CFC_node);
plot_node_CFC(CFC_node_test_4.CFC_node);


[r_cfc_a,plot_temp_cfc_all_a,index_pos_a,index_neg_a] = plot_rearrange_cfc(1,1,1,360,360);
% [r_cfc_b,plot_temp_cfc_all_b,index_pos_b,index_neg_b] = plot_rearrange_cfc(1,1,1,360,36);

com_element_pos_a = find_common_element(index_pos_a);
com_element_neg_a = find_common_element(index_neg_a);

union_element_pos_a = union_index(index_pos_a);
union_element_neg_a = union_index(index_neg_a);

time_table_pos_a = count_element_appear_time(index_pos_a);
time_table_neg_a = count_element_appear_time(index_neg_a);

time_diff = time_table_pos_a-time_table_neg_a;
time_diff_select_pos = find(time_diff(:,2)>60);
time_diff_select_neg = find(time_diff(:,2)<-40);

time_diff_value_pos = time_diff(time_diff_select_pos,2);
time_diff_value_neg = time_diff(time_diff_select_neg,2);

identify_index_internal = identify_internal_para(time_diff);
end


function [] = plot_node_CFC(CFC_node)
for figure_i = 1:floor(size(CFC_node,1)/36)
    figure % figure(figure_i)
    set(gcf, 'unit', 'centimeters', 'position', [10 1 45 35]);
    for subfigure_i = 1:6
        for subfigure_j = 1:6
            subplot(6,6,subfigure_i+(subfigure_j-1)*6)
            plot_CFC_node = CFC_node{subfigure_i+(subfigure_j-1)*6+(figure_i-1)*36,1};
            plot_CFC_node = plot_CFC_node - diag(diag(plot_CFC_node));
            imagesc(plot_CFC_node) % nthroot(,3)
            colorbar
            % clim([-0.1 0.1])
        end
    end
end
end


function [ave_node_cfc,select_ave_cfc] = plot_ave_CFC(CFC_node,wavelet)
% first representation
wavelet_num = size(wavelet,1);
node_num = size(CFC_node,1);
ave_node_cfc = zeros(wavelet_num,wavelet_num);
for i = 1:5:node_num
    ave_node_cfc = ave_node_cfc + CFC_node{i};
end
ave_node_cfc = ave_node_cfc / node_num;
ave_node_cfc = (ave_node_cfc - diag(diag(ave_node_cfc)));
% ave_node_cfc(abs(ave_node_cfc)<0.001)=0;
figure
imagesc(ave_node_cfc) % nthroot(ave_node_cfc,3)
clim([-max(max(ave_node_cfc)) max(max(ave_node_cfc))])
colormap("jet")
colorbar

% second representation
select_ave_cfc = zeros(wavelet_num,wavelet_num);
for j = 1:size(wavelet,1)
    ind = (diag(wavelet{j,1})>mean(wavelet{j,1})');
    for k = 1:size(wavelet{j,1},1)
        CFC_node{k}(j,:) = CFC_node{k}(j,:)*ind(k);
        CFC_node{k}(:,j) = CFC_node{k}(:,j)*ind(k);
    end
end
for l = 1:5:node_num
    select_ave_cfc = select_ave_cfc + CFC_node{l};
end
select_ave_cfc = select_ave_cfc / node_num;
select_ave_cfc = (select_ave_cfc - diag(diag(select_ave_cfc)));
% select_ave_cfc(abs(select_ave_cfc)<0.001)=0;
figure
imagesc(ave_node_cfc) % nthroot(select_ave_cfc,3)
clim([-max(max(select_ave_cfc)) max(max(select_ave_cfc))])
colormap("jet")
colorbar
end


function [] = plot_max_pooling_CFC(CFC_node)
max_node_cfc = zeros(360,360);
for i = 1:360
    for j = 1:360
        for k = 1:360
            if abs(max_node_cfc(j,k))<abs(CFC_node{i}(j,k))
                max_node_cfc(j,k) = CFC_node{i}(j,k);
            end
        end
    end
end
figure
% max_node_cfc(abs(max_node_cfc)<0.3)=0;

imagesc(max_node_cfc(1:360,10))
clim([-1 1])
colormap("jet")
colorbar
end



function [identify_index_internal] = identify_internal_para(time_diff)
identify_index_internal = 0;
loss = 0;
for index_internal = 1:5
    count_i = 0;
    temp_loss = 0;
    internal_seq = 1:index_internal:18;
    for i = 1:2:size(internal_seq,2)
        if i<=size(internal_seq,2)-2
            temp_seq_a = (internal_seq(i))^2;
            temp_seq_b = (internal_seq(i+1))^2;
            temp_seq_c = (internal_seq(i+2))^2;
            temp_loss_a = sum(time_diff(temp_seq_a:temp_seq_b,2));
            temp_loss_b = sum(time_diff(temp_seq_b:temp_seq_c,2));
            temp_loss = temp_loss + abs(temp_loss_a - temp_loss_b);
            count_i = count_i + 1;
        end
    end
    temp_loss = temp_loss % /count_i
    if loss<temp_loss
        identify_index_internal = index_internal;
        loss=temp_loss;
    end
end
end


function [time_table] = count_element_appear_time(index)
all_index = index{1,1};
for k=2:size(index,1)
    all_index = [all_index;index{k,1}];
end
time_table = tabulate(all_index(:));
end


function [union_element] = union_index(index)
union_element = union(index{1,1},index{2,1});
for k=3:size(index,1)
    union_element = union(union_element,index{k,1});
end
end


function [com_element] = find_common_element(index)
com_element = intersect(index{1,1},index{2,1});
for k=3:size(index,1)
    com_element = intersect(com_element,index{k,1});
    if isempty(com_element)
        break;
    end
end
end


function [r_cfc,plot_temp_cfc_all,index_pos,index_neg] = plot_rearrange_cfc(task_begin,task_end,plot_diag_only,wavelet_num,draw_index)
r_cfc = cell(360,1);
index_pos = cell(360,1);
index_neg = cell(360,1);
current_file = mfilename('fullpath');
[current_path, ~, ~] = fileparts(current_file);
for task_i = task_begin:task_end % 1:7
    cfc_path = [current_path, '\ave_cfc_wavelet_general_retest_', num2str(task_i), '_', num2str(wavelet_num),'_1.mat'];
    cfc = load(cfc_path);
    for subtask_i = 1:1 % size(cfc.ave_cfc,1)
        for k=1:360
            r_cfc{k,1} = zeros(wavelet_num, wavelet_num);
            for i=1:wavelet_num
                for j=1:wavelet_num
                    temp_cfc = cfc.ave_cfc{subtask_i,i+(j-1)*wavelet_num};
                    r_cfc{k,1}(i,j) = temp_cfc(k,1);
                end
            end
            index_pos{k,1}=find(r_cfc{k,1}(:,1)>0);
            index_neg{k,1}=find(r_cfc{k,1}(:,1)<0);
        end

        if plot_diag_only
            plot_temp_cfc_all = zeros(wavelet_num,wavelet_num);
            for figure_i = 1:1 % 10
                figure % figure(figure_i)
                set(gcf, 'unit', 'centimeters', 'position', [10 1 45 35]);
                for subfigure_i = 1:6
                    for subfigure_j = 1:6
                        subplot(6,6,subfigure_i+(subfigure_j-1)*6)
                        plot_temp_cfc = r_cfc{subfigure_i+(subfigure_j-1)*6+(figure_i-1)*36,1};

                        % plot_temp_cfc = plot_temp_cfc-diag(diag(plot_temp_cfc));
                        % plot_temp_cfc = (plot_temp_cfc-min(min(plot_temp_cfc)))/(max(max(plot_temp_cfc))-min(min(plot_temp_cfc)));
                        % plot_temp_cfc(plot_temp_cfc>0.3)=0.3;
                        plot_temp_cfc_all = plot_temp_cfc_all + plot_temp_cfc;
                        % thresold_value = mean(plot_temp_cfc,"all"); % 0.3;
                        % plot_temp_cfc(plot_temp_cfc>thresold_value) = 1;
                        % plot_temp_cfc(plot_temp_cfc<thresold_value) = 0;
                        imagesc(plot_temp_cfc(:,1:10))
                        colorbar
                        % plot(plot_temp_cfc(1:50,1))
                    end
                end
            end
            figure
            % all_thresold_value = mean(plot_temp_cfc_all,"all"); % 180
            % plot_temp_cfc_all(plot_temp_cfc_all>all_thresold_value) = 360;
            % plot_temp_cfc_all(plot_temp_cfc_all<all_thresold_value) = 0;
            % plot_temp_cfc_all = plot_temp_cfc_all + diag(360*ones(20,1));
            imagesc(plot_temp_cfc_all(1:draw_index,1:draw_index))
            colorbar
            % plot(plot_temp_cfc_all(:,1))
        else
            for figure_i = 1:16
                figure(figure_i)
                set(gcf, 'unit', 'centimeters', 'position', [10 1 45 35]);
                figcol = mod(figure_i-1,4);
                figrow = floor((figure_i-1)/4);
                for subfigure_i = 1:10
                    for subfigure_j = 1:10
                        subplot(10,10,subfigure_i+(subfigure_j-1)*10)
                        plot_temp_cfc = r_cfc(1+10*(subfigure_j-1+figrow*4):10+10*(subfigure_j-1+figrow*4),1+10*(subfigure_i-1+figcol*4):10+10*(subfigure_i-1+figcol*4));
                        plot_temp_cfc(1,1) = 0;
                        plot_temp_cfc(1,1) = max(max(plot_temp_cfc));
                        plot_temp_cfc = (plot_temp_cfc-min(min(plot_temp_cfc)))/(max(max(plot_temp_cfc))-min(min(plot_temp_cfc)));
                        imagesc(plot_temp_cfc)
                    end
                end
            end
        end
    end
end

end



function upd = textprogressbar(n, varargin)
% UPD = TEXTPROGRESSBAR(N) initializes a text progress bar for monitoring a
% task comprising N steps (e.g., the N rounds of an iteration) in the
% command line. It returns a function handle UPD that is used to update and
% render the progress bar. UPD takes a single argument i <= N which 
% corresponds to the number of tasks completed and renders the progress bar
% accordingly.
%                   
% TEXTPROGRESSBAR(...,'barlength',L) determines the length L of the
% progress bar in number of characters (see 'barsymbol' option). L must be
% a positive integer. 
% (Default value is 20 characters.)
%
% TEXTPROGRESSBAR(...,'updatestep',S) determines the minimum number of update
% steps S between consecutive bar re-renderings. The option controls how
% frequently the bar is rendered and in turn controls the computational 
% overhead that the bar rendering incurs to the monitored task. It is
% especially useful when bar is used for loops with large number of rounds
% and short execution time per round.
% (Default value is S=10 steps.)
%
% TEXTPROGRESSBAR(...,'startmsg',str) determines the message string to be
% displayed before the progress bar.
% (Default is str='Completed '.)
%
% TEXTPROGRESSBAR(...,'endmsg',str) determines the message string to be 
% displayed after progress bar when the task is completed.
% (Default is str=' Done.')
%
% TEXTPROGRESSBAR(...,'showremtime',b) logical parameter that controls
% whether an estimate of the remaining time is displayed.
% (Default is b=true.)
%
% TEXTPROGRESSBAR(...,'showbar',b) logical parameter that controls whether
% the progress bar is displayed. (Default is b=true.)
%
% TEXTPROGRESSBAR(...,'showpercentage',b) logical parameter that controls
% whether to display the percentage of completed items.
% (Default is true.)
%
% TEXTPROGRESSBAR(...,'showactualnum',b) logical parameter that controls
% whether to display the actual number of completed items.
% (Default is false.)
%
% TEXTPROGRESSBAR(...,'showfinaltime',b) logical parameter that controls
% whether to display the total run-time when completed.
% (Default is true.)
%
% TEXTPROGRESSBAR(...,'barsymbol',c) determines the symbol (character) to
% be used for the progress bar. c must be a single character.
% (Default is c='='.)
%
% TEXTPROGRESSBAR(...,'emptybarsymbol',c) determines the symbol (character)
% that is used to fill the un-completed part of the progress bar. c must be
% a single character.
% (Default is c=' '.)
%
% Example:
%
%   n = 150;
%   upd = textprogressbar(n);
%   for i = 1:n
%      pause(0.05);
%      upd(i);
%   end
%

    % Default Parameter values:
    defaultbarCharLen = 20;
    defaultUpdStep = 10;
    defaultstartMsg = 'Completed ';
    defaultendMsg = ' Done.';
    defaultShowremTime = true;
    defaultShowBar = true;
    defaultshowPercentage = true;
    defaultshowActualNum = false;
    defaultshowFinalTime = true;
    defaultbarCharSymbol = '=';
    defaultEmptybarCharSymbol = ' ';
    
    % Auxiliary functions for checking parameter values:
    ischarsymbol = @(c) (ischar(c) && length(c) == 1);
    ispositiveint = @(x) (isnumeric(x) && mod(x, 1) == 0 && x > 0);
 
    % Register input parameters:
    p = inputParser;
    addRequired(p,'n', ispositiveint);
    addParameter(p, 'barlength', defaultbarCharLen, ispositiveint)
    addParameter(p, 'updatestep', defaultUpdStep, ispositiveint)
    addParameter(p, 'startmsg', defaultstartMsg, @ischar)
    addParameter(p, 'endmsg', defaultendMsg, @ischar)
    addParameter(p, 'showremtime', defaultShowremTime, @islogical)
    addParameter(p, 'showbar', defaultShowBar, @islogical)
    addParameter(p, 'showpercentage', defaultshowPercentage, @islogical)
    addParameter(p, 'showactualnum', defaultshowActualNum, @islogical)
    addParameter(p, 'showfinaltime', defaultshowFinalTime, @islogical)
    addParameter(p, 'barsymbol', defaultbarCharSymbol, ischarsymbol)
    addParameter(p, 'emptybarsymbol', defaultEmptybarCharSymbol, ischarsymbol)
    
    % Parse input arguments:
    parse(p, n, varargin{:});
    n = p.Results.n;
    barCharLen = p.Results.barlength;
    updStep = p.Results.updatestep;
    startMsg = p.Results.startmsg;
    endMsg = p.Results.endmsg;
    showremTime = p.Results.showremtime;
    showBar = p.Results.showbar;
    showPercentage = p.Results.showpercentage;
    showActualNum = p.Results.showactualnum;
    showFinalTime = p.Results.showfinaltime;
    barCharSymbol = p.Results.barsymbol;
    emptybarCharSymbol = p.Results.emptybarsymbol;
    
    % Initialize progress bar:
    bar = ['[', repmat(emptybarCharSymbol, 1, barCharLen), ']'];
    
    nextRenderPoint = 0;
    startTime = tic;
    
    % Initalize block for actual number of completed items:
    
    ind = 1;
    
    % Start message block:
    startMsgLen = length(startMsg);
    startMsgStart = ind;
    startMsgEnd = startMsgStart + startMsgLen - 1;
    ind = ind + startMsgLen;
    
    % Bar block:
    barLen = length(bar);
    barStart = 0;
    barEnd = 0;
    if showBar
        barStart = ind;
        barEnd = barStart + barLen - 1;
        ind = ind + barLen;
    end
    
    % Actual Num block:
    actualNumDigitLen = numel(num2str(n));
    actualNumFormat = sprintf(' %%%dd/%d', actualNumDigitLen, n);
    actualNumStr = sprintf(actualNumFormat, 0);
    actualNumLen = length(actualNumStr);
    actualNumStart = 0;
    actualNumEnd = 0;
    if showActualNum
        actualNumStart = ind;
        actualNumEnd = actualNumStart + actualNumLen-1;
        ind = ind + actualNumLen;
    end
        
    % Percentage block:
    percentageFormat = sprintf(' %%3d%%%%');
    percentageStr = sprintf(percentageFormat, 0);
    percentageLen = length(percentageStr);
    percentageStart = 0;
    percentageEnd = 0;
    if showPercentage
        percentageStart = ind;
        percentageEnd = percentageStart + percentageLen-1;
        ind = ind + percentageLen;
    end
    
    % Remaining Time block:
    remTimeStr = time2str(Inf);
    remTimeLen = length(remTimeStr);
    remTimeStart = 0;
    remTimeEnd = 0;
    if showremTime
       remTimeStart = ind;
       remTimeEnd = remTimeStart + remTimeLen - 1;
       ind = ind + remTimeLen;
    end
    
    
    % End msg block:
    endMsgLen = length(endMsg);
    if showBar
        endMsgStart = barEnd + 1; % Place end message right after bar;
    else
        endMsgStart = startMsgEnd + 1;
    end
    endMsgEnd = endMsgStart + endMsgLen - 1;
    
    ind = max([ind, endMsgEnd]);
    
    % Determine size of buffer:
    arrayLen = ind - 1;
    array = repmat(' ', 1, arrayLen);
    
    % Initial render:
    array(startMsgStart:startMsgEnd) = sprintf('%s', startMsg);

    delAll = repmat('\b', 1, arrayLen);
    
        % Function to update the status of the progress bar:
        function update(i)
            
            if i < nextRenderPoint
                return;
            end
            if i > 0
                fprintf(delAll);
            end
            %pause(1)
            nextRenderPoint = min([nextRenderPoint + updStep, n]);
            
            if showremTime
                % Delete remaining time block:
                array(remTimeStart:remTimeEnd) = ' ';
            end

            if showPercentage
                % Delete percentage block:
                array(percentageStart:percentageEnd) = ' ';
            end
            
            if showActualNum
                % Delete actual num block:
                array(actualNumStart:actualNumEnd) = ' ';
            end
    
            if showBar
                % Update progress bar (only if needed):
                barsToPrint = floor( i / n * barCharLen );
                bar(2:1+barsToPrint) = barCharSymbol;
                array(barStart:barEnd) = bar;
            end
            
            % Check if done:
            if i >= n
                array(endMsgStart:endMsgEnd) = endMsg;
                array(endMsgEnd+1:end) = ' ';
                
                if showFinalTime
                    finalTimeStr = ...
                        sprintf(' [%d seconds]', round(toc(startTime)));
                    finalTimeLen = length(finalTimeStr);
                    if endMsgEnd + finalTimeLen < arrayLen
                       array(endMsgEnd+1:endMsgEnd+finalTimeLen) = ... 
                           finalTimeStr;
                    else
                       array = [array(1:endMsgEnd), finalTimeStr];
                    end
                end
                
                fprintf('%s', array);
                fprintf('\n');
                return;
            end
            
            if showActualNum
                % Delete actual num block:
                actualNumStr = sprintf(actualNumFormat, i);
                array(actualNumStart:actualNumEnd) = actualNumStr;
            end
            
            if showPercentage
                % Render percentage block:
                percentage = floor(i / n * 100);
                percentageStr = sprintf(percentageFormat, percentage);
                array(percentageStart:percentageEnd) = percentageStr;
            end
                
            % Print remaining time block:
            if showremTime
               t = toc(startTime);
               remTime = t/ i * (n-i);
               remTimeStr = time2str(remTime);
               array(remTimeStart:remTimeEnd) = remTimeStr;
            end
            fprintf('%s', array);
        end
    
    % Do the first render:
    update(0);
    
    upd = @update;
    
end

% Auxiliary functions

function timestr = time2str(t)

    if t == Inf
        timestr = sprintf(' --:--:--');
    else
        [hh, mm, tt] = sec2hhmmss(t);
        timestr = sprintf(' %02d:%02d:%02d', hh, mm, tt);
    end
end

function [hh, mm, ss] = sec2hhmmss(t)
    hh = floor(t / 3600);
    t = t - hh * 3600;
    mm = floor(t / 60);
    ss = round(t - mm * 60);
end
