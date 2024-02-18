function [current_path] = Filter_identify_general(filter_bandnumber,wavelet_sperate,normal_or_single,radius)
    % close all;
    % clear all;
    % clc;
    % wavelet_sperate = 1;
    % normal_or_single = 1;

    current_file = mfilename('fullpath');
    [current_path, ~, ~] = fileparts(current_file);
    [current_path, ~, ~] = fileparts(current_path);

    % rand('seed',0);
    % randn('seed',0);

    HCP_sMRI = struct;
    % filter_bandnumber = 3;
    subdict = cell(2,filter_bandnumber);

    for sample_i=1:1
        structnetwork_path = fullfile(current_path,'HCP_YA','MMP_DiffusionConnectivity_HCP_avg56.mat');
        HCP_network = load(structnetwork_path);
        W = HCP_network.SC_avg56;
        % W(W<0.01) = 0;
        W = W + diag(ones(360,1));
        [W,thr,knum] = density_thresholding(W,0.1,3);
        path = W^radius;
        path(path>0) = 1;
        W = W-diag(diag(W));

        if size(W,1) ~= size(W,2), error('Weight matrix W is not square'); end

        if (norm(W - W', 'fro') == 0)
            disp('The matrix W is symmetric');
        else
            W = max(W, W');;
        end

        % bool = sum(sum(abs(W - transpose(W))> eps(10) ))>0;
        % fprintf('W is direct: %d !', num2str(bool))
        
        node_num = size(W,1);
        HCP_sMRI(sample_i).W = W;
        % max_X = max(max(HCP_sMRI(sample_i).W(1:node_num,1:node_num)));
        % HCP_sMRI(sample_i).W = HCP_sMRI(sample_i).W/max_X;

        D = sum(W,2);

        % The first normalized
        ind = D>0;
        L = -W;
        L(ind, :) = bsxfun(@times, L(ind, :), 1./sqrt(D(ind)));
        L(:, ind) = bsxfun(@times, L(:, ind), 1./sqrt(D(ind))');
        L(1:node_num+1:end) = 1;

        % The second normalized
        % Ln = D^(-1/2) L D^(-1/2);

        HCP_sMRI(sample_i).D = diag(D);
        HCP_sMRI(sample_i).L = L; % diag(D)-W
        HCP_sMRI(sample_i).symL = (L + L')/2;

        [Phi_temp,value,~] = svd(HCP_sMRI(sample_i).symL);
        [E,inds] = sort(diag(value),'ascend');
        Phi_temp = Phi_temp(:,inds);

        signs=sign(Phi_temp(1,:));
        signs(signs==0)=1;
        U = Phi_temp*diag(signs);

        HCP_sMRI(sample_i).Ne = nnz(W)/2;
        HCP_sMRI(sample_i).mu = max(abs(U(:)));
        HCP_sMRI(sample_i).Phi{1}=U;
        HCP_sMRI(sample_i).Eigenvalue{1}=E;
        HCP_sMRI(sample_i).TransformedPhi{1}=HCP_sMRI(sample_i).Phi{1};

        % Design filter bank
        order=1000; % for filters % 10000
        lmax = max(HCP_sMRI(sample_i).Eigenvalue{1});
        if wavelet_sperate == 1
            shifted_ends = zeros(filter_bandnumber+2,1);
            shifted_ends(1,1) = abs(HCP_sMRI(sample_i).Eigenvalue{1}(1));
            for i=1:filter_bandnumber
                shifted_ends(i+1,1) = abs(HCP_sMRI(sample_i).Eigenvalue{1}(i));
            end
            shifted_ends(filter_bandnumber+2) = shifted_ends(filter_bandnumber+1);
        elseif wavelet_sperate == 0
            shifted_ends = zeros(filter_bandnumber+1,1);
            for i=1:filter_bandnumber
                shifted_ends(i+1,1) = lmax*(i/filter_bandnumber)^2;
            end
        else 
            shifted_ends = zeros(5,1); % 1-3, 330-331, 334, 345-346, 354
            shifted_ends(1,1) = 0; % +: 1, 2, 3, 18, 334, 354; -: 330, 331, 337, 345, 346.
            shifted_ends(2,1) = lmax*(15/40)^2;
            shifted_ends(3,1) = lmax*(25/40)^2;
            shifted_ends(4,1) = lmax*(33/40)^2;
            shifted_ends(5,1) = lmax*(40/40)^2;
        end

        eigenvalue_filter = cell(2,filter_bandnumber);
        approx_filter_for_figure={};
        count_no_one = 0;
        peak_value_limited = 0.9;
        for filter=1:filter_bandnumber
            if normal_or_single
                temp_filter = zeros(360, 1);
                temp_filter(filter,1) = 1; % HCP_sMRI(sample_i).Eigenvalue{1}(filter,1);
                eigenvalue_filter{sample_i,filter}=temp_filter;
                subdict{sample_i,filter}=HCP_sMRI(sample_i).Phi{1}*diag(temp_filter)*HCP_sMRI(sample_i).Phi{1}';
            else
                if wavelet_sperate == 1
                    % reduce_range = (shifted_ends(filter+2) - shifted_ends(filter));
                    [cheby_coeffs,jack_cheby_coeffs]=gsp_jackson_cheby_coeff(shifted_ends(filter), shifted_ends(filter+2),[0 lmax], order);
                elseif wavelet_sperate == 0
                    [cheby_coeffs,jack_cheby_coeffs]=gsp_jackson_cheby_coeff(shifted_ends(filter), shifted_ends(filter+1),[0 lmax], order);
                else
                    if filter<5
                        [cheby_coeffs,jack_cheby_coeffs]=gsp_jackson_cheby_coeff(shifted_ends(filter), shifted_ends(filter+1),[0 lmax], order);
                    end
                end
                approx_filter=@(x) gsp_cheby_eval(x,jack_cheby_coeffs,[0,lmax]);
                approx_filter_cheby=@(x) gsp_cheby_eval(x,cheby_coeffs,[0,lmax]);
                eigenvalue_filter{sample_i,filter} = approx_filter(HCP_sMRI(sample_i).Eigenvalue{1});
                subdict{sample_i,filter}=HCP_sMRI(sample_i).Phi{1}*diag(approx_filter(HCP_sMRI(sample_i).Eigenvalue{1}))*HCP_sMRI(sample_i).Phi{1}';
                % subdict{sample_i,filter}=subdict{sample_i,filter};
                subdict{sample_i,filter}=subdict{sample_i,filter}.*path;
                approx_filter_for_figure{filter}=approx_filter;
                if max(approx_filter(HCP_sMRI(sample_i).Eigenvalue{1}))<peak_value_limited
                    count_no_one = count_no_one + 1;
                    fprintf('The peak value of the %d-th sampling distirbution (wavelet filter) is %.4f!\n', filter, max(approx_filter(HCP_sMRI(sample_i).Eigenvalue{1})))
                end
            end
        end
        fprintf('There are %d sampling distirbutions (wavelet filters) whose peak is smaller than %.4f!\n', count_no_one, peak_value_limited)

        fprintf('Sample: %d\n',sample_i)
    end

    if normal_or_single
        subdict_filename = ['wavelet_',num2str(filter_bandnumber),'_single_original_pn.mat'];
        subdict_path = fullfile(current_path,'output','subdict',subdict_filename);
        save(subdict_path,'subdict','-v6')
        eigenvalue_filter_filename = ['eigenvalue_filter_',num2str(filter_bandnumber),'_single_original_pn.mat'];
        eigenvalue_filter_path = fullfile(current_path,'output','subdict',eigenvalue_filter_filename);
        save(eigenvalue_filter_path,'eigenvalue_filter','-v6')
    else
        subdict_filename = ['wavelet_',num2str(filter_bandnumber),'.mat'];
        subdict_path = fullfile(current_path,'output','subdict',subdict_filename);
        save(subdict_path,'subdict','-v6')
        eigenvalue_filter_filename = ['eigenvalue_filter_',num2str(filter_bandnumber),'.mat'];
        eigenvalue_filter_path = fullfile(current_path,'output','subdict',eigenvalue_filter_filename);
        save(eigenvalue_filter_path,'eigenvalue_filter','-v6')

        figure;
        param_filter.line_width=5;
        param_filter.npoints = 1000;

        lambdas = linspace(0,lmax,param_filter.npoints);
        Nf=numel(approx_filter_for_figure);
        fd=zeros(length(lambdas),Nf);
        for k=1:Nf
            fd(:,k)=approx_filter_for_figure{k}(lambdas);
        end

        plot(lambdas,fd,'LineWidth',param_filter.line_width);

        % lgd=legend({'$j=1$','$j=2$','$j=3$','$j=4$','$j=5$','$j=6$','$j=7$','$j=8$','$j=9$','$j=10$'},'Interpreter','latex','Location','northoutside');
        % legend('boxoff')
        % lgd.NumColumns = 3;
        set(gca, 'Fontname', 'Times New Roman','FontSize',24,'FontWeight','bold');
        xlabel('$\lambda$','Interpreter','LaTex','FontSize',24,'FontWeight','bold')
        box off;

        filter_figure_filename = ['approx_filter_',num2str(filter_bandnumber),'.png'];
        filter_figure_path = fullfile(current_path,'output','subdict',filter_figure_filename);
        exportgraphics(gcf,filter_figure_path,'ContentType','vector')
    end

    HCP_sMRI_path = fullfile(current_path,'output','subdict','HCP_sMRI.mat');
    save(HCP_sMRI_path,'HCP_sMRI','-v6')
end