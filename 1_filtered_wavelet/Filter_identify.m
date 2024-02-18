close all;
clear all;
clc;

rand('seed',0);
randn('seed',0);

subdict1=cell(135,1);
subdict2=cell(135,1);
subdict3=cell(135,1);

for sample_i=1:1
    % Graph
    G=gsp_bunny(sample_i);
    % G.N=90;
    % G.W=;
    % G.coords=;
    % G.A=sparse(G.W>0);
    G=gsp_compute_fourier_basis(G);

    % Parameters we do not usually change 
    param.order=50; % for density estimation
    param.num_vec=30; % number of random vectors used to estimate spectral density and non-uniform sampling distributions

    % Design filter bank
    order=100; %100; % for filters
    param.spacing = 0; % 0-even number of eigenvalue in each band; 1-logarithmic  
    num_bands=9;
    param_filter.plot_eigenvalues=0;
    param_filter.line_width=5;
    param_filter.show_sum=0;
    plot_param.show_edges=0;
    plot_param.vertex_size=100;

    % plot_param.climits=[-0.2,0.2];%[-0.025,0.025];

    [filter_bank2,shifted_ends2,band_ends2,G] = mcsfb_design_filter_bank(G,num_bands,param);

    filter1=3;
    [cheby_coeffs1,jack_cheby_coeffs1]=gsp_jackson_cheby_coeff(shifted_ends2(filter1), shifted_ends2(filter1+1),[0 G.lmax], order);
    approx_filter1=@(x) gsp_cheby_eval(x,jack_cheby_coeffs1,[0,G.lmax]);
    approx_filter1_cheby=@(x) gsp_cheby_eval(x,cheby_coeffs1,[0,G.lmax]);

    filter2=6;
    [cheby_coeffs2,jack_cheby_coeffs2]=gsp_jackson_cheby_coeff(shifted_ends2(filter2), shifted_ends2(filter2+1),[0 G.lmax], order);
    approx_filter2=@(x) gsp_cheby_eval(x,jack_cheby_coeffs2,[0,G.lmax]);
    approx_filter2_cheby=@(x) gsp_cheby_eval(x,cheby_coeffs2,[0,G.lmax]);

    filter3=9;
    [cheby_coeffs3,jack_cheby_coeffs3]=gsp_jackson_cheby_coeff(shifted_ends2(filter3), shifted_ends2(filter3+1)-1,[0 G.lmax], order);
    approx_filter3=@(x) gsp_cheby_eval(x,jack_cheby_coeffs3,[0,G.lmax]);
    approx_filter3_cheby=@(x) gsp_cheby_eval(x,cheby_coeffs3,[0,G.lmax]);

    subdict1{sample_i,1}=G.U*diag(approx_filter1(G.e))*G.U';

    subdict2{sample_i,1}=G.U*diag(approx_filter2(G.e))*G.U';

    subdict3{sample_i,1}=G.U*diag(approx_filter3(G.e))*G.U';
    fprintf('Sample: %d\n',sample_i)
%subdict2=G.U*diag(filter_bank{filter2}(G.e))*G.U';
end
% max_val=max(abs(atom3));
% plot_param.climits = [-max_val,max_val];

wavelet1 = subdict1{1,1};
wavelet2 = subdict2{1,1};
wavelet3 = subdict3{1,1};
save('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\HCP_dynamics\wavelet1.mat','wavelet1','-v6')
save('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\HCP_dynamics\wavelet2.mat','wavelet2','-v6')
save('D:\research\Nonliner Dimensional Reduction\5_IPMI_to_Journal\HCP_dynamics\wavelet3.mat','wavelet3','-v6')

approx_filter={approx_filter1;approx_filter2;approx_filter3};
figure;
gsp_plot_filter(G,approx_filter,param_filter);
%gsp_plot_filter(G,filter_bank{filter2},param_filter);
set(gca, 'Fontname', 'Times New Roman','FontSize',24,'FontWeight','bold');
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24,'FontWeight','bold')
box off;

figure
subdict1(logical(eye(size(subdict1))))=0;
subdict2(logical(eye(size(subdict2))))=0;
subdict3(logical(eye(size(subdict3))))=0;

c_min = min([subdict1(:);subdict2(:);subdict3(:)]);
c_max = max([subdict1(:);subdict2(:);subdict3(:)]);

tiledlayout(2,2,'TileSpacing','Compact');
nexttile
imagesc(subdict1,[c_min c_max])
set(gca, 'Fontname', 'Times New Roman','FontSize',20,'FontWeight','bold');
nexttile
imagesc(subdict2,[c_min c_max])
set(gca, 'Fontname', 'Times New Roman','FontSize',20,'FontWeight','bold');
nexttile
imagesc(subdict3,[c_min c_max])
set(gca, 'Fontname', 'Times New Roman','FontSize',20,'FontWeight','bold');

cb = colorbar;
cb.Layout.Tile = 'east';
% colormap(jet)
% max_sub=max(max(abs(subdict2)));
% caxis([-max_sub max_sub])

save('subdict1.mat','subdict1','-v6')
save('subdict2.mat','subdict2','-v6')
save('subdict3.mat','subdict3','-v6')

cb = colorbar;
cb.Layout.Tile = 'east';

num_filters=3;
param.order = 80; % used for density estimation and analysis filtering
param.replacement=0;
param.grid_order=100; %1000;
param.num_vec=30;

% Downsampling sets
signal=G.BOLDGraph;
param.cdf_method='kpm';
G=spectral_cdf_approx2(G, param);

% for each band, non-uniform adapted sampling with number of measurements
% adjusted for signal energy

% filter_type='wav_itersine';
% filts=spm_parseval_filters(G,num_filters,filter_type);
filts{1,1}=approx_filter1;
filts{2,1}=approx_filter2;
filts{3,1}=approx_filter3;
anal_coeffs=gsp_filter_analysis(G,filts,signal,param);
anal_coeffs_mat=reshape(anal_coeffs,G.N,num_filters);
param.jackson=0;
[param.signal_projections,filter_coeffs]=mcsfb_apply_filters(G,signal,filts,param);
r=gsp_cheby_opX(G,filter_coeffs);
nb_meas=zeros(num_filters,1);
num_vec=size(G.X,2);
weights=cell(num_filters,1);
for i=1:num_filters
    ri=r(:,(i-1)*num_vec+1:i*num_vec);
    nb_meas(i)=round(gsp_hutch(G,ri));
    norm_Uk= sum(ri.^2, 2);
    weights{i}=norm_Uk/sum(norm_Uk);
    if i>1
        weights{i}=norm_Uk.*log(1+abs(param.signal_projections(:,i)));
        weights{i}=weights{i}/sum(weights{i});
    end
end
proj_norms=sqrt(sum(param.signal_projections.^2));
nb_meas=nb_meas.*(log(1+proj_norms'));
% nb_meas(1)=nb_meas(1)*8; % extra for first band just for clarity of demo
total_samples=sum(nb_meas);
nb_meas=round(G.N/total_samples*nb_meas);
total_samples=sum(nb_meas);
if total_samples>G.N % eliminate from last band
    extra=total_samples-G.N;
    nb_meas(num_filters)=nb_meas(num_filters)-extra;
elseif total_samples<G.N % add more to first band
    nb_meas(1)=nb_meas(1)+G.N-total_samples;
end

centers=cell(num_filters,1);
for i = 1:num_filters % perform the sampling
    [~, selected] = build_sampling_matrix(G, weights{i}, nb_meas(i), param);
    centers{i} = selected;
end

figure
subdict = [subdict1(:,centers{1}) subdict2(:,centers{2}) subdict3(:,centers{3})];
imagesc(subdict)
set(gca, 'Fontname', 'Times New Roman','FontSize',20,'FontWeight','bold');
colormap(jet)
max_sub=max(max(abs(subdict)));
caxis([-max_sub max_sub])
colorbar;
save('subdict.mat','subdict','-v6')


plot_param.show_edges=0;
for i=1:num_filters
    Wi=G.W(centers{i},centers{i});
    Gi=gsp_graph(Wi,G.coords(centers{i},:),G.plotting.limits);
    figure;
    gsp_plot_signal(Gi,anal_coeffs_mat(centers{i},i),plot_param);
    view(0,90);
    set(gca,'FontSize',24);
    title(sprintf('coefficients {\\boldmath$\\alpha_{%d}$}',i),'interpreter','latex');
    colorbar off;
end
