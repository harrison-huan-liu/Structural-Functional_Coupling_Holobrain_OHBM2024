close all;
clear all;
clc;

rand('seed',0);
randn('seed',0);

subdict=cell(135,1);

for sample_i=1:135
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

    for lambda = 1:90
        approx_filter = zeros(90,1);
        approx_filter(lambda,1) = 1;
        subdict{sample_i,lambda}=G.U*diag(approx_filter)*G.U';
    end

    fprintf('Sample: %d\n',sample_i)
end

save('subdict.mat','subdict','-v6')
