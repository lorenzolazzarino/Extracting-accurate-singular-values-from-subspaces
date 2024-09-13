%% LANbound SCRIPT
% The aim of this script is provide the code to obtain the experiments used
% in the paper "Matrix Perturbation Analysis of Methods for Extracting 
% Singular Values from Approximate Singular Subspaces" 
% by L. Lazzarino, H. Al Daas, Y. Nakatsukasa.
%
% Throughout the script, together with the paper's notations,
% the following notation are used:
%
% *_exp = * is a variable in an experiment where singular values are
%           exponentially decaying
% *_alg = * is a variable in an experiment where singular values are
%           algebraically decaying
% *_*_ov = *_* is a variable in an experiment with oversample (if the 
%              variable name do not contain _ov then it is a variable in 
%              an experiment without oversample
%
% The script is divided in the following sections:
% - PROBLEM GENERATION: Setting the parameter and creating the matrix A
% - EXPERIMENT 1: Generalized Nystrom without Oversampling
% - EXPERIMENT 2: Generalized Nystrom with Oversampling
% - EXPERIMENT 3: Methods Comparison
% - EXPERIMENT 4: Computability
% - FUNCTIONS: (1) Generate Approximated Subspaces, (2) Generalized Nystrom
%              (3) Rayleigh-Ritz Method, (4) HMT Method
%
% Note: The function LANbound is needed to run this script.


% ------------------------- Set Plot parameters
    MS = 'Markersize'; LW = 'linewidth'; FS = 'fontsize';  IN = 'Interpreter';
    ms = 10; ms1=6; fs = 16; fs1 = 10; lw = 2; lw1 = 1; in = 'latex';
    
    % If save_plot is set to be 1, each figure below will be automatically
    % saved as a jpeg in a folder named Figures, THAT HAS TO EXISTS IN 
    % CURRENT PATH. 
    save_plot = 0;
    % ATTENTION: If a file with the below filenames already exists then 
    % it will be overwritten!!! To avoid it, change the save_number
    % parameter, that will be added at the end of the filenames so to avoid
    % to loose previously saved plots!
    save_number = 1;
%% PROBLEM GENERATION
    
    m = 1000;   % Row size of matrix A
    n = 1000;   % Column size of matrix A
    
% ------------------------- Construct Exact (left/right) Singular Vectors
    
    [Uex,~] = qr(randn(m));        % Haar matrix
    [Vex,~] = qr(randn(n));        % Haar matrix
     
    
% ------------------------- Construct Exact Singular Values
    
    % Algebraic decay:
    D_alg = (1./(1:n)).^4; 
    D_alg = D_alg(:); 
    D_alg = sort(D_alg,'descend'); 
    
    Dex_alg = [diag(D_alg); zeros(m-n,n)];  % Exact singular values matrix
    
    A_alg = Uex*Dex_alg*Vex';               % Form A from its exact SVD
    
    % Exponential decay:
    D_exp = logspace(-30,0,n)'; 
    D_exp = sort(D_exp,'descend'); 
    
    
    Dex_exp = [diag(D_exp); zeros(m-n,n)];  % Exact singular values matrix
    A_exp = Uex*Dex_exp*Vex';               % Form A from its exact SVD
    


%% EXPERIMENT 1: GENERALIZED NYSTROM WITHOUT OVERSAMPLING
% In this experiment we consider two approximations tV and tU of the leading
% singular subspaces of the same column size, compute the Generalized
% Nystrom approximation and the correspondent bound. We do it for both 
% decays.

    r = 200;    % Column Size of approximate singular subspaces
  
% ------------------------- Generate Approximate Singular Subspaces

    % See function genSubsp below - These approximations will be used for
    % all experiments without oversampling below
    [tV_exp,tU_exp] = genSubsp(A_exp,r);
    [tV_alg,tU_alg] = genSubsp(A_alg,r);

% ------------------------- Exponential decay: Compute GN and Bound    
   
    svGN_exp = GN(A_exp,tV_exp,tU_exp);
    [BoundGN_exp, WeylGN_exp] = LANbound(A_exp,D_exp,tV_exp,tU_exp);
    
% ------------------------- Algebraic decay: Compute GN and Bound    
    
    svGN_alg = GN(A_alg,tV_alg,tU_alg);
    [BoundGN_alg, WeylGN_alg] = LANbound(A_alg,D_alg,tV_alg,tU_alg);
    
% ------------------------- Plot Results
    figure()    % Exponential decay
        semilogy(abs(D_exp(1:r)-svGN_exp(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_exp,'b-',LW,lw)
        semilogy(BoundGN_exp,'g-',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl',...
                'Bound',FS,fs,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('GN, exponential sing. decay, no-oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
            filename = sprintf('Figures/GN_exp_%d.jpg',save_number);
            % Save subplot as JPEG
            saveas(gca, filename, 'jpeg');
        end

    figure()    % Algebraic decay
        semilogy(abs(D_alg(1:r)-svGN_alg(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_alg,'b-',LW,lw)
        semilogy(BoundGN_alg,'g-',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl',...
                'Bound',FS,fs,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('GN, algebraic sing. decay, no-oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/GN_alg_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
%% EXPERIMENT 2: GENERALIZED NYSTROM WITH OVERSAMPLING
% In this experiment we consider two approximations tV and tU of the leading
% singular subspaces. tV is a n x r matrix, while tU is a m x (r+l). We 
% compute the Generalized Nystrom approximation and the correspondent 
% bound. We do it for both decays. We eventually consider an
% heuristic strategy to improve the sharpness of the bound.
    
    r = 200;    % Column Size of tV
    rl = 300;   % rl:= r+l = 1.5r = Column Size of tU
    
    
% ------------------------- Generate Approximate Singular Subspaces

    % See function genSubsp below - These approximations will be used for
    % all experiments with oversampling below
    [tV_exp_ov,tU_exp_ov] = genSubsp(A_exp,r,rl);
    [tV_alg_ov,tU_alg_ov] = genSubsp(A_alg,r,rl);
    
% ------------------------- Exponential decay: Compute GN and Bound 

    svGN_exp_ov = GN(A_exp,tV_exp_ov,tU_exp_ov);
    [BoundGN_exp_ov, WeylGN_exp_ov] = LANbound(A_exp,D_exp,tV_exp_ov,tU_exp_ov);
    
% ------------------------- Algebraic decay: Compute GN and Bound 

    svGN_alg_ov = GN(A_alg,tV_alg_ov,tU_alg_ov);
    [BoundGN_alg_ov, WeylGN_alg_ov] = LANbound(A_alg,D_alg,tV_alg_ov,tU_alg_ov);
    
% ------------------------- Plot Results

    figure()  % Exponential decay
        semilogy(abs(D_exp(1:r)-svGN_exp_ov(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_exp_ov,'b-',LW,lw)
        semilogy(BoundGN_exp_ov,'g-',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Bound',...
                FS,fs,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('GN, exponential sing. decay, oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/GN_exp_ov_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end

    figure()    % Algebraic decay
        semilogy(abs(D_alg(1:r)-svGN_alg_ov(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_alg_ov,'b-',LW,lw)
        semilogy(BoundGN_alg_ov,'g-',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Bound',...
                FS,fs,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('GN, algebraic sing. decay, oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/GN_alg_ov_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
    
%------------------------- Add Improved Bound
    
    % Exponential decay: Improved Bound
    ImprBoundGN_exp_ov = LANbound(A_exp,D_exp,tV_exp_ov,tU_exp_ov,'ImproveBound','Y');
    
    % Algebraic decay: Improved Bound
    ImprBoundGN_alg_ov = LANbound(A_alg,D_alg,tV_alg_ov,tU_alg_ov,'ImproveBound','Y');
    
% ------------------------- Plot Results - Improved Bound

    figure()    % Exponential decay
        semilogy(abs(D_exp(1:r)-svGN_exp_ov(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_exp_ov,'b-',LW,lw)
        semilogy(BoundGN_exp_ov,'g-',MS,ms,LW,lw)
        semilogy(ImprBoundGN_exp_ov,'m--',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Bound',...
                'Improved Bound',FS,fs,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('GN, exponential sing. decay, oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/GN_exp_ov_impr_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end

    figure()    % Algebraic decay
        semilogy(abs(D_alg(1:r)-svGN_alg_ov(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_alg_ov,'b-',LW,lw)
        semilogy(BoundGN_alg_ov,'g-',MS,ms,LW,lw)
        semilogy(ImprBoundGN_alg_ov,'m--',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Bound',...
                'Improved Bound',FS,fs,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('GN, algebraic sing. decay, oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/GN_alg_ov_impr_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
    
%% EXPERIMENT 3: METHOD COMPARISON
% In this experiment, given the approximate singular subspaces computed above,
% we extract singular values using: Generalized Nystrom (we consider
% results obtained above), Rayleigh-Ritz method, HMT method, and SVD 
% method. We do it for both decays, and with and without oversample.

% EXPONENTIAL DECAY ------------------------------------------------

    % Rayleigh-Ritz without oversampling
    svRR_exp = RR(A_exp,tV_exp,tU_exp);
    [BoundRR_exp, WeylRR_exp] = LANbound(A_exp,D_exp,tV_exp,tU_exp,'Method','RR');

    % Rayleigh-Ritz with oversampling
    svRR_exp_ov = RR(A_exp,tV_exp_ov,tU_exp_ov);
    [BoundRR_exp_ov, WeylRR_exp_ov] = LANbound(A_exp,D_exp,tV_exp_ov,tU_exp_ov,'Method','RR');
    
    % HMT method
    svHMT_exp = HMT(A_exp,tV_exp);
    [BoundHMT_exp, WeylHMT_exp] = LANbound(A_exp,D_exp,tV_exp);
    
    % SVD method
    svSVD_exp = svd(A_exp*tV_exp);
    [BoundSVD_exp, WeylSVD_exp] = LANbound(A_exp,D_exp,tV_exp,[],'Method','SVD');
    
% ------------------------- Plot Results

    figure()    % No oversample for every method
        % GN
        semilogy(abs(D_exp(1:r)-svGN_exp(1:r)),'ro',MS,ms1)
        hold on
        grid on
        yline(WeylGN_exp,'r--',LW,lw1)
        semilogy(BoundGN_exp,'r--',MS,ms1,LW,lw1)
        
        % HMT
        semilogy(abs(D_exp(1:r)-svHMT_exp(1:r)),'b+',MS,ms1)
        yline(WeylHMT_exp,'b-.',LW,lw1)
        semilogy(BoundHMT_exp,'b-.',MS,ms1,LW,lw1)
        
        % RR
        semilogy(abs(D_exp(1:r)-svRR_exp(1:r)),'gs',MS,ms1)
        yline(WeylRR_exp,'g-',LW,lw1)
        semilogy(BoundRR_exp,'g-',MS,ms1,LW,lw1)
        
        % SVD
        semilogy(abs(D_exp(1:r)-svSVD_exp(1:r)),'k.',MS,ms1)
        yline(WeylSVD_exp,'k:',LW,lw1)
        semilogy(BoundSVD_exp,'k:',MS,ms1,LW,lw1)
        
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl for GN','Bound for GN',... 
            '$|\sigma_i(A) - \sigma_i(A_{HMT})|$','Weyl for HMT','Bound for HMT',...
            '$|\sigma_i(A) - \sigma_i(A_{RR})|$','Weyl for RR','Bound for RR', ...
            '$|\sigma_i(A) - \sigma_i(A_{SVD})|$','Weyl for SVD','Bound for SVD',...
            FS,fs1,IN,in,'Location','Best','NumColumns',2)

        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Comparison, exponential sing. distr., no-oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Comparison_exp_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
  %      
    figure()  % Oversample where meaningful
        % GN
        semilogy(abs(D_exp(1:r)-svGN_exp_ov(1:r)),'ro',MS,ms1)
        hold on
        grid on
        yline(WeylGN_exp_ov,'r--',LW,lw1)
        semilogy(BoundGN_exp_ov,'r--',MS,ms1,LW,lw1)
        semilogy(ImprBoundGN_exp_ov,'r--',MS,ms1,LW,lw)
        
        % HMT
        semilogy(abs(D_exp(1:r)-svHMT_exp(1:r)),'b+',MS,ms1)
        yline(WeylHMT_exp,'b-.',LW,lw1)
        semilogy(BoundHMT_exp,'b-.',MS,ms1,LW,lw1)
        
        % RR
        semilogy(abs(D_exp(1:r)-svRR_exp_ov(1:r)),'gs',MS,ms1)
        yline(WeylRR_exp_ov,'g-',LW,lw1)
        semilogy(BoundRR_exp_ov,'g-',MS,ms1,LW,lw1)
        
        % SVD
        semilogy(abs(D_exp(1:r)-svSVD_exp(1:r)),'k.',MS,ms1)
        yline(WeylSVD_exp,'k:',LW,lw1)
        semilogy(BoundSVD_exp,'k:',MS,ms1,LW,lw1)
        
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl for GN','Bound for GN','Improved Bound for GN',... 
            '$|\sigma_i(A) - \sigma_i(A_{HMT})|$','Weyl for HMT','Bound for HMT',...
            '$|\sigma_i(A) - \sigma_i(A_{RR})|$','Weyl for RR','Bound for RR', ...
            '$|\sigma_i(A) - \sigma_i(A_{SVD})|$','Weyl for SVD','Bound for SVD',...
             FS,fs1,IN,in,'Location','Northwest','NumColumns',2)

        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Comparison, exponential sing. distr., oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Comparison_exp_ov_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
        
% ALGEBRAIC DECAY --------------------------------------------------

    % Rayleigh-Ritz without oversampling
    svRR_alg = RR(A_alg,tV_alg,tU_alg);
    [BoundRR_alg, WeylRR_alg] = LANbound(A_alg,D_alg,tV_alg,tU_alg,'Method','RR');
    
    % Rayleigh-Ritz with oversampling
    svRR_alg_ov = RR(A_alg,tV_alg_ov,tU_alg_ov);
    [BoundRR_alg_ov, WeylRR_alg_ov] = LANbound(A_alg,D_alg,tV_alg_ov,tU_alg_ov,'Method','RR');
    
    % HMT method
    svHMT_alg = HMT(A_alg,tV_alg);
    [BoundHMT_alg, WeylHMT_alg] = LANbound(A_alg,D_alg,tV_alg);
    
    % SVD method
    svSVD_alg = svd(A_alg*tV_alg);
    [BoundSVD_alg, WeylSVD_alg] = LANbound(A_alg,D_alg,tV_alg,[],'Method','SVD');
    
% ------------------------- Plot Results

    figure()    % No oversample for every method
        % GN
        semilogy(abs(D_alg(1:r)-svGN_alg(1:r)),'ro',MS,ms1)
        hold on
        grid on
        yline(WeylGN_alg,'r--',LW,lw1)
        semilogy(BoundGN_alg,'r--',MS,ms1,LW,lw1)
        
        % HMT
        semilogy(abs(D_alg(1:r)-svHMT_alg(1:r)),'b+',MS,ms1)
        yline(WeylHMT_alg,'b-.',LW,lw1)
        semilogy(BoundHMT_alg,'b-.',MS,ms1,LW,lw1)
        
        % RR
        semilogy(abs(D_alg(1:r)-svRR_alg(1:r)),'gs',MS,ms1)
        yline(WeylRR_alg,'g-',LW,lw1)
        semilogy(BoundRR_alg,'g-',MS,ms1,LW,lw1)
        
        % SVD
        semilogy(abs(D_alg(1:r)-svSVD_alg(1:r)),'k.',MS,ms1)
        yline(WeylSVD_alg,'k:',LW,lw1)
        semilogy(BoundSVD_alg,'k:',MS,ms1,LW,lw1)
        
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl for GN','Bound for GN',... 
            '$|\sigma_i(A) - \sigma_i(A_{HMT})|$','Weyl for HMT','Bound for HMT',...
            '$|\sigma_i(A) - \sigma_i(A_{RR})|$','Weyl for RR','Bound for RR', ...
            '$|\sigma_i(A) - \sigma_i(A_{SVD})|$','Weyl for SVD','Bound for SVD',...
            FS,fs1,IN,in,'Location','Best','NumColumns',2)

        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Comparison, algebraic sing. distr., no-oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Comparison_alg_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
        
    figure()  % Oversample where meaningful
        % GN
        semilogy(abs(D_alg(1:r)-svGN_alg_ov(1:r)),'ro',MS,ms1)
        hold on
        grid on
        yline(WeylGN_alg_ov,'r--',LW,lw1)
        semilogy(BoundGN_alg_ov,'r--',MS,ms1,LW,lw1)
        semilogy(ImprBoundGN_alg_ov,'r--',MS,ms1,LW,lw);
        
        
        % HMT
        semilogy(abs(D_alg(1:r)-svHMT_alg(1:r)),'b+',MS,ms1)
        yline(WeylHMT_alg,'b-.',LW,lw1)
        semilogy(BoundHMT_alg,'b-.',MS,ms1,LW,lw1)
        
        
        % RR
        semilogy(abs(D_alg(1:r)-svRR_alg_ov(1:r)),'gs',MS,ms1)
        yline(WeylRR_alg_ov,'g-',LW,lw1)
        semilogy(BoundRR_alg_ov,'g-',MS,ms1,LW,lw1)
        
        % SVD
        semilogy(abs(D_alg(1:r)-svSVD_alg(1:r)),'k.',MS,ms1)
        yline(WeylSVD_alg,'k:',LW,lw1)
        semilogy(BoundSVD_alg,'k:',MS,ms1,LW,lw1)
        
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl for GN','Bound for GN','Improved Bound for GN',... 
            '$|\sigma_i(A) - \sigma_i(A_{HMT})|$','Weyl for HMT','Bound for HMT',...
            '$|\sigma_i(A) - \sigma_i(A_{RR})|$','Weyl for RR','Bound for RR', ...
            '$|\sigma_i(A) - \sigma_i(A_{SVD})|$','Weyl for SVD','Bound for SVD',...
            FS,fs1,IN,in,'Location','Best','NumColumns',2)

        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Comparison, algebraic sing. distr., oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Comparison_alg_ov_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end

%% EXPERIMENT 4: COMPUTABILITY
% In this experiment, we explore the computability of the bound by comparing
% the Backward bounds for the Generalized Nystrom Approximation with the
% forward bounds computed above. We do it for both decays and
% with/without oversampling.

% EXPONENTIAL DECAY ------------------------------------------------


    % Backward Bound - Case without oversample case
    Bound_exp_bb = LANbound(A_exp,svGN_exp,tV_exp,tU_exp,'Orientation','Backward');

    % Backward Bound - Case with oversample 
    Bound_exp_ov_bb = LANbound(A_exp,svGN_exp_ov,tV_exp_ov,tU_exp_ov,'Orientation','Backward');

    % Improved Backward Bound - Case with oversample 
    ImprBound_exp_ov_bb = LANbound(A_exp,svGN_exp_ov,tV_exp_ov,tU_exp_ov,'Orientation','Backward','ImproveBound','Y');
    
    % Approximated Backward Bound - Case without oversample case
    Bound_exp_Abb = LANbound(A_exp,svGN_exp,tV_exp,tU_exp,'Orientation','ApproxBackward');

    % Approximated Backward Bound - Case with oversample 
    Bound_exp_ov_Abb = LANbound(A_exp,svGN_exp_ov,tV_exp_ov,tU_exp_ov,'Orientation','ApproxBackward');

    % Improved Approximated Backward Bound - Case with oversample 
    ImprBound_exp_ov_Abb = LANbound(A_exp,svGN_exp_ov,tV_exp_ov,tU_exp_ov,'Orientation','ApproxBackward','ImproveBound','Y');

%% ------------------------- Plot Results
    figure()    % Without Oversample
        semilogy(abs(D_exp(1:r)-svGN_exp(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_exp,'b-',LW,lw)
        semilogy(BoundGN_exp,'g-',MS,ms,LW,lw)
        semilogy(Bound_exp_bb,'k:',MS,ms,LW,lw)
        semilogy(Bound_exp_Abb,'m--',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Forward Bound',...
                'Backward Bound','Approximated Backward Bound',...
                FS,13,IN,in,'Location','Best')

        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Computability, exponential sing. distr., no-oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Computability_exp_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
        
    figure()    % With Oversample
        semilogy(abs(D_exp(1:r)-svGN_exp_ov(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_exp_ov,'b-',LW,lw)
        semilogy(BoundGN_exp_ov,'g-',MS,ms,LW,lw)
        semilogy(Bound_exp_ov_bb,'k:',MS,ms,LW,lw)
        semilogy(Bound_exp_ov_Abb,'m--',MS,ms,LW,lw)
        semilogy(ImprBoundGN_exp_ov,'b-.',MS,ms,LW,lw)
        semilogy(ImprBound_exp_ov_bb,'r:',MS,ms,LW,lw)
        semilogy(ImprBound_exp_ov_Abb,'m--',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Forward Bound',...
            'Backward Bound','Approximated Backward Bound',...
            'Improved Forward Bound','Improved Backward Bound',...
            'Improved Approximated Backward Bound',...
            FS,fs1,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Computability, exponential sing. distr., oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Computability_exp_ov_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end

% ALGEBRAIC DECAY ------------------------------------------------


    % Backward Bound - Case without oversample case
    Bound_alg_bb = LANbound(A_alg,svGN_alg,tV_alg,tU_alg,'Orientation','Backward');
    
    % Backward Bound - Case with oversample 
    Bound_alg_ov_bb = LANbound(A_alg,svGN_alg_ov,tV_alg_ov,tU_alg_ov,'Orientation','Backward');
    
    % Improved Backward Bound - Case with oversample 
    ImprBound_alg_ov_bb = LANbound(A_alg,svGN_alg_ov,tV_alg_ov,tU_alg_ov,'Orientation','Backward','ImproveBound','Y');
    
    % Approximated Backward Bound - Case without oversample case
    Bound_alg_Abb = LANbound(A_alg,svGN_alg,tV_alg,tU_alg,'Orientation','ApproxBackward');
    
    % Approximated Backward Bound - Case with oversample 
    Bound_alg_ov_Abb = LANbound(A_alg,svGN_alg_ov,tV_alg_ov,tU_alg_ov,'Orientation','ApproxBackward');
    
    % Improved Approximated Backward Bound - Case with oversample 
    ImprBound_alg_ov_Abb = LANbound(A_alg,svGN_alg_ov,tV_alg_ov,tU_alg_ov,'Orientation','ApproxBackward','ImproveBound','Y');
    
% ------------------------- Plot Results

    figure()  % Without Oversample
        semilogy(abs(D_alg(1:r)-svGN_alg(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_alg,'b-',LW,lw)
        semilogy(BoundGN_alg,'g-',MS,ms,LW,lw)
        semilogy(Bound_alg_bb,'k:',MS,ms,LW,lw)
        semilogy(Bound_alg_Abb,'m--',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Forward Bound',...
                'Backward Bound','Approximated Backward Bound',...
                FS,fs,IN,in,'Location','Best')
        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Computability, algebraic sing. distr., no-oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Computability_alg_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end
        
    figure()  % With Oversample
        semilogy(abs(D_alg(1:r)-svGN_alg_ov(1:r)),'r.',MS,ms)
        hold on
        grid on
        yline(WeylGN_alg_ov,'b-',LW,lw)
        semilogy(BoundGN_alg_ov,'g-',MS,ms,LW,lw)
        semilogy(Bound_alg_ov_bb,'k:',MS,ms,LW,lw)
        semilogy(Bound_alg_ov_Abb,'m--',MS,ms,LW,lw)
        semilogy(ImprBoundGN_alg_ov,'b-.',MS,ms,LW,lw)
        semilogy(ImprBound_alg_ov_bb,'r:',MS,ms,LW,lw)
        semilogy(ImprBound_alg_ov_Abb,'m--',MS,ms,LW,lw)
        legend('$|\sigma_i(A) - \sigma_i(A_{GN})|$','Weyl','Forward Bound',...
                'Backward Bound','Approximated Backward Bound',...
                'Improved Forward Bound','Improved Backward Bound',...
                'Improved Approximated Backward Bound',...
                FS,fs1,IN,in,'Location','Best')

        xlabel('$i$',FS,fs,IN,in); set(gca,FS,fs);
        title('Computability, algebraic sing. distr., oversample',FS,fs)

        % Save subplot (if specified)
        if save_plot == 1
                filename = sprintf('Figures/Computability_alg_ov_%d.jpg',save_number);
                % Save subplot as JPEG
                saveas(gca, filename, 'jpeg');
        end

%% FUNCTIONS
% Below all functions used in the script.    
    
%------------------------- Generate approximate singular subspaces
    function [tV,tU] = genSubsp(A,r,rl)
        
        if nargin < 2
            error('Not enough Input')
        elseif (nargin < 3)    
            rl = r;
        elseif (nargin == 3) && (nargout == 1)
            warning('Parameter rl not used')
        end
        
        [m,n] = size(A);
        
        % Consider Gaussian Matrices
        Omega1 = randn(m,r); % Corresponding to tV
        Omega2 = randn(n,rl); % Corresponding to tU

        % Perform one power iteration (more power iteration could be used
        % to obtain even better approximations)
        AOmega2 = A*Omega2;
        Omega1A = Omega1'*A;

        % Compute Approximate leading singular subspaces
        [tU,~] = qr(AOmega2,0);
        [tV,~] = qr(Omega1A',0); 
    
    end
    
%------------------------- Generalized Nystrom
    function [svGN,AGN] = GN(A,tV,tU,compute_pseudo)
    % compute_pseudo is a parameter to select how to compute the
    % pseudoinverse, it can be set to be: 0 (default) for QR factorization
    % of the matrix to invert, 1 for the use of lsqminnorm, 2 for the use
    % of the MatLab backslash.

        if (nargin < 3 || isempty(tU) || isempty(tV))
            error('Not Enough Inputs')
        elseif nargin<4
            compute_pseudo = 0; % Default value for compute_pseudo
        end
        
        AtV = A*tV;
        [~,R1] = qr(AtV,0);
        tUA = tU'*A;
        [~,R2] = qr(tUA',0);

        % Form Generalized Nystrom Matrix 
        if compute_pseudo == 0
            [Q,R] =qr(tU'*AtV,0);
            P = Q'*R2';
            AGN = (R1/R)*P;
        elseif compute_pseudo == 1
            AGN = R1*(lsqminnorm(tU'*AtV,R2'));
        elseif compute_pseudo == 2
            AGN = (R1/(tU'*AtV))*R2';
        end
        
        % Compute singular values
        svGN = svd(AGN);
    
    end
    
%------------------------- Rayleigh-Ritz Method
    function [svRR,ARR] = RR(A,tV,tU)
        
        if (nargin < 3 || isempty(tU) || isempty(tV))
            error('Not Enough Inputs')
        end
    
        % Form Rayleigh-Ritz Matrix 
        AtV = A*tV;
        ARR = tU'*AtV; 

        % Compute singular values
        svRR = svd(ARR);
    
    end
    
%------------------------- HMT Method
    function [svHMT,AHMT] = HMT(A,tV)
            
        if (nargin < 2 || isempty(tV))
            error('Not Enough Inputs')
        end
        
        % Form HMT Matrix 
        AtV = A*tV;
        [Q,~] = qr(AtV,0);
        AHMT = Q'*A;

        % Compute singular values
        svHMT = svd(AHMT);
    
    end
    
    
    

