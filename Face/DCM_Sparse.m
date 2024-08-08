function [ output ] = DCM_Sparse(DCM, X, Y, options, args)
% [ output ] = DCM_Sparse(DCM, X, Y, options, args)
% 
% Variational Bayesian inversion of a linear DCM with regression DCM, 
% including sparsity constrains on the connectivity parameters. The
% function implements the VB update equations derived in Frasle et al.
% (2018).
%
% This file is modified by Haifeng Wu. Based on the original sparse 
% variational inference, high-pass and low-pass filtering are added, 
% imaginary numbers are changed to real numbers (enhance computational 
% efficiency) , a cost function (i.e. an error function )are added, and 
% GLM estimates are embedded. (2023)
% 
%   Input:
%   	DCM             - model structure
%       X               - design matrix (predictors)
%       Y               - data
%       options         - estimation options
%       args            - arguments
%
%   Output:
%       output          - output structure
% 
%   Reference:
%       Frasle, S., Lomakina, E.I., Kasper, L., Manjaly Z.M., Leff, A., 
%       Pruessmann, K.P., Buhmann, J.M., Stephan, K.E., 2018. A generative 
%       model of whole-brain effective connectivity. NeuroImage 179, 505-529.
%       https://doi.org/10.1016/j.neuroimage.2018.05.058
%
% ----------------------------------------------------------------------
% 
% Authors: Stefan Fraessle (stefanf@biomed.ee.ethz.ch), Ekaterina I. Lomakina
%          Haifeng Wu (whf5469@gmail.com)
% 
% Copyright (C) 2016-2023 Translational Neuromodeling Unit
%                         Institute for Biomedical Engineering
%                         University of Zurich & ETH Zurich
%                         School of Electrical & Information Technology
%                         Yunnan Minzu University
%
% ----------------------------------------------------------------------


% number of regions
[ny,nr] = size(DCM.Y.y);

% precision limit
pr = 10^(-5);

% add confound regressor dimensions

    DCM.b(:,:,end+1)  = DCM.b(:,:,1);
    DCM.cb(:,end+1)   = zeros(1,size(DCM.cb,1));

% get the priors
[m0, l0, a0, b0] = DCM_get_prior_all(DCM, options);

D  =  size (m0, 2);     

% define a random order for the regions
ord_regions = 1:nr;

% results arrays (per iteration)
mN   = zeros(nr,D);
sN   = cell(nr,1);
aN   = zeros(nr,1);
bN   = zeros(nr,1);
z    = zeros(nr,D);
logF = zeros(1,nr);

if ( ~isfield(args,'diagnostics_logF') ), args.diagnostics_logF = 0; end

% whether high-pass filter or low-pass filter
if options.filter_flag == 'H'
    
    % high-pass filter
    Yf = HPFilter (Y, args.Freq) ;
    Xf = HPFilter (X, args.Freq) ;
    for k = 1 : nr
        YF{k} = Yf(:,k) ;
        XF{k} = Xf(:,:,k) ;
    end

elseif options.filter_flag == 'L'
    
    % low-pass filter
    r_dt     = DCM.Y.dt/DCM.U.dt;
    h_fft    = fft( DCM.h );
    u        = full(DCM.U.u);
    nu       = size (u, 2);
    u        = ifft(fft(u).* repmat(h_fft, 1, nu));
    DCM.U.X0 = ones(size(DCM.U.u,1),1);
    u        = [u, DCM.U.X0];
    u        = u(1 : r_dt : end, : );
    u_fft    = fft(u);
    h_fft    = fft( DCM.h(1 : r_dt : end, : ) );
    [~, idx] = tapas_rdcm_filter(fft(DCM.Y.y), u_fft, h_fft, ny, options);
    for k = 1 : nr
        idf   = 2 : round( sum( idx(:,k) )/2 ) ;
        YTmp  = Y(idx(:,k),  k) ;        
        YF{k} = [real( YTmp(idf) ); imag( YTmp(idf) )];
        XTmp  = X(idx(:,k),:,k);
        XF{k} = [real( XTmp(idf,:) ); imag( XTmp(idf,:) )];
    end

end

reverseStr = '';   
for k = ord_regions
  
    % output progress of regions
    if options.search.flag == 0
        msg = sprintf('Processing region: %d/%d', k, length(ord_regions));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    % results arrays (per region)
    mN_iter             = cell(options.iter,1);
    sN_iter             = cell(options.iter,1);
    aN_iter             = cell(options.iter,1);
    bN_iter             = cell(options.iter,1);
    z_iter              = cell(options.iter,1);
    l0_iter             = cell(options.iter,1);
    logF_iter           = NaN(options.iter,1);
    log_lik_iter        = NaN(options.iter,1);
    log_p_weight_iter	= NaN(options.iter,1);
    log_p_prec_iter     = NaN(options.iter,1);
    log_p_z_iter        = NaN(options.iter,1);
    log_q_weight_iter	= NaN(options.iter,1);
    log_q_prec_iter     = NaN(options.iter,1);
    log_q_z_iter        = NaN(options.iter,1);
    
    % iterate over the number of re-initializations
    for iter = 1:options.iter

        % filter unnecessary data
%         idx_y = ~isnan(Y(:,k));
        X_r   = XF{k};
        Y_r   = YF{k};

        % effective number of data points
%         N_eff = sum(idx_y);
        N_eff = size(Y_r,1);
        
        % get the priors per regions
        l0_r = diag(l0(k,:));
        m0_r = m0(k,:)';

        % set p0 (feature informativeness)
        try
            p0 = ones(D,1)*args.Pro;
        catch
            p0 = ones(D,1)*0.5;
        end
        
        % ensure self-connectivity
        p0(k) = 1.0;
        
        p0(end-(size(DCM.cb,2)-1):end-1) = args.Pro *DCM.cb(k,1:end-1);

        
        % results array per regions
        sN_r = zeros(size(l0_r));
        mN_r = zeros(size(m0_r));
        
        
        %% initalize
        
        % estimate variables X'X and X'Y per region
        W = X_r'*X_r;
        v = X_r'*Y_r;
        
        % intialize z, t & aN per regions
        z_r   = p0;
        t     = a0/b0;
        aN_r  = a0 + N_eff/(2*DCM.r_dt);
        
        
        % cycle stops after 500 iterations;
        count = 500;
        
        
        %% updating mN, sN

        % define random matrix Z
        Z = diag(z_r);
        
        % expectation (over Z) of ZX'XZ
        G = Z*W*Z; 
        G(eye(D)>0) = z_r.*diag(W);

        % set old F
        logF_old = -inf;
        
        % convergence criterion
        convergence = 0;

        % loop until convergence
        while ~convergence
            
            
            % update covariance and mean of parameters
            sN_r = inv(t*G + l0_r);
            mN_r = sN_r*(t*Z*v + l0_r*m0_r);


            %% updating z

            % estimate some factors
            Wb = W.*(mN_r*mN_r');
            Ws = W.*sN_r;
            A  = Wb + Ws;

            % estimate g
            g = log(p0./(1-p0)) + t*mN_r.*v + t*diag(A)/2;
            
            % define a random order
            ord = randperm(D);
            
            % iterate through all variables
            for i = ord
                z_r(i) = 1;
                g(i) = g(i) - t*z_r'*A(:,i);
                z_r(i) = 1/(1+exp(-g(i)));
            end


            %% updating mN, sN

            % build random matrix Z
            Z = diag(z_r);

            % re-estimate expectation (over Z) of ZX'XZ
            G = Z*W*Z; 
            G(eye(D)>0) = z_r.*diag(W);

            %% updating bN

            % update bN per region
            QF = Y_r'*Y_r/2 - mN_r'*Z*v + mN_r'*G*mN_r/2 + trace(G*sN_r)/2;
            bN_r = b0 + QF;

            % update tau
            t = aN_r/bN_r;
            
            
            % check for sparsity (because of small values)
            mN_r(abs(mN_r)<10^(-5)) = 0;
            z_r(mN_r==0)            = 0;

            % get the "present" connections
            z_idx  = (z_r>pr^2).*(z_r<1)>0;


            %% compute model evidence

            % compute the components of the model evidence
            log_lik = N_eff*(psi(aN_r) - log(bN_r))/2 - N_eff*log(2*pi)/2 - t*QF/2;
            log_p_weight = 1/2*tapas_rdcm_spm_logdet(l0_r) - D*log(2*pi)/2 - (mN_r-m0_r)'*l0_r*(mN_r-m0_r)/2 - trace(l0_r*sN_r)/2;
            log_p_prec = a0*log(b0) - gammaln(a0) + (a0-1)*(psi(aN_r) - log(bN_r)) - b0*t;
            log_p_z = sum(log(1-p0(z_idx)) + z_r(z_idx).*log(p0(z_idx)./(1-p0(z_idx))));
            log_q_weight = 1/2*tapas_rdcm_spm_logdet(sN_r) + D*(1+log(2*pi))/2;
            log_q_prec = aN_r - log(bN_r) + gammaln(aN_r) + (1-aN_r)*psi(aN_r);
            log_q_z = sum(-(1-z_r(z_idx)).*log(1-z_r(z_idx)) - z_r(z_idx).*log(z_r(z_idx)));

            % compute the negative free energy per region
            logF_temp = real(log_lik + log_p_prec + log_p_weight + log_p_z + log_q_prec + log_q_weight + log_q_z);
            
            % check whether convergence is reached
            if ( ((logF_old - logF_temp)^2 < pr^2) )
                convergence = 1;
            end
            
            % store old negative free energy
            logF_old = logF_temp;

            % decrese cycle counter and check for end
            count = count - 1;
            
            % end optimization when number of iterations is reached
            if count<0
                break;
            end
        end


        % check for sparsity (because of small values)
        mN_r(abs(mN_r)<10^(-5)) = 0;
        z_r(mN_r==0)            = 0;
        
        % get the "present" connections
        z_idx  = (z_r>pr^2).*(z_r<1)>0;


        %% compute model evidence

        % compute the components of the model evidence
        log_lik_iter(iter)      = N_eff*(psi(aN_r) - log(bN_r))/2 - N_eff*log(2*pi)/2 - t*QF/2;
        log_p_weight_iter(iter) = 1/2*tapas_rdcm_spm_logdet(l0_r) - D*log(2*pi)/2 - (mN_r-m0_r)'*l0_r*(mN_r-m0_r)/2 - trace(l0_r*sN_r)/2;
        log_p_prec_iter(iter)   = a0*log(b0) - gammaln(a0) + (a0-1)*(psi(aN_r) - log(bN_r)) - b0*t;
        log_p_z_iter(iter)      = sum(log(1-p0(z_idx)) + z_r(z_idx).*log(p0(z_idx)./(1-p0(z_idx))));
        log_q_weight_iter(iter) = 1/2*tapas_rdcm_spm_logdet(sN_r) + D*(1+log(2*pi))/2;
        log_q_prec_iter(iter)   = aN_r - log(bN_r) + gammaln(aN_r) + (1-aN_r)*psi(aN_r);
        log_q_z_iter(iter)      = sum(-(1-z_r(z_idx)).*log(1-z_r(z_idx)) - z_r(z_idx).*log(z_r(z_idx)));

        % compute the negative free energy per region
        logF_iter(iter) = real(log_lik_iter(iter) + log_p_prec_iter(iter) + log_p_weight_iter(iter) + log_p_z_iter(iter) + log_q_prec_iter(iter) + log_q_weight_iter(iter) + log_q_z_iter(iter));
        
        % asign the iteration-specific values
        mN_iter{iter} = mN_r;
        z_iter{iter}  = z_r;
        l0_iter{iter} = l0_r;
        sN_iter{iter} = sN_r;
        aN_iter{iter} = aN_r;
        bN_iter{iter} = bN_r;
        
        % clear the variables
        clear mN_r z_r l0_r sN_r aN_r bN_r
        
    end
    
    % get the logF for all iterations
    [~, best] = max(logF_iter);
    
    % store region-specific parameters
    mN(k,:)             = real(mN_iter{best});
    z(k,:)              = real(z_iter{best});
    sN{k}            	= real(sN_iter{best});
    aN(k)               = real(aN_iter{best});
    bN(k)               = real(bN_iter{best});
    logF(k)             = logF_iter(best);
    
    % get the predicted signal from the GLM (in frequency domain)
    yd_pred_dcm_fft(:,k)	= X(:,:,k) * mN_iter{best};
    
    % store the region-specific components of logF (for diagnostics)
    if ( args.diagnostics_logF )
        logF_term.log_lik(k)        = log_lik_iter(best);
        logF_term.log_p_weight(k)   = log_p_weight_iter(best);
        logF_term.log_p_prec(k)     = log_p_prec_iter(best);
        logF_term.log_p_z(k)        = log_p_z_iter(best);
        logF_term.log_q_weight(k)   = log_q_weight_iter(best);
        logF_term.log_q_prec(k)     = log_q_prec_iter(best);
        logF_term.log_q_z(k)        = log_q_z_iter(best);
    else
        logF_term = [];
    end
    
    
    % clear the region-specific parameters
    clear mN_iter z_iter l0_iter sN_iter aN_iter bN_iter logF_iter
    clear log_lik_iter log_p_weight_iter log_p_prec_iter log_p_z_iter log_q_weight_iter log_q_prec_iter log_q_z_iter
    
    % -- compute error --
    if options.type == 'S'
        y_Ob = DCM.y_GLM_mean(:,k);
    elseif options.type == 'E'
        y_Ob = DCM.y_GLM_mean(1:end-1,k);
    end

    % get the error between prediction and GLM estimate
    y_Error_k(k,1) = norm( ifft( yd_pred_dcm_fft(:,k) ) - ...
        y_Ob ) / std (DCM.y_GLM_mean(1:end-1,k)) ;

   end

%%% for debugging %%
% for k = 1 : nR 
%     ytmp2(k,:,:) = y_Error(k,:,:)*std (y_GLM(1:end-1,k));
%     ytmp         = squeeze(y_Error(k,:,:));
%     [lin,col]    = find ( ytmp == min(min(ytmp)) );
% %     [lin,col]    = find ( y_Error(k,:,:) == min(min(y_Error(k,:,:))) );
%     Pro_best (k) = Pro(lin);
%     Freq_best(k) = Freq(col);
%     AC_best(k,:) = AC(k,:,lin,col);
% end
% 
% for p = 1 : length(Pro)
%     for f = 1 : length(Freq)

y_Error = sum (y_Error_k);

%     end
% end
% [lin,col]    = find ( y_Error_k == min(min(y_Error_k)) );
% AC_best      = AC(:,:,lin,col);

% -- for debugging --
% for nR = 1 : nr 
%     Beta(:,nR) = DCM.xY(nR).Beta;
% end
% y_GLM = DCM.SPM.xX.X(:,1:end-1) * Beta (1:end-1,:);
% y_mean = y_GLM - repmat (mean(y_GLM),ny,1);
% region  = 2;
% Y_Right = X(:,:,region)*AC_best(region,:).' ;
% DY_left = fft( -diff(y_mean(1:end,region)) );
% y_Right = ifft(Y_Right);
% 
% figure (2);
% subplot(4,1,1);
% plot (abs(Y_Right),'b');
% hold on;
% plot (abs(DY_left),'r');
% hold on;
% plot (abs(DYBold_left(:,region)),'g');

% subplot(4,1,2);
% plot (real(Y_Right),'b');
% hold on;
% plot (real(DY_left),'r');
% 
% 
% subplot(4,1,3);
% plot (imag(Y_Right),'b');
% hold on;
% plot (imag(DY_left),'r');
% 
% 
% subplot(4,1,4);
% plot (y_mean(1:end-1,region),'r');
% hold on;
% plot (y_Right,'b');


% write the results to the output file
output{1} = DCM_store_parameters(DCM, mN, sN, aN, bN, logF, logF_term, z, args);

% store the priors
output{1}.priors.m0 = m0;
output{1}.priors.l0 = l0;
output{1}.priors.a0 = a0;
output{1}.priors.b0 = b0;
output{1}.priors.p0 = args.Pro;

% store the rDCM variant
output{1}.inversion = 'tapas_rdcm_sparse';

% store the true and predicted temporal derivatives (in frequency domain)
output{1}.Y.yd_source_fft   = Y;
output{1}.Y.yd_pred_dcm_fft = yd_pred_dcm_fft;

output{1}.Error.y_Error_k = y_Error_k;
output{1}.Error.y_Error   = sum( y_Error_k );


end

function YF = HPFilter (Y,fr) 
% This function perform a high-pass filter
    L_y = size (Y,1);
    y_O = Y (fr : floor(L_y/2)+1,:,:);
    YF  = [real(y_O); imag(y_O)] ; 
end

