function [X, Y, DCM] = DCM_Create_Regressors(DCM)
% [X, Y, DCM] = DCM_Create_Regressors(DCM)
% 
% Transforms intial DCM signal into a set of regressors X (design matrix)
% and Y (data)
% 
%   Input:
%   	DCM         - model structure
%
%   Output:
%   	X           - design matrix (predictors)
%       Y           - data
%   	DCM         - either model structure or a file name
%
 
% ----------------------------------------------------------------------
% 
% Authors: Haifeng Wu (whf5469@gmail.com)
% 
% Copyright (C) 2023 School of Electricical & Information Technology
%                         Yunnan Minzu University
% ----------------------------------------------------------------------


%% preparations


NR   = DCM.n;         % the number of regions
u    = full(DCM.U.u); % experimental stimulation
nu   = size (u, 2);   % the number of stimulations
r_dt = DCM.Y.dt/DCM.U.dt; % the coefficient for downsampling

% whether GLM estimation are perfomred or not 
if  isfield(DCM,'SPM') && isfield(DCM,'xY') 
    
    % GLM has performed by SPM software and Beta needs to be extracted
    h = DCM.SPM.xBF.bf; % Hemodynamic response function
    X = DCM.SPM.xX.X;
    for nR = 1 : NR 
        Beta(:,nR) = DCM.xY(nR).Beta; 
    end
    y_GLM = X(:,1:end-1) * Beta (1:end-1,:); % GLM output
    x     = u * Beta(1:end-1,:); % neuronal activity
    X_O  = [y_GLM X] ;           % designed matrix on time-domain
else
    
    % GLM need to perform from a function, GLM_filter( )
    DCM.Y.y_GLM = GLM_filter (DCM.Y.y,DCM); % GLM output
    y_GLM       = DCM.Y.y_GLM;
    h        = DCM.h;                      % Hemodynamic response function
    h_fft    = fft( h );
    u        = ifft(fft(u).* repmat(h_fft, 1, nu));
    DCM.U.X0 = ones(size(DCM.U.u,1),1);    % for confounds
    u        = [u, DCM.U.X0];
    u        = u(1 : r_dt : end, : );
    X_O      = [y_GLM u/(r_dt)];           % designed matrix on time-domain
end

Ny    = size(y_GLM,1);
Nu    = size(DCM.b,3);
X_O   = repmat (X_O, [1,1,NR]);

DCM.cb = DCM.c; % generate a new matrix of cb, involving b 

for nR = 1 : NR
    
    if sum(DCM.b(nR,:,:),'all') ~= 0 % judge what region exists bi-linearity
       
        if  NR - sum( DCM.cb(nR,:) ) < sum( DCM.b(nR,:,:), "all" )
            error (['Coefficient C can not be used for ' ...
                    'the estimation of B because of its space!']);
        else         
            
            for nu = 1: Nu
                
                if sum( DCM.b(nR,:,nu) ) ~= 0 % judge what input exist bi-linearity
                    % -- generate the item of bi-linearity UXH --
                    bCol = find ( DCM.b(nR,:,nu) );
                    ux   = repmat( u(:,nu), 1, length(bCol) ).* x(:,bCol);
                    % ux   = repmat( u(:,nu), 1, Nu ).* x(:,bCol);
                    uxh  = convmv ( ux , h ) ;
                    UXH  = uxh(1:r_dt: Ny*r_dt,:,:);
                    % --------------------------------------
                    
                   cCol = find ( ~ DCM.cb(nR,:)  ); % Find where the part of
                                                     % coefficient c being 0 
                                                     % will be replaced with UXH
                    for c = 1:length(bCol)
                        X_O(:,NR+cCol(c),nR) = UXH(:,c); % The designed matrix's colums where
                                                         % the c matrix's colums are 0 
                                                         % will be replaced with UXH
                        
                        DCM.cb(nR,cCol(c))   = 1; % update the new matrix cb
                    end
                     
                end
            
            end
             
        
        end
              
    end
end


   
% whether large TR or small TR
if  DCM.dt_flag == 'L'
    
    %  fft( dy ) = fft( y(t) - y(t+1) )
    Y              = ( fft(-diff( y_GLM ) ) );
    X_O            = X_O - repmat( mean(X_O), Ny,1);
    X              = fft( X_O(1:end-1,:,:) ) ;
    
elseif DCM.dt_flag == 'S'

    %  if dy <--> DY, then jw * Y = DY 
    ic    = sqrt(-1);
    coef  = exp(2 * pi * ic * (0 : Ny - 1)' / Ny) - 1;
    y_fft = fft( y_GLM );
    Y     = repmat(coef, 1, NR).*  y_fft / DCM.Y.dt;
    u_fft = fft(u);
    X     = [y_fft u_fft/(r_dt)];
    X     = repmat (X, [1,1,NR]);

end

DCM.r_dt       = r_dt;
DCM.y_GLM      = y_GLM;
DCM.y_GLM_mean = y_GLM - repmat(mean(y_GLM),Ny,1);

% %---- for debugging ---------
% y_mean= y_GLM - repmat(mean(y_GLM),Ny,1);
% region = 1;
% % AC     = [0.7706	0.4874	0        1.4087 0      0
% %           -0.06	    0.5581	-0.5676  0.5879 0.4863 0
% %                0	0.3218	0.2098   0      0      0];
% AC       = [DCM.a DCM.c];
% Y_Right = X(:,1:end-1,region)*AC(region,:).' ;
% Y_Right = Y;
% y_Right = ifft(Y_Right);
% 
% Y       = ( fft(-diff( DCM.Y.y) ) );
% DY_left = Y;
% 
% 
% 
% 
% figure (1);
% subplot(4,1,1);
% plot (abs(Y_Right(:,region)),'b');
% hold on;
% plot (abs(DY_left(:,region)),'r');
% % hold on;
% % plot (abs(DYBold_left(:,region)),'g');
% 
% 
% subplot(4,1,2);
% plot (real(Y_Right),'b');
% hold on;
% plot (real(DY_left(:,region)),'r');
% % hold on;
% % plot (real(DYBold_left(:,nr)),'g');
% 
% subplot(4,1,3);
% plot (imag(Y_Right),'b');
% hold on;
% plot (imag(DY_left(:,region)),'r');
% % hold on;
% % plot (imag(DYBold_left(:,nr)),'g');
% 
% subplot(4,1,4);
% plot (y_mean(1:end-1,region),'b');
% hold on;
% plot (y_Right,'r');


% -- high pass filter for debuging---
% L_y                                = size (Y,1);
% idx_freq                           = [70; 70; 70;]';
% idx                                = logical( zeros (L_y,1) );
% for nR = 1 : NR 
%     idx (idx_freq(nR): end-idx_freq(nR)+2,nR) = 1;
%     Y(~idx(:,nR),nR)                          = NaN;  
% end
% --


end


function y = convmv(A,h)
% This function is for a convolution between a vecotr and a matrix

for col = 1 : size (A,2)
    y(:,col) = conv(A(:,col),h);
end

end
