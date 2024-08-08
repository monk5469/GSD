function [output, options] = DCM_Estimate(DCM, options)
% [output, options] = DCM_Estimate(DCM, options)
% 
% Main analysis function, which calls the subfunctions necessary to run a
% regression DCM analysis.
% 
%   Input:
%   	DCM             -  model structure 
%       options         - estimation options (if empty would be set to default)
%
%   Output:
%       output          - output structure
%       options         - estimation options
% ----------------------------------------------------------------------
% 
% Authors: Haifeng Wu (whf5469@gmail.com)
% 
% Copyright (C) 2023 School of Electricial & Information Technology
%                         Yunnan Minzu University



% display
fprintf('\n========================================================\n')
fprintf('Regression dynamic causal modeling \n')
fprintf('========================================================\n\n')


% create the regressors
[X, Y, DCM] = DCM_Create_Regressors(DCM);


% display start
fprintf('Run model inversion\n')



% whether searche the best value or not
args.search.flag = options.search.flag;
if args.search.flag == 0
    args.Pro    = options.p0;
    args.Freq   = options.freq_cutoff;

    % get the regression result
    output = DCM_Sparse(DCM, X, Y, options, args);

else
    reverseStr = '';
    for f = 1:length(options.search.freq_all)    
       
        args.Freq = options.search.freq_all(f) ;

        for p = 1:length(options.search.p0_all)
            
            args.Pro = options.search.p0_all(p);
            msg = sprintf('Processing: Freq %d/%d, Pro %d/%d', ...
                   f, length(options.search.freq_all), ...
                   p, length(options.search.p0_all));
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

            % get the regression result
            output_all (p,f) = DCM_Sparse(DCM, X, Y, options, args);
            
            % get a cost function, i.e. an error function
            Error_all  (p,f) = output_all{p,f}.Error.y_Error;

        end
    end

    % find optimal hyperparameter 
    [lin,col]    = find ( Error_all == min(min(Error_all)) );
   
    % asign results
    output = output_all{lin,col};

end        

% display finalizing
fprintf('\nFinalize results\n')

% store version
output{1}.ver     = '2023_v1106';

end

