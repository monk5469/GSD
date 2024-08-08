%---------The  author of the original code-----------------
%%Prof Richard (Rik) Henson, MA, MSc, PhD, FBA
%%Professor of Cognitive Neuroscience, Department of Psychiatry, University of Cambridge
%%Director, Cambridge Centre for Ageing and Neuroscience (CamCAN)
%%President-Past, British Neuroscience Association
%%Fellow of the British Academy

%%MRC Cognition and Brain Sciences Unit
%%University of Cambridge
%%15 Chaucer Road
%%Cambridge, CB2 7EF
%%England
%%URL: http://www.mrc-cbu.cam.ac.uk/people/rik.henson/personal
%%Source code URL：https://github.com/SMScottLee/Face_DCM_fMRI
%---------------------------------------------------------------------------

function p=t2p(t,df,tf)

if nargin < 3
    tf = 1; % tail-flag, ie one-tailed by default (assumes sign of t matches direct of hypothesis)
end

try
    p = spm_Tcdf(t,df);
catch
    p = tcdf(t,df);
end

if tf == 1
    warning('One-tailed: Assuming sign of t matches direction of hypothesis!\n')
    p = 1 - p;
else
    f = find(p>0.5);
    p(f) = (1 - p(f));
    p = p * 2;
end

