function mib = twoStepFirstBind(inarg,intervals)
%
% in beta. see also
% exp1_global_active_fraction_nonspecific_background_mxl.m
%
% the form of the curve to be fit is:
%   af*(1 - exp(-k1*t)).*(1 - exp(-k2*t));
% but we are actually fitting:
%   k1k2/(k1 - k2)*( exp(-k2*t) - exp(-k1*t) );  

% af = inarg(1);
% k1 = inarg(2);
% k2 = inarg(3);

k1 = inarg(1);
% k2 = inarg(2);

% Nt = 130;

% probvector = af/Nt*k1*k2/(k1 - k2).*( exp(-k2.*intervals) - exp(-k1.*intervals));  
% probvector = ( exp(-k2.*intervals) - exp(-k1.*intervals));  
% if k1 ~= k2
%     probvector = k1*k2/(k1 - k2).*( exp(-k2.*intervals) - exp(-k1.*intervals)); 
% else
    probvector = k1^2.*intervals.*exp(-k1.*intervals);
% end

prodprob=sum(log(probvector));           % Take product of all probabilities;
mib = prodprob;  
