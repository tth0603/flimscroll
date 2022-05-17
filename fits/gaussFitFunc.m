function pc=gaussFitFunc(inarg,data)
%
% function gaussFitFunc(inarg,data)
%
% This function will be called from using the
% lscurvefit() function to a gaussian.  The form of the
% function is:
% 1/sqrt(2*pi*sigma^2)*exp( -( (x-mean).^2/(2*sigma^2))
%
% where the input arguements are given in inarg according to
% inarg  == [ amp mean sigma  ]
% data == x values which are the centers of the bins for the histogram
%
%Call using:
%   out =lsqcurvefit('gaussFitFunc',[10 0.0001 0.1],xValues,yValues);
%       where yValues are the probabilities from teh histogram
%   out == a vecter of amplitude, mean and width


% pc= inarg(1)*exp( -( (data-inarg(2)).^2/(2*inarg(3)^2) ) );     % inputs: [amp mu sigma]... dumb
pc = 1/sqrt(2*pi*inarg(2)^2)*exp( -( (data-inarg(1)).^2/(2*inarg(2)^2) ) );     % inputs: [mu sigma]

%to fit data with a gaussian and plot:
%out =lsqcurvefit('gaussFitFunc',[10 0.0001 0.1],pointXvalues,normFactor);
% x = [-0.3:0.001:0.3];
% fitCurve = out(1) * exp(-( (x-out(2)).^2/(2*out(3)^2) ) );
% hold on; plot(x,fitCurve,'r-');shg
