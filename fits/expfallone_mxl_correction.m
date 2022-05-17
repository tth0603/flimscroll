function pc=expfallone_mxl_correction(argouts,intervals,tm,tx)
%
% function expfallone_mxl_correction(argouts,intervals,tm,tx)
%
% Companion to expfallone_mxl.   This will calculate the number of events
% that we must add to our cumulative distribution to correct for missed
% events.  That is, if our cumulative distribution is
% cumul=[ time   (# events <= time) ] , then we must add to cumul(:,2) the
% number of events that occur beyond time tx in order compare with the fit.
% This program will output the structure
% pc.CorrectionFactor       % = (  1/A  ) 
% pc.Nzero                  % = actual # of events in data set (observed+not observed)  
% pc.Offset            % # remaining events at time = tx
%
%  so that we need to correct our cumul distribution using
% cumul(:,2)= cumul(:,2)+ pc.Offset;
% 
% Model distribution
%(  1/ A )*( 1/tau1 *exp(-intervals/tau1) );       
%
%  where A = ( exp(-tm/tau1) - exp(-tx/tau1) ) 
% 
%  (see B18p36)
%        
% argouts = [ tau1],  starting fit parameters in the distribution, 
%                defined in the above equation for the distribution

% intervals == vector list of intervals measured (e.g. residence times) that
%            should be part of the distribution 
% tm == minimum interval length that can be resolved in the experiment
% tx== maximum interval length that can be resolved in the experiment
%

N=length(intervals);                    % Number of events measured in the experiment
                                % Calculate the probability of each
                                % measured interval being part of the
                                % hypothesized distribution
tau1=abs(argouts(1));
                                       % Expression in terms of time constants
%A=( exp(-tm/tau1) - exp(-tx/tau1) );
%B = ( exp(-tm/tau2) - exp(-tx/tau2) );
%probability_vector=(  1/( a*A + (1-a)*B  )*...
%                                 ( a/tau1*exp(-intervals/tau1)+(1-a)/tau2*exp(-intervals/tau2) )+...
%                                   0);
                                        % Expression in terms of rates
A=( exp(-tm*tau1) - exp(-tx*tau1) );
%B = ( exp(-tm*tau2) - exp(-tx*tau2) );

pc.CorrectionFactor=(  1/A );
pc.Nzero=pc.CorrectionFactor*N;
pc.Offset=(exp(-tx*tau1) )*pc.Nzero;
