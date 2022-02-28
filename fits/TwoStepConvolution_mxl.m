function pc=TwoStepConvolution_mxl(inarg,intervals,tm,tx)
%
% function TwoStepConvolution_mxl(inarg,intervals,tm,tx)
%
% Will use a MAXIMUM LIKELIHOOD algorithm for fitting a distribution
% function.  In this instance we will fit a model consisting of a two step 
% process, i.e. a convolution of two exponentials with rates k1 (step one)
% and k2 (for step two).  We must also input the minimum
% resolution time for the distribution, tm and the maximum time interval tx.
%  
%(  k1*k2/( A - B )*...
%          ( exp(-intervals*k2)-exp(-intervals*k1) );
%
%  where A = k1*( exp(-tm*k2) - exp(-tx*k2) ) 
%    and B = k2*( exp(-tm*k1) - exp(-tx*k1) )            see B30p43
% 
%        
% inarg = [ k1 k2],  starting fit parameters in the distribution, 
%                defined in the above equation for the distribution
% intervals == vector list of intervals measured (e.g. residence times) that
%            should be part of the distribution 
% tm == minimum interval length that can be resolved in the experiment
% tx== maximum interval length that can be resolved in the experiment
%
% call via fminsearch('TwoStepConvolution_mxl',inargzero,[],intervals,tm,tx)

%N=length(intervals);                    % Number of events measured in the experiment
                                % Calculate the probability of each
                                % measured interval being part of the
                                % hypothesized distribution

k1=abs(inarg(1));
k2=abs(inarg(2));
A = k1*( exp(-tm*k2) - exp(-tx*k2) );
B = k2*( exp(-tm*k1) - exp(-tx*k1) );

                                       % Expression in terms of time constants
%A=( exp(-tm/tau1) - exp(-tx/tau1) );
%B = ( exp(-tm/tau2) - exp(-tx/tau2) );
%probability_vector=(  1/( a*A + (1-a)*B  )*...
%                                 ( a/tau1*exp(-intervals/tau1)+(1-a)/tau2*exp(-intervals/tau2) )+...
%                                   0);
                                        % Expression in terms of rates

probability_vector=(k1*k2/( A - B ))*( exp(-intervals*k2)-exp(-intervals*k1) );
                                
                  %               86.4617*exp(-intervals/1.3642) + 14.1050*exp(-intervals/5.9616);
               % Corrects for (reanalysis) nonspecific binding b18p54 top of page
              % 82.94*exp(-intervals/1.45) + 12.22*exp(-intervals/6.13);  
              % Corrects for nonspecific binding b18p42 top of page
           

prodprob=sum(log(probability_vector));           % Take product of all probabilities;
pc=-prodprob;                         % The fminsearch will then minimize the
                                      % -prodprob (i.e. maximize prodprob)
                                      % so we will maximize the likelihood
                                      % that the intervals vector
                                      % represents the said distribution
                                      % (The sum of the log of the vector entries
                                      % is the same as the log of the
                                      % product of the vector entries.  We
                                      % use the sum( log() ) because taking
                                      % the product of the entries yields a
                                      % number too small for the computer
                                      % to handle.
