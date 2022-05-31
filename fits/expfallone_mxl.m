function pc=expfallone_mxl(tau,intervals,tm,tx)
%
% function expfallone_mxl(tau,intervals,tm,tx)
%
% Will test out a maximum likelihood algorithm for fitting a distribution
% function.  in this instance we will try just a simple exponential with a
% single fit parameter tau.  User must also input the minimum resolution
% time for the distribution, tm and the maximum time interval tx.
%
%       ( 1/( exp(-tm/tau) - exp(-tx/tau) ) )* (1/tau)*exp(-intervals/tau)
%
% tau == exponential time constant   e.g.exp(-intervals/tau)
% intervals = vector list of intervals measured (e.g. residence times) that
%            should be part of the distribution 
%            (1/tau)*exp(-intervals/tau)
% N == length(intervals),  number of binding events in our list
% tm== minimum interval length that can be resolved in the experiment
% tx== maximum interval length that can be resolved in the experiment
%
% call via fminsearch('expfallone_mxl',tauzero,[],intervals,tm,tx)
%N=length(intervals);
                                  % Form the vector indicating the probabilities that
                                  % the intervals are part of the
                                   % trial distribution.
%probvector=( 1/( exp(-tm/tau) - exp(-tx/tau)) )*(1/tau)*exp(-intervals/tau); 
tau=abs(tau);
probvector=( 1/( exp(-tm*tau) - exp(-tx*tau)) )*(tau)*exp(-intervals*tau); 
%probvector=( N/exp(-tm/tau) )*(1/(1-exp(-28/tau) ))*(1/tau)*exp(-intervals/tau);     

prodprob=sum(log(probvector));           % Take product of all probabilities;
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
