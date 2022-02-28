function pc=TwoStepConvolution_mxl_eval(k1,k2,intervals,tm,tx)
%
% function TwoStepConvolution_mxl_eval(inarg,intervals,tm,tx)
%
% This is a tool for use with TwoStepconvolution_mxl function.  The current
% function will generate a vector of values for the two step distribution
% function using the functional form:
%  
%(  k1*k2/( A - B )*...
%          ( exp(-intervals*k2)-exp(-intervals*k1) );
%
%  where A = k1*( exp(-tm*k2) - exp(-tx*k2) ) 
%    and B = k2*( exp(-tm*k1) - exp(-tx*k1) )            see B30p43
% 
% Use this to generate a table of probabilities that can then be used with
% probability_steps() to generate a set of 'fake data' for testing out the
% TwoStepConvolution_mxl() function.
%        
% inarg = [ k1 k2],  starting fit parameters in the distribution, 
%                defined in the above equation for the distribution
% intervals == vector list of intervals at which the function will be
%              evaluated.  These values should run between tm and tx
% tm == minimum interval length that can be resolved in the experiment
% tx== maximum interval length that can be resolved in the experiment
%

%N=length(intervals);                    % Number of events measured in the experiment
                                % Calculate the probability of each
                                % measured interval being part of the
                                % hypothesized distribution


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
                                
pc=probability_vector;
