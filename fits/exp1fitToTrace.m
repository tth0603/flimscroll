function    pc=exp1fitToTrace(inarg,data)
%pc=exp1fitToTrace(inarg,data,yInt,yShift)
% This function uses a maximum likelihood algorithm to fit a lifetime distribution, 
% taking into account all observations, whether or not a dissociation event
% was detected.  In this instance we will fit a single exponential decay.
% User must input two sets of intervals, the dwell-time intervals ('dwellts')
% from all dissociation events and the times of all observations for which 
% a dissociation event was not detected ('obsts'). The minimum resolution 
% time for the distribution tm is also taken as an input.
% 
% The contribution to the likelihood from the dissociation interval data
% is:
% (1/tau)*exp (-dwellts/tau)/exp(-tm/tau)
%
% The contribution to the likelihood from the non-dissociated events is:
% exp (-obsts/tau)/exp(-tm/tau)
%
% inarg = tau, fit parameters in the distribution, defined in
%                the above equation for the distribution  
% dwellts = vector list of intervals measured (e.g. residence times) that
%            should be part of the distribution 
%            (1/tau)*exp(-t/tau)
% obsts = vector list of intervals measured that did not end in a
%            dissociation event and therefore either contribute to 
%            estimating the lower limit on the dissociation distribution  
%            or to estimating the stable fraction
% tm== minimum interval length that can be resolved in the experiment

%
% call via fminsearch('exp1_mxlall',inargzero,[],dwellts,obsts,tm)

                                  % Form the vector indicating the probabilities that
                                  % the intervals are part of the
                                   % trial distribution.

tau=abs(inarg);

probvector=exp(-data/tau); 

prodprob=sum(log(probvector));   
% Take product of all probabilities;
pc=-prodprob;                         
% The fminsearch will then minimize the
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