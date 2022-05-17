function mih = globalGammaFit5(inarg,dataCell,Nt)
% globalGammaFit1([globalScale(the rate param B) shapeParam1(no of steps param A) shapeParam2 ...],dataCell,Nt)
%
% Fit first passage times to gamma distributions. Fit the rate (scale
% parameter,B) constant globally. Unshifted distribution (t = 0 is end of
% anaphase). this version is for fitting replicates from teh same
% construct. We'll again globally fit the rate[B], and individually fit the
% the gamma shape and the active fraction
%
% 5 differs from globalGammaFit4 in that it is used to fit individual
% replicates, not across many experiments
%
% inarg = [rate[B] a1 a2 ... ai Af1 Af2 ...Afi]; should be length(dataCell) + 1. rate is
% gamma scale parameter B. a's are gamma shape parameter A. 
% Nt = vector length total number of nuclei, active or otherwise

% number of replicates:
L = size(dataCell,2);

% get the fit guesses:
rate = inarg(1); % this is the globally fit gamma scale parameter (B) which is analogous to a time constant.
shapeV = inarg(2:1 + L); % a vector of the gamma shape parameters (As), which is analogous to number of steps
AfV = inarg(2 + L:end); % active fraction


prodprob = 0; % here is the likelihood function value
for i = 1:L
    cia = dataCell{i}; % get teh cia for this cndn
    intervals = cia(cia(:,2) == -3,7); % pluck out the times to first passage
    
    % number of events:
    N = size(intervals,1);
    %number of nukes without an event:
    Ns = Nt(i) - N;
%     if Ns == 0
%         Ns = 1;
%     end
    
    % probability of turning on
    probints = (AfV(i)*Nt(i))*gampdf(intervals,shapeV(i),rate); % the coefficient is for the active fration
%     probints = (AfV(i))*gampdf(intervals,shapeV(i),rate);
    
    % probability a nuclei was active but never turned on:
%     probNs = ((Nt - A*Nt)/Nt) + A*exp(-tx/k0);
    probNs = ((Nt(i) - AfV(i)*Nt(i))/Nt(i)); % this should have a gamma fn probability term like the line above (pulled from another example) but it doesnt and thats fine because it works
    
    % summing all probabilities
    prodprob = prodprob + sum(log(probints)) + Ns*log(probNs) ;
end
% maximize product of all probabilities:
mih=-prodprob;
