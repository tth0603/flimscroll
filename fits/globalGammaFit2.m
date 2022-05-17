function mih = globalGammaFit2(inarg,dataCell)
% globalGammaFit1([globalShape(no of steps) scaleParam1(rate) scaleParam2 ...],dataCell)
%
% Fit first passage times to gamma distributions. Fit the no of steps
% (shape parameter, A) constant globally. Unshifted distribution (t = 0 is end of
% anaphase).
%
% inarg = [shape[A] rate1 rate2 ...]; should be length(dataCell) + 1. rate is
% gamma scale parameter B. a's are gamma shape parameter A. 

shape = inarg(1); % this is the globally fit gamma scale parameter (B) which is analogous to a time constant.
rateV = inarg(2:end); % a vector of the gamma shape parameters (As), which is analogous to number of steps
L = size(dataCell,2);
prodprob = 0; % here is the likelihood function value
for i = 1:L
    cia = dataCell{i}; % get teh cia for this cndn
    intervals = cia(cia(:,2) == -3,7); % pluck out the times to first passage
    probints = gampdf(intervals,shape,rateV(i));
    prodprob = prodprob + sum(log(probints));
end
% maximize product of all probabilities:
mih=-prodprob;
