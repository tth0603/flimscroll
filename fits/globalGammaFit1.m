function mih = globalGammaFit1(inarg,dataCell)
% globalGammaFit1([globalScale(the rate param B) shapeParam1(no of steps param A) shapeParam2 ...],dataCell)
%
% Fit first passage times to gamma distributions. Fit the rate (scale
% parameter,B) constant globally. Unshifted distribution (t = 0 is end of
% anaphase).
%
% inarg = [rate[B] a1 a2 ... ai]; should be length(dataCell) + 1. rate is
% gamma scale parameter B. a's are gamma shape parameter A. 

rate = inarg(1); % this is the globally fit gamma scale parameter (B) which is analogous to a time constant.
shapeV = inarg(2:end); % a vector of the gamma shape parameters (As), which is analogous to number of steps
L = size(dataCell,2);
prodprob = 0; % here is the likelihood function value
for i = 1:L
    cia = dataCell{i}; % get teh cia for this cndn
    intervals = cia(cia(:,2) == -3,7); % pluck out the times to first passage
    probints = gampdf(intervals,shapeV(i),rate);
    prodprob = prodprob + sum(log(probints));
end
% maximize product of all probabilities:
mih=-prodprob;
