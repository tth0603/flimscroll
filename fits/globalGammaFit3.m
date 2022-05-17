function mih = globalGammaFit3(inarg,dataCell)
% globalGammaFit1([globalMu(locations param) globalScale(the rate, B) shapeParam1(no of steps, A) shapeParam2 ...],dataCell)
%
% Fit first passage times to gamma distributions. Fit the rate (scale
% parameter,B) constant globally. This time the distributions are shifted,
% represented by the location parameter mu, which is fit globally
%
% fits to: p = ( ((t - mu)/B)^(A - 1)*exp( -(t - mu)/B ) )/B*gamma(A);
%
% A (gamma) is shape (no of steps), B (beta) is scale (rate), mu is location
%
% inarg = [mu rate[B] a1 a2 ...]; should be length(dataCell) + 2. rate is
% gamma scale parameter B. a's are gamma shape parameter A. 
% if inarg(2) < 0
%     keyboard
% end
inarg = abs(inarg); % there are better ways to ensure the inputs are positive, and still the outputs can be negative but whatever
mu = inarg(end); % this is the globally fit position parameter. represents a shift along the x axis
% mu = 0;
rate = inarg(1); % this is the globally fit gamma scale parameter (B) which is analogous to a time constant.
shapeV = inarg(2:end - 1); % a vector of the gamma shape parameters (As), which is analogous to number of steps
L = size(dataCell,2);
prodprob = 0; % here is the likelihood function value
for i = 1:L
    cia = dataCell{i}; % get teh cia for this cndn
    intervals = cia(cia(:,2) == -3,7); % pluck out the times to first passage
    probints = ( ((intervals - mu)./rate).^(shapeV(i) - 1).*exp( -(intervals - mu)./rate ) )./(rate*gamma(shapeV(i)));
    prodprob = prodprob + sum(log(probints));
end
% maximize product of all probabilities;
mih=-prodprob;
