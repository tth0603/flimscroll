function mib = twoStepFirstPassage(params,vectorOfOnTimes)

%
k1 = params(1);
k2 = params(2);

probabilityVector = 1-(k1*k2/(k1 - k2).*((1/k2).*exp(-k2.*vectorOfOnTimes)-(1/k1).*exp(-k1.*vectorOfOnTimes)));

prodprob = sum(log(probabilityVector));
mib = -prodprob;

end
