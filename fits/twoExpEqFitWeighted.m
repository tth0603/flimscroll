function out = twoExpEqFitWeighted(inarg,wtFractions,mutFractions,wtInts,mutInts,wtErr,mutErr)
%FAIL
%inarg = [tau Awt0 Amut0];

%Usage:
%out = fminsearch('twoExpEqFitWeighted',inarg,[],wtFractions,mutFractions,wtInts,mutInts,wtErr,mutErr)
wWt = 1./wtErr.^2;
wMut = 1./mutErr.^2;
Q1 = inarg(2)*exp(inarg(1).*wtInts')-wtFractions'.*wWt';
Q2 = inarg(3)*exp(mutInts*inarg(1))-mutFractions*wMut;
%Q1 = Awt*exp(tau.*wtInts')-wtFractions';
%Q2 = Amut*exp(mutInts*tau)-mutFractions;

probProd = sum(log(Q1)) + sum(log(Q2));
%probProd = sum([log(Q1);log(Q2)]);      %this and the above are equivalent

out = probProd;
end