function out = twoExpEqFit(inarg,wtFractions,mutFractions,wtInts,mutInts)
%inarg = [tau Awt0 Amut0];
%
%Usage:
%out = fminsearch('twoExpEqFit',inarg,[],wtFractions,mutFractions,wtInts,mutInts)
%
Q1 = inarg(2)*exp(inarg(1).*wtInts)-wtFractions;
Q2 = inarg(3)*exp(mutInts*inarg(1))-mutFractions;

probProd = sum(log(Q1)) + sum(log(Q2));
%probProd = sum([log(Q1);log(Q2)]);      %this and the above are equivalent

out = probProd;
end