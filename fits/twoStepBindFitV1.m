function mib = twoStepBindFitV1(inarg,intervals)
% See t1p65 & t1p65a.m
%
% INPUTS
%   inarg = [A k0 k1];
%   intervals = '-3' events shifted so the first starts at the origin
%
%USAGE
%   fminsearch('twoStepBindFitV1',[0.87 0.0020 0.0013],[],intervals);
%
% Timothy Harden 2020
%

A = inarg(1);
k0 = inarg(2);
k1 = inarg(3);

probability_vector = A*( k0*k1/(k0 - k1).*( (1/k0)*exp(-k0.*intervals) - (1/k1)*exp(-k1.*intervals) ) + 1 );
% probability_vector = A*k0*k1/(k0 - k1 + 0.000001).*( (1/k0)*exp(-k0.*intervals) - (1/k1)*exp(-k1.*intervals) + 1 ); %to prevent inf 

%only two params
% k0 = inarg(1);
% k1 = inarg(2);
% probability_vector = 0.87*k0*k1/(k0 - k1).*( (1/k0)*exp(-k0.*intervals) - (1/k1)*exp(-k1.*intervals) ) + 1;


prodprob=sum(log(probability_vector));           % Take product of all probabilities; (via LJF)
mib=-prodprob; 