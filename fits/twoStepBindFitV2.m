function mib = twoStepBindFitV2(inarg,intervals,tm,tx)
% 
% garbage. see v3 
%
% See initialFraction4.m. fits a two step initial passage where each of the
% steps has the same rate. 
%
% INPUTS
%   inarg = [Af k0]; 
%   intervals = '-3' events shifted so the first starts at the origin
%
%USAGE
%   fminsearch('twoStepBindFitV1',[0.87 0.0020 0.0013],[],intervals);
%
% Timothy Harden 2020
%

A = inarg(1);
k0 = (inarg(2));

probability_vector = A*k0^2*intervals.*exp(-k0*intervals); %rates

prodprob=sum(log(probability_vector));           % Take product of all probabilities; (via LJF)
mib=prodprob; %this is wrong, but this with initialFraction4 works so I'm keeping this for now. See initailFraction5/twoStepBindFitV3 for a better way