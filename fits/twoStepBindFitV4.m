function mib = twoStepBindFitV4(inarg,intervals,tm,tx,Nt)
% mib = twoStepBindFitV4(inarg,intervals,tm,tx,Nt)
%
% fits a two step initial passage where each of the steps are indie. we try
% something to account for the ints in which we do not see a subsequent
% activation for (the ones that occur at the end of an experiment).
%
% Fit to: (k1*k2)/(k2 - k1)*(exp(-t*k1) - exp(-t*k2)) (rates, not time constants)

% INPUTS
%   inarg = [ap Tau0]; note that Af = ap^2/(1+ap^2)
%   intervals = '-3' events shifted so the first starts at the origin
%
%USAGE
%   fitParams = fminsearch('twoStepBindFitV3',fitInputs,[],first,tm,tx,mx);
%
% Timothy Harden 2020
%

A = inarg(1)^2/(1+(inarg(1))^2);
k1 = abs(inarg(2));
k2 = abs(inarg(3));

% number of events:
N = size(intervals,1);
%number of nukes without an event:
Ns = Nt - N;

% probability of events we see:
% probability_vector = ( 1/( exp(-tm/k0) - exp(-tx/k0)) )*(A*Nt)*(1/k0)^2*intervals.*exp(-intervals/k0); %this is from v3 as a template
probability_vector = (A*Nt)*(k1*k2)/(k2 - k1)*(exp(-intervals*k1) - exp(-intervals*k2)); %rate constants

% probability a nuclei was active but never turned on:
% probNs = ((Nt - A*Nt)/Nt) + A*exp(-tx/k0); %from LJF exp1_global_active_fraction_nonspecific_background_mxl.m. See his email on 20190826
% (inactive #)/(total #)+  (# active w/o landing prior to tx)/(total AOI #)
probNs = ((Nt - A*Nt)/Nt) + A*(k1*k2)/(k2 - k1)*(exp(-tx*k1) - exp(-tx*k2));

prodprob=sum(log(probability_vector)) + Ns*log(probNs) ;           % Take product of all probabilities; (via LJF)
mib=-prodprob; %negative to maximize the probability vector