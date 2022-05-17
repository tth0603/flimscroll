function mib = twoStepBindFitV5(inarg,intervals,emptyInts,Nt)
% mib = twoStepBindFitV4(inarg,intervals,tm,tx,Nt)
%
% fits a two step initial passage where each of the steps are indie. we try
% something to account for the ints in which we do not see a subsequent
% activation for (the ones that occur at the end of an experiment).
%
% Fit to: A*exp(t*tau1) + (1 - A)*exp(t*tau2) [time constants]
%
% THis was a failure
%
% INPUTS
%   inarg = [ap Tau0]; note that Af = ap^2/(1+ap^2)
%   intervals = '-3' events shifted so the first starts at the origin
%
%USAGE
%   fitParams = fminsearch('twoStepBindFitV3',fitInputs,[],first,tm,tx,mx);
%
% Timothy Harden 2020
%

Af = inarg(1)^2/(1+(inarg(1))^2);
A = inarg(2)^2/(1+(inarg(2))^2);
tau1 = abs(inarg(3));
tau2 = abs(inarg(4));

% number of events:
N = size(intervals,1);
%number of nukes without an event:
Ns = Nt - N;

% probability of events we see:
% probability_vector = ( 1/( exp(-tm/k0) - exp(-tx/k0)) )*(A*Nt)*(1/k0)^2*intervals.*exp(-intervals/k0); %this is from v3 as a template
probability_vector = (Af)*(A*(1/tau1)*exp(intervals/tau1) + (1 - A)*(1/tau2)*exp(intervals/tau2)); %time constants

% probability a nuclei was active but never turned on:
% probNs = ((Nt - A*Nt)/Nt) + A*exp(-tx/k0); %from LJF exp1_global_active_fraction_nonspecific_background_mxl.m. See his email on 20190826
% (inactive #)/(total #)+  (# active w/o landing prior to tx)/(total AOI #)
% probNs = ((Nt - Af*Nt)/Nt) + Af*(A*exp(tx/tau1) + (1 - A)*exp(tx/tau2));
probNs = (A*exp(emptyInts/tau1) + (1 - A)*exp(emptyInts/tau2));

% prodprob=sum(log(probability_vector)) + Ns*log(probNs) ;           % Take product of all probabilities; (via LJF)
prodprob=sum(log(probability_vector)) + log(probNs) ;  
mib=-prodprob; %negative to maximize the probability vector