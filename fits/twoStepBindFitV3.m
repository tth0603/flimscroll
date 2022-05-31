function mib = twoStepBindFitV3(inarg,intervals,tm,tx,Nt)
% mib = twoStepBindFitV3(inarg,intervals,tm,tx,Nt)
%
% fits a two step initial passage where each of the steps has the same
% rate. Here we account for limited observation window and properly fit the
% fraction of active nuclei by determining the probability a nuclei didn't
% run on given the fitted rate constant.
%
% Fit to: A*1/tau0*t*exp(-time/tau0)
%
% INPUTS
%   inarg = [ap Tau0]; note that Af = ap^2/(1+ap^2)
%   intervals = '-3' events shifted so the first starts at the origin
%   Nt = total number of nukes
%
%USAGE
%   fitParams = fminsearch('twoStepBindFitV3',fitInputs,[],first,tm,tx,mx);
%
% Timothy Harden 2020
%

A = inarg(1)^2/(1+(inarg(1))^2);
k0 = abs(inarg(2));

% number of events:
N = size(intervals,1);
%number of nukes without an event:
Ns = Nt - N;

% probability of events we see:
probability_vector = ( 1/( exp(-tm/k0) - exp(-tx/k0)) )*(A*Nt)*(1/k0)^2*intervals.*exp(-intervals/k0); %time constants
% the term out frot corrects for the limited time window of observation

% probability a nuclei was active but never turned on:
probNs = ((Nt - A*Nt)/Nt) + A*exp(-tx/k0); %from LJF exp1_global_active_fraction_nonspecific_background_mxl.m. See his email on 20190826
% (inactive #)/(total #)+  (# active w/o landing prior to tx)/(total AOI #)
% (Af*Nt - Nz)/Nt )

prodprob=sum(log(probability_vector)) + Ns*log(probNs) ;           % Take product of all probabilities; (via LJF)
mib=-prodprob; %negative to maximize the probability vector