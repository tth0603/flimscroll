function mih = initialFraction6(dataVector,normalization,fitInputs,tm,tx)
% mih = initialFraction6(dataVector,normalization,fitInputs,tm,tx)
% 
% Same as initialFraction5, but only makes the fit. Useful. Used for
% fitting bootstraps--firstEventBootstrapFitter.m
%
% INPUTS
%   dataVector == e.g. from firstEventBootstrapGenerator.m
%   normalization == a scaler. the number we are normalizing the initial
%   bind curve to. if '1' we normalize to the number of events in cia, if
%   empty '[]' we do not normalize. 
%   fitInputs == model initial guesses. [Af k0] 
%   tm == minimum time resolution
%   tx == maximum time an event can take place
%
% OUTPUS
%   MIH = a two member struc with members: Af (scaler), k0 (scaler)
%
% Timothy Harden 2020

%determine normalization
if normalization == 1
    mx = size(first,1);
elseif ~isempty(normalization)
    mx = normalization;
else 
%     mx = length(first);
    mx = 1;
end

%fit:
fitParams = fminsearch('twoStepBindFitV3',fitInputs,[],dataVector,tm,tx,mx);
ap = fitParams(1);
Af = ap^2/(1 + ap^2);
k0 = fitParams(2);

mih.Af = Af;
mih.k0 = k0;

    
     
