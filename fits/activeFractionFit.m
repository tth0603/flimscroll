function mih = activeFractionFit(varargin)
% out = activeFractionFit(iniRateBSvector,cia,N,nukesMat)
%
% This goes along with firstEventBootstrapFitter to determine the error
% on the active fraction of nukes under a two-equal-steps model. See
% t1p54e.m for usage. 
%
% INPUTS
%
% OUTPUTS
%
% USAGE
%
% Timothy Harden 2020

%initiate
rateV = varargin{1};
cia = varargin{2}; %to get time of first event
N = varargin{3};
nukesMat = varargin{4}; % to get the length of the experiment

%get number of BSs
numBS = size(rateV,1);

%get teh length:
firstTimes = cia(cia(:,2) == -3,7); %a vector
L = ceil(max(nukesMat(:,1))) - min(firstTimes) + 1; %a scaler
% L = 3000; %make time constant
% for each BS rate, simulate an event lenght:
dim = [N,1];
px = 0:3*L; %this needs to be large, but if you make it too large the ourcome doesn't change but the compute time gets long
fraction = zeros(numBS,1);
for i = 1:numBS
    tau0 = rateV(i);
    p = (1/tau0)^2*px.*exp(-px/tau0);
    simEventLength = randpdf(p,px,dim); % a scritp I got from the web that works well
    % see how many of our simulated 
    successfulEvents = length(simEventLength(simEventLength < L));
    % fraction of nukes that turn on:
    fraction(i) = successfulEvents/N;
end
mih.fractionV = fraction;
mih.mean = mean(fraction);
mih.std = std(fraction);
mih.L = L;

% timeV = 
% (k0)^2*t*exp(-t*k0)

