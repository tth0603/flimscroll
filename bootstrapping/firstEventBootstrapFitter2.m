function mih = firstEventBootstrapFitter2(varargin)
%   mih = firstEventBootstrapFitter2(bsCell,fitInputs)
%
% fits a matrix of bootstraps to a gamma initial bind model. 
%
% INPUTS
%   bsCell == each col is a BS. from firstEventBootstrapGenerator2.m
%   fitInputs == initial guesses
%
% OUTPUTS
%   mib == an numberOfBS x 2 matrix where the first columns is Af and the
%   second column is tau0
%
% Timothy Harden 2020

% initialize
bsCell = varargin{1};
inarg = varargin{2};
% the BS number
bsNum = size(bsCell{1},2);
% the number of conditions
cndns = size(bsCell,2);

% Set up some things
NtV = zeros(cndns,1);
for i = 1:cndns
    %get the number of total nukes for each condition
    NtV(i) = size(bsCell{i},1);
end

fitMat = zeros(bsNum,11);
for i = 1:bsNum
    repCell = {};
    for j = 1:cndns
        v = bsCell{j}(:,i); %snag the intervals for this bootstrap
        repCell{j} = v(v > 0); % get rid of the zeros, place in cell
    end
    fitMat(i,:) = fminsearch('globalGammaFit4v2',inarg,[],repCell,NtV);
end

mih.fitMat = fitMat;
mih.mean = mean(fitMat);
mih.std = std(fitMat);











