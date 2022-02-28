function mih = firstEventBootstrapFitter(varargin)
%   mih = firstEventBootstrapFitter(bsMat,normalization,fitInputs)
%
% fits a matrix of bootstraps to a two step initial bind model. 
%
% INPUTS
%   bsMat == each col is a BS
%   normalization == scaler. what to normalize to. can be empty, and the curve is not normalized
%   fitInputs == initial guesses
%
% OUTPUTS
%   mib == an numberOfBS x 2 matrix where the first columns is Af and the
%   second column is tau0
%
% Timothy Harden 2020

%initiate:
bsMat = varargin{1};
normalization = varargin{2};
fitInputs = varargin{3};
numBS = size(bsMat,2);

fitMat = zeros(numBS,1);
for i = 1:numBS
    fitParams = initialFraction6(bsMat(:,i),normalization,fitInputs,30,3000);
    fitMat(i) = fitParams.k0;
end
mih = fitMat;
