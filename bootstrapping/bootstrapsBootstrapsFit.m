function mib = bootstrapsBootstrapsFit(varargin)
%   output = bootstrapsBootstraps(dataMatrix,dataColumn,numOfBS,fit,fitParams)
% My latest and greatest bootstrap script. will spit out and plot the 5th &
% 95th percentile bootstrap vector so as to get error bars. Used in
% t1p54b.m
%
% INPUTS
%   dataMatrix == some 2D matrix with x axis values as rows
%   dataColumn == a scaler indicating the column index for your data
%   numOfBS == number of bootstrap data replicates
%   fit == a scaler. '1' for single exp fit. '2' for double
%   fitParams == a vector of: min times resolution, max time resolution,
%   initial guesses for fit:
%       EX: 
%       [1 3600 0.1] For single exp
%       [1 3600 0.9 0.01 0.001] for dbl exp
%
% OUTPUTS
%   a structure with the following members:
%   BSparams ==  an numOfBS x N matrix of fit params where N is the number of fit
%   params for the specified fit. N = 1 for single exp, 3 for dbl exp
%       EX of cols:
%       [k1] for sngle exp
%       [a k1 k2] for dbl exp
%   
%
% USAGE
%   paramMat = bootstrapsBootstraps(dataMat,6,100,2,[1 3600 0.9 0.01 0.001])
%
% Timothy Harden 2020

%initiate
dataMat = varargin{1};
col = varargin{2};
intervals = dataMat(:,col);
numOfBS = varargin{3};
n = size(dataMat,1);
mx = ceil(max(intervals)); %the longest event in data
fit = varargin{4};
fitV = varargin{5};
%bootstrap the data:
tbl=[[1:n]' 1/n*ones(n,1)];  %this is from LJF's method of picking with replacement
%may want to change this to a parfor loop:
% bsMat = zeros(n,numOfBS); %old way
indx = zeros(n,numOfBS);
for i = 1:numOfBS
%     indx = probability_steps(tbl,n); %old way.
%     bs = intervals(indx); 
%     bsMat(:,i) = bs;
    indx(:,i) = probability_steps(tbl,n);
end
bsMat = intervals(indx);
%fit the bootstraps:
if fit == 2 %for dbl exp fit:
    BSparams = zeros(numOfBS,3); %bc we have three params for a 2 exp fit, a, k1 & k2
    for i = 1:numOfBS
        [pStruc,~] = survivalPlotFromVector(bsMat,i,fit,fitV(1),fitV(2),[],fitV(3:end));
        BSparams(i,1) = pStruc.a;
        BSparams(i,2) = pStruc.k1;
        BSparams(i,3) = pStruc.k2;
    end
end

% if trhe fit sucks, get rid of it:
% BSparams(BSparams(:,1) > 0.999,:) = [];

mib.BSparams = BSparams;
mib.mean = mean(BSparams);
mib.std = std(BSparams);
mib.sem = std(BSparams)/sqrt(numOfBS);

%this is code copied over from bootstrapsBootstraps. May be useful if we
%decide to use this script to plot error on fits
%     
% sintsMat = zeros(mx + 1,numOfBS); %plus one bc we start from zero time
% for j = 1:numOfBS
%     sints = [];
%     for i = 0:mx %zero bc time
%         logik = bsMat(:,j) > i;
%         sints = [sints; sum(logik)];
%     end
%     sintsMat(:,j) = sints;
% end
% %so this is a bit tricky to wrap my mind around, but i'll try:
% % for each time (row) we will sort all the different BS vector values in
% % ascending order. then we can pluckl out the 5th & 95th percentile at each
% % time. make sense?
% B = sort(sintsMat,2);
% 
% %pluck out the 5th & 95th:
% lowerI = floor(0.05*numOfBS); %5th index
% if lowerI == 0
%     lowerI = 1;
% end
% upperI = ceil(0.95*numOfBS);
% lower = B(:,lowerI);
% upper = B(:,upperI);
% 
% %time vector:
% t = [0:mx]';
% 
% % plot them:
% figure(figNum);
% % plotHandle = stairs([t t],[log(lower/n) log(upper/n)],'LineWidth',1);
% % plotHandle = stairs([t t],[log(lower/n) log(upper/n)],'LineWidth',1);
% % add shading
% lower2 = nonzeros(lower);
% upper2 = nonzeros(upper); %must get rid of zeros to make shading work bc ln leads to -inf vals
% % now pad lower2 with 1's so its the same size as upper:
% pad = ones(size(upper2,1) - size(lower2,1),1);
% lower2 = [lower2; pad];
% bnd = size(upper2,1);
% lower2 = lower2';
% upper2 = upper2'; %bc they need to be row v's for some reason
% ts = [0:bnd - 1, fliplr(0:bnd - 1)]; %not sure why we gotta flip these things or if we really need fliplr and not just flip but whatever
% shade = [log(lower2/n), fliplr(log(upper2/n))];
% plotHandle = fill(ts,shade,colr); 
% set(plotHandle,'facealpha',0.25,'edgealpha',0.25,'edgecolor','none');
% 
% %% output
% mib = [lower upper];
































