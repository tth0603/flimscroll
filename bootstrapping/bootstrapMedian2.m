function out = bootstrapMedian2(vector,numberOfStrapOns)
%   This will find the error on theh median of an distribution
%   INPUTS:
%   vector == a vector with the points to be included in the bootstrap
%   
%sigmaLifeTime=vector;     %gives duration of event in time. This works for the '-3' events that are in out.sigmaCIA, but may need to be edited
n=length(vector);
tbl=[[1:n]' 1/n*ones(n,1)];  

bsMat=zeros(n,numberOfStrapOns);
for i = 1:numberOfStrapOns;
    indx=probability_steps(tbl,n);
    bs=vector(indx); 
    bsMat(:,i)=bs;
    %bootstraps=[bootstraps;bs];     %these are the values plucked from the experimental data daisy chained together
    %bsToPlot = [bsToPlot bs];
end
medVector = median(bsMat);
out.test = medVector;
out.averageMedian = mean(medVector);
out.medianStd = std(medVector);

    
% out.medians=bsMedian;
% out.averageMedian=mean(bsMedian);
% out.std=std(bsMedian);

% bsFitTimeConst=1./bsFitRate;
% out.stdTime=std(bsFitTimeConst);
% out.bootstraps=bsToPlot;
end
