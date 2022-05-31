function mib = bootstrapFromVector(vector,numberOfStrapOns)
% mib = bootstrapFromVector(vector,numberOfStrapOns)
%
%   This will find the error on a single exponential fit through boot
%   strapping. it takes a while to run, so keep boot straps to 100 or so.
%   You may have to screw with fit params below.
%
%   INPUTS
%   vector == a vector with the points to be included in the bootstrap
%   numvberOfStrapOns == a scaler, the number of BSs
%   OUTPUTS
%   a structure containing: bs - a vector of the fit rate for each bs;
%   stdRate - the stdv for the fitted rates; stdTime - the stdv of the fit
%   time constants; bootstraps - a mat of all teh bootstrap vectors as
%   columns.
%   
% Harden 2015
sigmaLifeTime=vector;     %gives duration of event in time. This works for the '-3' events that are in out.sigmaCIA, but may need to be edited
n=length(vector);
tbl=[[1:n]' 1/n*ones(n,1)];  

%fit params
tau=0.1;           %screw with this one
tm=1; %dont make this much greater than one
tx=max(sigmaLifeTime); %leave this untouvhed

bsFitRate=ones(numberOfStrapOns,1);
tic;
bsToPlot = zeros(n,numberOfStrapOns);
for i = 1:numberOfStrapOns
    indx=probability_steps(tbl,n);
    bs=sigmaLifeTime(indx); 
    bsFitRate(i)=abs(fminsearch('expfallone_mxl',tau,[],bs,tm,tx));
    bsToPlot(:,i) = bs;
end
toc

mib.bs=bsFitRate;
mib.stdRate=std(bsFitRate);

bsFitTimeConst=1./bsFitRate;
mib.stdTime=std(bsFitTimeConst);
mib.bootstraps=bsToPlot;
end
