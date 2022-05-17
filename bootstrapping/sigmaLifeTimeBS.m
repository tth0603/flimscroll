function out = sigmaLifeTimeBS(cia,numberOfStrapOns)
%   INPUTS:
%   cia == a cia like array, usu from 'sigmaRetenDistribution' function.
%   use 'out.sigmaCIA' tro ensure you are using sigma events that preceed
%   a trxn event.
sigmaLifeTime=cia(:,5);     %gives duration of event in time. This works for the '-3' events that are in out.sigmaCIA, but may need to be edited
n=length(cia);
tbl=[[1:n]' 1/n*ones(n,1)];  

%fit params
tau=0.1;           
tm=1;
tx=max(sigmaLifeTime);

bsFitRate=ones(numberOfStrapOns,1);
for i = 1:numberOfStrapOns;
    indx=probability_steps(tbl,n);
    bs=sigmaLifeTime(indx); 
    bsFitRate(i)=abs(fminsearch('expfallone_mxl',tau,[],bs,tm,tx));
end
out.bs=bsFitRate;
out.stdRate=std(bsFitRate);

bsFitTimeConst=1./bsFitRate;
out.stdTime=std(bsFitTimeConst);
end
