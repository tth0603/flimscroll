function mib = bootstrapExp1mxlAll(dwellts,obsts,numberOfStrapOns)
%   This will find the error on a single exponential fit to data where some
%   of the observations have been truncated by, say, the length of the
%   observation window. Used to find errors on SMfig3 fo sigma paper
%   INPUTS:
%   vector == a vector with the points to be included in the bootstrap
%   
n1 = length(dwellts);
n2 = length(obsts);
newdwellts = [dwellts zeros(n1,1)];
newobsts = [obsts ones(n2,1)];
sigmaLifeTime=[newdwellts;newobsts];     %gives duration of event in time. This works for the '-3' events that are in out.sigmaCIA, but may need to be edited
n=(n1+n2);
tbl=[[1:n]' 1/n*ones(n,1)];  

%fit params
tau=200;        %this fit param is in units of s, not s^-1     
tm=16;

bsFitRate=ones(numberOfStrapOns,1);
bsToPlot=[];
for i = 1:numberOfStrapOns;
    indx=probability_steps(tbl,n);
    bs=sigmaLifeTime(indx,1:2);
    for j = 1:n
        if bs(j,2)==1
            bsobsts(j) = bs(j,1);
        else
            bsdwellts(j) = bs(j,1);
        end
    end 
    bsFitRate(i)=fminsearch('exp1_mxlall',tau,[],bsdwellts,bsobsts,tm);
    %bootstraps=[bootstraps;bs];     %these are the values plucked from the experimental data daisy chained together
    %bsToPlot = [bsToPlot bs];
end
mib.bs=bsFitRate;
mib.stdTime=std(bsFitRate);

bsFitTimeConst=1./bsFitRate;
mib.stdRate=std(bsFitTimeConst);
%out.bootstraps=bsToPlot;
end