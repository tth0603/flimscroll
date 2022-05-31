function out = fracRetainedThroughoutBS(retentionMat,numberOfStrapOns)
%   This bootstrapping routine is to be used in conjunction with
%   sigmaRetentionDistribution.m
%   INPUTS:
%   sigmaCIA == cia-like array from above .m file output.
%   oligoCIA == same same but different
%   fps == scaler from above .m
%   retentionMat == a mat from 'sigmaRetenDistribution' function.
%   use 'out.retention.' See the help file of sigmaRetentionDistribution.m
%   for member details.
%   OUTPUT
%   .vector = vector of bootstrapped percent retained
%   .std = error on amt retained
%   .mean = mean percent retained
%   .mean = mean percent retained
%USAGE:
%   re16RetenBSthroughout = fracRetainedThroughoutBS(res16struc.retention,10);
n=length(retentionMat);
tbl=[[1:n]' 1/n*ones(n,1)];  

fracRetainedThroughout=ones(numberOfStrapOns,1);
for i = 1:numberOfStrapOns
    indx=probability_steps(tbl,n);
    oligoDept=retentionMat(indx,8);     % Since we are only looking at retained sigma events, we only care about the relative sigma and probe departure
    sigmaDept=retentionMat(indx,9);
    releaseVect=oligoDept-sigmaDept; 
    logi=(releaseVect<=0);      % If oligo depart < sigma depart sigma is retained and this returns affirmative logical
    retention=releaseVect(logi);   
    fracRetainedThroughout(i)=sum(logi)/length(releaseVect);
end
out.oligoDept=oligoDept;
out.releaseVect=releaseVect;
out.logi=logi;
out.vector=fracRetainedThroughout;
out.std=std(fracRetainedThroughout);
out.mean=mean(fracRetainedThroughout);
end
