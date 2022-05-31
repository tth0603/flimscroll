function out = fractionOfSigmaRetentionBS(sigmaCIA,oligoCIA,fps,numberOfStrapOns)
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
%USAGE:
%   re16outRetenBS = fractionOfSigmaRetentionBS(res16struc.sigmaCIA,res16struc.oligoCIA,res16struc.fps,1000);
n=length(sigmaCIA);
tbl=[[1:n]' 1/n*ones(n,1)];  

delT=oligoCIA(:,2)/fps-sigmaCIA(:,5);

amtRetainedVect=ones(numberOfStrapOns,1);
for i = 1:numberOfStrapOns
    indx=probability_steps(tbl,n);
    bs=delT(indx);
    logi=(bs<0);
    %bsAmtRetained=sum(logi);      %change this to 'bsAmtRetained=sum(logi)/n;' to make this a fraction 
    bsAmtRetained=sum(logi)/n;
    amtRetainedVect(i)=bsAmtRetained;
end
out.vector=amtRetainedVect;
out.stdMean=std(amtRetainedVect)/sqrt(numberOfStrapOns);
out.std=std(amtRetainedVect);
out.mean=mean(amtRetainedVect);
end
