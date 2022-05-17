%function simSigRetenData(n,appK,probeLifeTimes)
%number of observations:

n=100;
%simulated apparent off rate:
%appK=0.005;
appK=1/212;
releasedMat=[];
retainedMat=[];
eventMat=[];
for i = 1:2
%simulate event intervals (analogous to sigma lifetimes):
simInts=oneStepRxnMC(appK,n);
%for observation times plucked from the experimental data:
probeLifetimes = res16strucV3.oligoCIA(:,5);
%probeLifetimes = res16strucV3.retention(:,7);
L = length(probeLifetimes);
obsTsStruc = bootstrapFromVector(probeLifetimes,ceil(n/L));
obsTs = obsTsStruc.bootstraps(1:n);

%a simulated data matrix. columns are: 
%1:observation times (obsTs, ie probelifetimes) 
%2:event intervals (ie sigma lifetimes) 
%3:timeSigmaDeparture OR timeProbeDepart, whichever is less 
%4:'1' if timeSigmaDepart>timeProbeDepart OR '0' if timeSigmaDepart<timeProbeDepart 
simData=[obsTs simInts];
logi1=simData(:,1)<simData(:,2);
for i=1:n
    if simData(i,1)<=simData(i,2);
        simData(i,3)=simData(i,1);
    else simData(i,3)=simData(i,2);
    end
end
simData=[simData logi1];    %this gives col 4
        
logi4=simData(:,4)==0;
dwellts = simData(logi4,3);
releasedEvents=dwellts;
dwell_rec = simData(logi4,1);
logi3=simData(:,4)==1;
obsts = simData(logi3,3);
retainedEvents=obsts;

%maybe this is wrong. kMaybe you want to bin each BS using 'bins' and
%'histc' then include them in a mat.... In fact you need to do this,
%otherwise you cannot concatonate the matrices
releasedMat=[releasedMat releasedEvents];
retainedMat=[retainedMat releasedEvents];
eventMat = [eventMat releasedEvents;releasedEvents];
end