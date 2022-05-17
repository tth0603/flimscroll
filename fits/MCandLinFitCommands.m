%here is making simulated data based on the rate of sigma photobleaching.
%We complare the sim data to the oligo distribution to see if there is an
%appreciable off rate from teh elongation complex
logi=(cia(:,1)==1|cia(:,1)==3);
oligocia=cia(logi,:);
oligoevents=res4struc.retention(:,7);
n=length(oligoevents);

% here is a hist for the distribution of transcript events:
bins=[0:25:800];
delTcounts=histc(oligoevents,bins);
figure(30);hold on;bar(bins,delTcounts/n,'histc');shg

 mc=oneStepRxnMC(0.0021,n);
 trxnSim=[oligoevents mc mc-oligoevents];
 logi2=trxnSim(:,3)>0;
 %percent retained:
 sum(logi2)/n
 % the first time I ran this, considering a PB off rate of 0.0024 (thats 0.0039 of the fit off rate at this exposure minus the esitmated yint of the pb off rate curve, 0.0015) 
 % 76% percent of the simulated datae events
 % outlived the oligo transcript events. This is higher than our observed
 % value of ~57%, but it's in the same ballpark. It's not unreasonable
 
 %let's bootstrap the above:
 simFrac=zeros(100,1);
 for i = 1:100
  mc=oneStepRxnMC(0.0021,n);
 trxnSim=mc-oligoevents;
 logi2=trxnSim(:)>0;
 simFrac(i)=sum(logi2)/n;
 end
 simFracMean=mean(simFrac)
 simFracStd=std(simFrac)
 %plot(fracExp(2),simFracMean,'ro');
 errorbar(fracExp(2),simFracMean,simFracStd,'ko');
 
 %fitting to a linear regrssion:
fracExp=[0.9591 0.1149 0.0635];      %x values
rates=[0.0176 0.0039 0.0027];       %y values
[p,errEst] = polyfit(fracExp,rates,1);
fit=polyval(p,fracExp,errEst);
figure(10);plot(fracExp,fit,'-',fracExp,rates,'*');

rates=[0.0027 0.0039 0.0176];       %from single exponential fits to sigma off rate at various exposures [16s 4s 1s]
fracExp=[0.9591 0.1149 0.0635];
%for res 1s (continuous observation):
res1k=rates(3)-0.0018;
res1oligoRetd=res1struc.retention(:,7);
n=length(res1oligoRetd);
 simFrac=zeros(100,1);
 for i = 1:100
  mc=oneStepRxnMC(res1k,n);
 trxnSim=mc-res1oligoRetd;
 logi2=trxnSim(:)>0;
 simFrac(i)=sum(logi2)/n;
 end
 res1simFracMean=mean(simFrac)
 res1simFracStd=std(simFrac)
 %plot(fracExp(2),simFracMean,'ro');
 errorbar(fracExp(1),res1simFracMean,res1simFracStd,'ko');
 
 %for res 16s (continuous observation):
res1k=rates(1)-0.0018;
res16oligoRetd=res16struc.retention(:,7);
n=length(res16oligoRetd);
 simFrac=zeros(100,1);
 for i = 1:100
  mc=oneStepRxnMC(res16k,n);
 trxnSim=mc-res16oligoRetd;
 logi2=trxnSim(:)>0;
 simFrac(i)=sum(logi2)/n;
 end
 res16simFracMean=mean(simFrac)
 res16simFracStd=std(simFrac)
 %plot(fracExp(2),simFracMean,'ro');
 errorbar(fracExp(3),res16simFracMean,res16simFracStd,'ko');
