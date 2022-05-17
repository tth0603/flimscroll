%sigma off rates from ECs (Considering all retention events either isgma or probe departure, whichever is
%less [1s 1/8s 1/16s]:
rates=[0.02 0.0123 0.0117];
%fps:
ts=[0.9591 0.1149 0.0635];
figure(35);plot(ts,rates,'r.',ts,fit,'-');

[p,errEst] = polyfit(ts,rates,1);
fit=polyval(p,ts,errEst);
figure(36);plot(ts,fit,'-',ts,rates,'*');
%to plot the errors (calc'd using singleExpRateBootstrap.m)
err=[0.0078 0.0032 0.0021];
hold on
figure(37);errorbar(ts,rates,err);



%sigma off rates from ECs (only those events where sigma leaves before probe) [1s 1/8s 1/16s]:
rates=[0.0262 0.0181 0.0086]; %these are from the reciprical of the mean given by survivalPlotFromVector because the exp fits are wonky ands need to be revised, the rates from the fits are: [0.0148 0.0110 0.0085]
%fps:
ts=[0.9591 0.1149 0.0635];
figure(45);plot(ts,rates,'r.',ts,fit,'-');

[p,errEst] = polyfit(ts,rates,1);
fit=polyval(p,ts,errEst);
figure(45);plot(ts,fit,'-',ts,rates,'*');
%to plot the errors (calc'd using singleExpRateBootstrap.m)
err=[0.0078 0.0032 0.0021];
hold on
figure(45);errorbar(ts,rates,err,'r.');

fracRetained=[19/154 26/157 54/193];
ts=[0.9591 0.1149 0.0635];
inFracRetained=1./fracRetained;
[p,errEst] = polyfit(ts,inFracRetained,1)
%p =[-4.1154 7.4664] %(slope, y int);
t=0:0.01:1;
fracRetainedfit=1./(p(1)*t+p(2));
figure(91);plot(t,fracRetainedfit);shg

fracRetainedfit=polyval(p,ts,errEst);
figure(91);plot(ts,fracRetainedfit,'-',ts,inFracRetained,'*');

fracRetainedThroughout=[0.263 0.577 0.5741]; 
ts=[0.9591 0.1149 0.0635];
invFracRetainedThroughout=1./fracRetainedThroughout;
[p,errEst] = polyfit(ts,invFracRetainedThroughout,1)
t=0:0.01:1;
fracRetainedfit=1./(p(1)*t+p(2));
figure(91);plot(t,fracRetainedfit);shg
