%WT survival curve ploit with a two EXP fit

load t:\matlab\data\b1p63a.dat -mat %Interval Data File for control, 
                                    
cia=Intervals.CumulativeIntervalArray; %inspect Intervals variable using 
                                       %Intervals.CumulativeIntervalArrayDescription
logik=cia(:,1)==1;

intervals=cia(logik,4);
%disp('Sample Intervals:');
%disp(intervals(1:10));  %Defines the intervals data to be plotted


%lg=intervals<400;
%intervals=intervals(lg);         %removes single outlier at 617.  May want to 
                            %undue this
ap=1;
tau1=1/10;
tau2=1/100;  %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=0.5;
tx=max(intervals);
inarg=[ap tau1 tau2];  %tm is (1/2)*frame length, tm is the longest event 
                        %in the data set

params=fminsearch('expfalltwo_mxl',inarg,[],intervals,tm,tx);

ap=params(1);
a=1/(1+ap^2)
k1=params(2)
k2=params(3)  %defines rate constants for fit

sints=sort(intervals);  %this sorts the intervals data by length

tfit=0:0.5:620;

SPfit=a*exp(-tfit*k1)+(1-a)*exp(-tfit*k2);  %This is the survival 
                                            %curve model
                                                     
SPx=a*exp(-tx*k1)+(1-a)*exp(-tx*k2);  
SPm=a*exp(-tm*k1)+(1-a)*exp(-tm*k2);  %survival probability evaluated at 
                                      % tm and tx

for i=1:478
    SP(i)=(1-i/478)*(SPm-SPx)+SPx;    
end
                                       %This is what the data set will be 
                                      %plotted against.  Instead of 
                                      %simply time, it is plotted against 
                                      %the time intervals of events, then 
                                      %scaled to account for the fact that
                                      %we are describing the entire
                                      %phenomonon, but are limited by time
                                      %resolution and observation time. 


figure(4);plot(tfit,log(SPfit));
hold on;stairs(sints,log(SP),'r');shg
%To plot this simultaneously with the control curve and fit, uncomment the 
%below command and run this script along with b1p63k.m

%figure(10);plot(tfitmut2,log(SPfitmut2),sintsmut2,log(SPmut2),tfit,log(SPf
%it),sints,log(SP));shg

