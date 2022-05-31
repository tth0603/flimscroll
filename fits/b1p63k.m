%MUT DNA survival plot with a 2 EXP fit

load t:\matlab\data\b1p62c.dat -mat %Interval Data File, 

ciamut2=Intervals.CumulativeIntervalArray; %inspect Intervals variable using 
                                       %Intervals.CumulativeIntervalArrayDescription
logik=ciamut2(:,1)==1;

intervals=ciamut2(logik,4);
%disp('Sample Intervals:');
%disp(intervals(1:10));  %Defines the intervals data to be plotted


%lg=intervals<420;
%intervals=intervals(lg);        %removes two outliers above 500.  May want to 
                            %play with this. Changes the lenth of
                            %'sints' from 633 to 630.

ap=1;
tau1=1/10;
tau2=1/100;  %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=0.5;
%tx=max(intervals);
tx=3600;
inarg=[ap tau1 tau2];  %tm is (1/2)*frame length, tm is the longest event 
                        %in the data set

params=fminsearch('expfalltwo_mxl',inarg,[],intervals,tm,tx);

ap=params(1);
a=1/(1+ap^2)
k1=params(2)
k2=params(3)  %defines rate constants for fit

sintsmut2=sort(intervals);  %this sorts the intervals data by length

tfitmut2=0:0.5:900;

SPfitmut2=a*exp(-tfitmut2*k1)+(1-a)*exp(-tfitmut2*k2);  %This is the survival 
                                            %curve model
                                                     
SPx=a*exp(-tx*k1)+(1-a)*exp(-tx*k2);  
SPm=a*exp(-tm*k1)+(1-a)*exp(-tm*k2);  %survival probability evaluated at 
                                      % tm and tx

for i=1:633
    SPmut2(i)=(1-i/633)*(SPm-SPx)+SPx;    
end
                                       %This is what the data set will be 
                                      %plotted against.  Instead of 
                                      %simply time, it is plotted against 
                                      %the time intervals of events, then 
                                      %scaled to account for the fact that
                                      %we are describing the entire
                                      %phenomonon, but are limited by time
                                      %resolution and observation time. 
figure(9);plot(tfitmut2,log(SPfitmut2));
hold on;stairs(sintsmut2,log(SPmut2),'g');shg
%To plot this simultaneously with the Mutant curve and fit, uncomment the 
%below command and run this script along with b1p63i.m

%figure(10);plot(tfitmut2,log(SPfitmut2),sintsmut2,log(SPmut2),tfit,log(SPf
%it),sints,log(SP));shg