load C:/matlab/data/b1p121c.dat -mat %Interval Data File, MUT DNA, 1 EXP fit

cia=Intervals.CumulativeIntervalArray; %inspect Intervals variable using 
                                       %Intervals.CumulativeIntervalArrayDescription
logik=cia(:,1)==0;

intervals=(cia(logik,5));
L=length(intervals);
%disp('Sample Intervals:');
%disp(intervals(1:10));  %Defines the intervals data to be plotted

%lg=intervals<200;
%intervals=intervals(lg);         
                            %removes two outliers above 300.  May want to 
                            %play with this. Changes the lenth of
                            %'sints' from 290 to 288.

ap=1;
a=1/(1+ap^2);
tau=0.01;    %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=0.5;
tx=max(intervals);      %tm is (1/2)*frame length, tm is the longest event 
                      %in the data set

params=fminsearch('expfallone_mxl',tau,[],intervals,tm,tx);

k=params(1);  %gives rate constants for fit
disp('Rate Constant (1/s):');
disp(k);
t=1./k;
disp('Time Constant (s):');
disp(t);

sints=sort(intervals);  %this sorts the intervals data by length

tfit=0:0.5:max(intervals);

SPfit=exp(-tfit*k);  %This is the survival curve model
                                                     
SPx=exp(-tx*k);  
SPm=exp(-tm*k);  %survival probability evaluated at tm and tx

%for i=1:length(sints)
    %SPmut(i)=(1-i/length(sints))*(SPm-SPx)+SPx;    
%end
                                       %This is what the data set will be 
                                      %plotted against.  Instead of 
                                      %straight time, it is plotted against 
                                      %the time intervals of events, then 
                                      %scaled to account for the fact that
                                      %we are describing the entire
                                      %phenomonon, but only observing a
                                     %small window of time 
sints=[];
for i=0:max(intervals)
logik=intervals>i;
sints=[sints;i sum(logik)];
end
figure(8);plot(tfit,log(SPfit),sints(:,1),log(sints(:,2)/L));shg

%figure(6);plot(tfit,log(SPfit),sints,log(SPmut));shg