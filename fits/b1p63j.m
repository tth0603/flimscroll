load t:\matlab\data\b1p62c.dat -mat %Interval Data File, MUT DNA, 1 EXP fit

ciamut=Intervals.CumulativeIntervalArray; %inspect Intervals variable using 
                                       %Intervals.CumulativeIntervalArrayDescription
logik=ciamut(:,1)==1;

intervals=ciamut(logik,4);
disp('Sample Intervals:');
disp(intervals(1:10));  %Defines the intervals data to be plotted

%lg=intervals<200;
%intervals=intervals(lg);         
                            %removes two outliers above 300.  May want to 
                            %play with this. Changes the lenth of
                            %'sints' from 290 to 288.

ap=1;
a=1/(1+ap^2);
tau=1/10;    %tau refers to rate constants, change to time constants by    
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

tfitmut1=0:0.5:200;

SPfitmut1=exp(-tfitmut1*k);  %This is the survival curve model
                                                     
SPx=exp(-tx*k);  
SPm=exp(-tm*k);  %survival probability evaluated at tm and tx

for i=1:length(sints)
    SPmut(i)=(1-i/length(sints))*(SPm-SPx)+SPx;    
end
                                       %This is what the data set will be 
                                      %plotted against.  Instead of 
                                      %straight time, it is plotted against 
                                      %the time intervals of events, then 
                                      %scaled to account for the fact that
                                      %we are describing the entire
                                      %phenomonon, but only observing a
                                      %small window of time 


figure(6);plot(tfitmut1,log(SPfitmut1),sints,log(SPmut));shg