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
tx=3600;
inarg=[ap tau1 tau2];  %tm is (1/2)*frame length, tm is the longest event 
                        %in the data set

params=fminsearch('expfalltwo_mxl',inarg,[],intervals,tm,tx);

ap=params(1);
a=1/(1+ap^2)
k1=params(2)
k2=params(3)  %defines rate constants for fit

tfitmut2=0:0.5:900;

SPfitmut2=a*exp(-tfitmut2*k1)+(1-a)*exp(-tfitmut2*k2);  %This is the survival 
                                            %curve model
sintsmut2=[];                                                     
for indx = 0:886
    logik=intervals>=indx;
    sintsmut2=[sintsmut2;indx sum(logik)];
end
    
figure(4);plot(tfitmut2,log(SPfitmut2),'r',sintsmut2(:,1),log(sintsmut2(:,2)/633),'g');shg

%To plot this simultaneously with the Mutant curve and fit, uncomment the 
%below command and run this script along with b1p63i.m

%figure(10);plot(tfitmut2,log(SPfitmut2),sintsmut2,log(SPmut2),tfit,log(SPf
%it),sints,log(SP));shg