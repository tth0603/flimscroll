%This will give boot strapped data for the MUT 2 EXP fits

load t:\matlab\data\b1p63a.dat -mat

cia=Intervals.CumulativeIntervalArray;
logik=cia(:,1)==1;
intervals=cia(logik,4);

%help probability_steps          %This is the Bootstrap function

wttbl=[[1:478]' 1/478*ones(478,1)];  %Makes a 2xN array of numbers with weights 1/N.
indx=probability_steps(wttbl,478);   %Here N is 633. 'indx draws N numbers from 'tbl'
tstint=intervals(indx);            %with weight 1/N. tstint is an example of a single bootstrapped data set.

argouts=zeros(100,3);             %allocates space for 'argouts'

for indx=1:100                    %generates 1000 bootstrapped data sets, and the corresponding observable
ind=probability_steps(wttbl,478);    %rate constants
tstint=intervals(ind);
argouts(indx,:)=fminsearch('expfalltwo_mxl',[0.85 0.07 0.01],[],tstint,0.5,2500);
if indx/50==round(indx/50)         %not sure why Larry through this in, perhaps I should play with it a bit
indx;                               % 
end
end
