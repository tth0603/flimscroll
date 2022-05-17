%This will give boot strapped data for the MUT 2 EXP fits

load t:\matlab\data\b1p62c.dat -mat

ciamut=Intervals.CumulativeIntervalArray;
logik=ciamut(:,1)==1;
intervals=ciamut(logik,4);

%help probability_steps          %This is the Bootstrap function

tbl=[[1:633]' 1/633*ones(633,1)];  %Makes a 2xN array of numbers with weights 1/N.
indx=probability_steps(tbl,633);   %'N' is the number of events. Here N is 633. 'indx' draws N numbers from 'tbl'
tstint=intervals(indx);            %with weight 1/N. tstint is an example of a single bootstrapped data set.

argouts=zeros(100,3);             %allocates space for 'argouts'

for indx=1:100                    %generates 1000 bootstrapped data sets, and the corresponding observable
ind=probability_steps(tbl,633);    %rate constants
tstint=intervals(ind);
argouts(indx,:)=fminsearch('expfalltwo_mxl',[0.85 0.07 0.01],[],tstint,0.5,2500);
if indx/50==round(indx/50)         %not sure why Larry threw this in, perhaps I should play with it a bit
indx;                               %LAST I TRIED THIS IS DID NOT WORK!!! works when probability in 'tbl' is 1/290, but not 1/633 for some reason  
end
end

%[bootstat bootsam]=bootstrp(100,[],intervals);  %creates a n x m matrix.  
                                            %Where n is the sample size and
                                            %m is the number of samples.
                                            %Here m is 100
%intervals=bootsam;