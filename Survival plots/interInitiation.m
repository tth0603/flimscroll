function output = interInitiation(cia,figureNumber,fit)
%output = interInitiation(cia,figureNumber,fit)
%this function plots the time between the start of successive 'hi' events.
%   cia: load CumulativeIntervalArray file
%   figureNumber: The fig number you are plotting in.
%   fit: ***OPTIONAL*** if fit == 1, a single exponential, if fit == 2, a double
%   exponential, otherwise, no fit.
%   output == (a structure with the following: 
%   [rateConstant timeConstant average standardDeviation(sigma) numberOfEvents(n) intervals(a struc of the plotted event times) sints(a struc of the sorted event times)]
%   
if nargin==2        %This allows the optional input
    fit=[];
end

logik=(cia(:,1)==0)|(cia(:,1)==1);      %here we cook up our list of intervals between the start of successive '1' events.  i.e. 'inter-initiation.'

data=[cia(logik,5) cia(logik,7) cia(logik,1)];

L=length(data);

ints=[];
for i=1:L-1;
    if (data(i,2)==data(i+1,2)) && (data(i,3)==1);
    ints(i)=data(i,1)+data(i+1,1);
    end
end
ints=ints(ints~=0);

intervals=ints';

output.avg=mean(intervals);
output.sigma=std(intervals);
output.n=length(intervals);

%%%Starts 3 part if statement to address fitting with a single, double, or
%%%no fit:
%%%For a single exp fit:
if fit==1       %an example of a proper way to fit survival curves with a single exponential max likelihood
      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
tau=0.001;    %tau refers to rate constants, change to time constants by    
             %edit expfalltwo_mxl       
tm=10;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
tx=max(intervals);      %tm is about the frame length, of the lower event threshold 
                        %tx is more of less the longest possible observable
                        %event.

params=fminsearch('expfallone_mxl',tau,[],intervals,tm,tx);
corr=expfallone_mxl_correction(tau,intervals,tm,tx);

rateConst=params(1);  %gives rate constants for fit
timeConst=1./rateConst;

tfit=0:max(intervals);
SPfit=exp(-tfit*rateConst);

sints=[];
for i=1:max(intervals);
logik=intervals>i;
sints=[sints;i sum(logik)];
end

figure(figureNumber);
plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'b',tfit,log(SPfit),'c');shg   %The intervals plotted here are corrected for finite observation time
output.rateConst=rateConst;
output.timeConst=timeConst;
output.avg=mean(intervals);
output.sigma=std(intervals);
output.n=n;
output.intervals=intervals;

%%%Part 2 of if statement:
elseif fit==2       %an example of a proper way to fit survival curves with a double exponential max likelihood
    ap=2.3;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
    tau1=0.003;
    tau2=0.0006;

    tm=10;      %%%%%%%%%%%%%%%You will need to fuck wiht these numbers%%%%%%%%%%%%%%%%%%
    tx=max(intervals);
    inarg=[ap tau1 tau2];
    
    params=fminsearch('expfalltwo_mxl',inarg,[],intervals,tm,tx);   %calcs the rates and amplitudes of the mxl fit
%%%To correct for a finite observation time, uncomment 'corr' and comment
%%%'corr.Nzero' and 'corr.offset'%%%
    %corr=expfalltwo_mxl_correction(inarg,intervals,tm,tx);      %correction factors for the fact that want to characterize the entire phenomena, but only observe it for a part of the time
        %'corr' is a struc: 
        % corr.CorrectionFactor       % = (  1/( a*A + (1-a)*B ) 
        % corr.Nzero                  % = actual # of events in data set (observed+not observed)  
        % corr.Offset            % # remaining events at time = tx
        
ap=params(1);
a=1/(1+ap^2);
k1=params(2);
k2=params(3);  %defines rate constants for fit

tfit=[0:max(intervals)];

SPfit=a*exp(-tfit*k1)+(1-a)*exp(-tfit*k2);  %This is the survival 
                                            %curve model
sints=[];
for i=1:max(intervals);
logik=intervals>i;
sints=[sints;i sum(logik)];
end

corr.Offset=0;
corr.Nzero=max(sints(:,2));

figure(figureNumber);
plot(sints(:,1),log((sints(:,2)+corr.Offset)/corr.Nzero),'b',tfit,log(SPfit),'c');shg
output.amplitude=a;
output.k1=params(2);
output.k2=params(3);
output.intervals=intervals;
output.sints=sints;

%%%Part 3 of if statement (no fit):
else
    sints=[];
for i=1:max(intervals);
logik=intervals>i;
sints=[sints;i sum(logik)];
end

avg=mean(sints);
sigma=std(sints);
l=max(sints(:,2));

figure(figureNumber);
plot(sints(:,1),log(sints(:,2)/l),'b');shg
output.avg=avg;
output.sigma=sigma;
output.n=l;
output.intervals=intervals;
output.sints=sints;
end


end
