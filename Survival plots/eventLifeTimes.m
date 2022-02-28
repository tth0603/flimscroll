function output = eventLifeTimes(cia,lowLengthThrd,figureNumber)
%Plug in the unedited cia, this returns a cia of purely 1 & 3 evnts and a survival plot
%of those events lifetimes on a linear time scale.
%output = eventLifeTimes(cia,lowLengthThrd,figureNumber)
%   cia: load CumulativeIntervalArray file, unedited
%   figureNumber: The fig number you are plotting in.
%   lowLengthThrd: threshold below which events will be excluded
%   output == (a structure with the following: 
%   [rateConstant average standardDeviation(sigma) numberOfEvents(n) intervals(a struc of the plotted event times) sints(a struc of the sorted event times)]
%   USAGE:
%   out = eventLifeTimes(cia,10,38);
logik=(cia(:,1)==-3)|(cia(:,1)==1);        %find the up events
cia2=cia(logik,:);

logik2=(cia2(:,5)>lowLengthThrd);
cia3=cia2(logik2,:);        %get rid of the short events, a cia-like array
    
m=max(cia3(:,5));

sints=[];           
for i=1:m;
    logik=cia3(:,5)>i;
    sints=[sints;i sum(logik)];
end

figure(figureNumber);
plot(sints(:,1),log(sints(:,2)/length(cia3)),'b');shg
output.avg=mean(cia3(:,5));
output.sigma=std(cia3(:,5))/sqrt(length(cia3));
output.n=length(cia3);
output.cia3=cia3;       %a cia-like array of up events
output.sints=sints;
end