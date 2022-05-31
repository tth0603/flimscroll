function out = sigmaBootstrapon(WTdata,MUTdata,numberOfStrapOns,figNum)
%Takes in tweo binary data sets for sigma retention (1) or release (0), then
%bootstraps the two sets in pairs and asks which (the WT or MUT data) yield
%a higher ratio of retention
%   WTdata/MUTdata == nx2 matrix of elongation events. Create by running 'sigmaRetenDistribution'
%   [expl: stat=sigmaRetenDistribution(mut2sigcia,mutoligocia,[],1,52)]
%   This returns a strucure with member 'delT,' a 2xn mat. delT is the
%   input
%USAGE:
%       out = sigmaBootstrapon(WTdata,MUTdata,numberOfStrapOns)
WTn=length(WTdata);
WTtbl=[[1:WTn]' 1/WTn*ones(WTn,1)];  

MUTn=length(MUTdata);
MUTtbl=[[1:MUTn]' 1/MUTn*ones(MUTn,1)]; 

WTbsMat=[];
MUTbsMat=[];
tic
for i = 1:numberOfStrapOns;
%tic;
WTindx=probability_steps(WTtbl,WTn);   
WTbs=WTdata(WTindx);                %one bootstrap, I think. This should be the same as delT from sigmaRetenDistribution
bins = [-1000:100:2500];
WTnormalization=histc(WTbs(:,1),bins)/WTn;     %normalized counts to be plotted in PDF
WTbsMat=[WTbsMat;WTnormalization'];

MUTindx=probability_steps(MUTtbl,MUTn);   
MUTbs=MUTdata(MUTindx);    
MUTcounts=histc(MUTbs(:,1),bins)/MUTn;
MUTbsMat=[MUTbsMat;MUTcounts'];
end
toc;
out.WTbsMat=WTbsMat;
% 
% WTavgCounts=mean(WTbsMat);
% WTxLower = [-1000 -500 -200 -100:100:2500];
% WTxUpper = [-500 -200 -100:100:2600];
% WTy = [mean(WTavgCounts(1:5)) mean(WTavgCounts(6:8)) WTavgCounts(9:36)];
% hold on;barGraph(WTxLower,WTxUpper,WTy,'r',figNum)
% 
% MUTavgCounts=mean(MUTbsMat);
% MUTxLower = [-1000 -500 -100:100:2500];
% MUTxUpper = [-500 -100:100:2600];
% MUTy = [mean(MUTavgCounts(1:5)) mean(MUTavgCounts(6:9)) MUTavgCounts(10:36)];
% barGraph(MUTxLower,MUTxUpper,MUTy','k',figNum)

%to put error bars on the original plots results:
WTstd=std(WTbsMat);
WTxLower = [-1000 -500 -200 -100:100:2500];
WTxUpper = [-500 -200 -100:100:2600];
WTerror = [mean(WTstd(1:5)) mean(WTstd(6:8)) WTstd(9:36)];
WTnormalization=histc(WTdata(:,1),bins)/WTn;
WTy = [mean(WTnormalization(1:5))' mean(WTnormalization(6:8))' WTnormalization(9:36)'];
WTyUpper = WTy + 2*WTerror;
WTyLower = WTy - 2*WTerror;
%hold on;barGraph(WTxLower,WTxUpper,WTyUpper,'r',figNum);barGraph(WTxLower,WTxUpper,WTyLower,'b',figNum)
WTmidX=(WTxLower+WTxUpper)/2;
hold on;figure(figNum);scatter(WTmidX,WTyUpper,'.','r');
hold on;figure(figNum);scatter(WTmidX,WTyLower,'.','r');


MUTstd=std(MUTbsMat);
out.MUTstd=MUTstd;
MUTxLower = [-1000 -500 -100:100:2500];
MUTxUpper = [-500 -100:100:2600];
MUTerror = [mean(MUTstd(1:5)) mean(MUTstd(6:9)) MUTstd(10:36)];
MUTnormalization=histc(MUTdata(:,1),bins)/MUTn;
MUTy= [mean(MUTnormalization(1:5))' mean(MUTnormalization(6:9))' MUTnormalization(10:36)'];
MUTyUpper = MUTy + 2*MUTerror;
MUTyLower = MUTy - 2*MUTerror;
MUTmidX=(MUTxLower+MUTxUpper)/2;
hold on;figure(figNum);scatter(MUTmidX,MUTyUpper,'.','k');
hold on;figure(figNum);scatter(MUTmidX,MUTyLower,'.','k');