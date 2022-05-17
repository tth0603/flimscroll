function out = sigmaBinaryBootstrapon(WTdata,MUTdata,numberOfStrapOns)
%Takes in tweo binary data sets for sigma retention (1) or release (0), then
%bootstraps the two sets in pairs and asks which (the WT or MUT data) yield
%a higher ratio of retention
%   WTdata/MUTdata == a binary row vector. 'n' is the number of elongation
%   events in the data. create by running 'sigmaRetenDistribution'
%   [expl: stat=sigmaRetenDistribution(mut2sigcia,mutoligocia,[],1,52)]
%   This returns a strucure with member 'delT,' a 2xn mat. you can then
%   turn this matrix into a binary matrix by:
%       WTdata = stat.delT(:,1)<0;
%USAGE:
%       out = sigmaBootstrapon(WTdata,MUTdata,numberOfStrapOns);

diffVect=[];
WTn=length(WTdata);
WTtbl=[[1:WTn]' 1/WTn*ones(WTn,1)];  

MUTn=length(MUTdata);
MUTtbl=[[1:MUTn]' 1/MUTn*ones(MUTn,1)]; 
tic
for i = 1:numberOfStrapOns;
%tic;
WTindx=probability_steps(WTtbl,WTn);   
WTtstint=WTdata(WTindx);    
WTs=sum(WTtstint);
WTratio=WTs/WTn;

MUTindx=probability_steps(MUTtbl,MUTn);   
MUTtstint=MUTdata(MUTindx);    
MUTs=sum(MUTtstint);
MUTratio=MUTs/MUTn;

digg=WTratio-MUTratio;
diffLogi=digg>0;
diffVect=[diffVect;diffLogi];
%toc;
end
toc;
out.tic=toc;
out.diffVect=diffVect;
out.WTfrac=sum(diffVect)/numberOfStrapOns;


% argouts=zeros(numberOfStrapOns,3);             
% 
% for indx=1:numberOfStrapOns                    
% ind=probability_steps(tbl,n);    %rate constants
% WTtstint=intervals(ind);
% argouts(indx,:)=fminsearch('expfalltwo_mxl',[0.85 0.07 0.01],[],WTtstint,0.5,2500);
% if indx/50==round(indx/50)         
% indx;                               
                                           