function [mih,dataHandle] = firstPassagePlot(data,col,tMx,normMx,figNum)
% [mib,dataHandle,fitHandle] = firstPassagePlot(data,col,tMx,normMx,figNum)
% 
% INPUTS
%   data == a 2D matrix or vector
%   col == the column within DATA to be plotted
%   tMx == the maximum time of the experiment
%   normMx == what to normalize the plot to. e.g. the number of nuclei. set
%       equal to 1 if unnormalized
%   figNum == a scaler

%grab the first event
first = data(:,col); %a vector

%cumulative intervals:
sints = [];
for i = 0:ceil(tMx) + 1 %here we take the length of the experiment (in seconds; during NC14; + 1 in case events last until the end of acquisition) from the nukes mat
    logi = first(:) < i;
    sints = [sints;i sum(logi)];
end
%normailization:
sints(:,2) = sints(:,2)./normMx;
%plot fit:
figure(figNum);
dataHandle = stairs(sints(:,1),sints(:,2),'LineWidth',1);shg

% outputs
mih = sints;