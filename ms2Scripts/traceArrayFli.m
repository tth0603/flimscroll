function traceArrayFli(nukesMat,figNum)
% This script will plot the integrated spot intenstiy traces from all
% nuclei with a spot from the nukesMat member of the analData strucure
% flimscroll2 outputs using matlab subplots. Similar to
% specifiedTraceArrayFli.m 
%
% This script generates a lot of plots. To close the use: CLOSE ALL at the
% command line
%
% INPUTS:
%   nukesMat == matrix member of the same name in analData struc.
%       EX: nukesMat = analData.nukesMat;
%   figNum == A scaler. The first figure to be plotted in. Starting from
%       figNum the script will use as many figures as needed to display
%       traces from all nukes in nukesVect.
%
% Timothy Harden 20200202

%as an imperfect measure of baseline intensity, us the dimmest spot
%measured:
minInt = min(nukesMat(:,10));
%get rid of all the NaNs and replace with our min intensity
nukesMat(isnan(nukesMat)) = minInt;
%get the nukes with at a spot in at least one frame
nukeV = unique(nukesMat(:,3)); %all nuke numbers
NOI = [];
newNukesMat = [];
for i = nukeV'
    logi = nukesMat(:,3) == i;
    mat = nukesMat(logi,:);
    if sum(mat(:,4)) > 0
        NOI = [NOI; i];
        newNukesMat = [newNukesMat; mat];  %if the current nuke shows even one spot, we add that data to a new mat
    end
end
nukesMat = newNukesMat;
% number of traces
L = length(NOI);
% set up plotting stuff:
%%%%%%%%%% change this for 2x2, 3x3 arrays, etc...
arrayDim = 2; 
%%%%%%%%%%
figTotal = ceil(L/arrayDim^2);
%set up where to put axes labels
if arrayDim == 2
    axesLabel = 3;
elseif arrayDim == 3
    axesLabel = 7;
else
    axesLabel = arrayDim^2;
end

for i = 0:figTotal - 1
    hold on;figure(figNum+i);
    for j = 1:arrayDim^2
        %set up subplot
        subplot(arrayDim,arrayDim,j);
        %complicated way to pick out each nuke:
        try
            nuke = NOI(i*arrayDim^2 + j);
        catch
            return
        end
        %pick out the data for each nuke:
        logi = nukesMat(:,3) == nuke;
        mat = nukesMat(logi,:);
        %grab the time values
        ts = mat(:,1);
        %grab intensity values
        ys = mat(:,10);
        plot(ts,ys,'Color',[0.23 0.44 0.33]);shg
        title(sprintf('%d',nuke));
        if j == axesLabel
            xlabel('Time (s)','FontSize',16)    %Puts one x axis label per figure
            ylabel('integrated intensity (au)','FontSize',16)
        end   
    end
end

