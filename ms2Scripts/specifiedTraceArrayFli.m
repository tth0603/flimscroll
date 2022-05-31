function specifiedTraceArrayFli(nukesMat,nukeVect,figNum)
% This script will plot many integrated intenstiy traces in several plots
% using subplots from the nukesMat member of the analData strucure
% flimscroll2 outputs. You must specify which nukes you'd like to plot.
%
% INPUTS:
%   nukesMat == matrix member of the same name in analData struc.
%       EX: nukesMat = analData.nukesMat;
%   nukeVect == a vector of integers of the nuclei numbers of which you'd
%       like to plot.
%       EX: v = [100:108];
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
% number of traces
L = length(nukeVect);
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
            nuke = nukeVect(i*arrayDim^2 + j);
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
