function fullTraceArrayFli(analData,figNum)
% This script will plot the integrated spot intenstiy traces from all
% nuclei with a spot from the nukesMat member of the analData strucure
% flimscroll2 outputs using matlab subplots. Similar to
% specifiedTraceArrayFli.m 
% Now includes integrated intensity for times between spots. requires an
% updated analData. See t1p54
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

%nukesMat:
nukesMat = analData.nukesMat; %[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.relativeFrameNumber 12.NC 13.nuclearAPpos]
noSpotMat = analData.noSpots; %[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.noSpotXpos 6.noSpotYpos 7.gaussAmp 8.sigma 9.offset 10.intInt 11.relativeFrameNumber]'

%get frame and Nuke info
frameV = unique(noSpotMat(:,2)); 
nukeV = unique(noSpotMat(:,3)); %select nukes from noSpotMat
frameNum = size(frameV,1);
nukeNum = size(nukeV,1);
%cycle through all frames, then all nukes to replace NaNs in nuke mat with
%noSpot info:
for ii = 1:frameNum
    frame = frameV(ii);
    for j = 1:nukeNum
        nuke = nukeV(j);
        %get the nukesMat indicie
        nukeLogi = nukesMat(:,2) == frame & nukesMat(:,3) == nuke;
        % if there isnt a spot, add no spot info
        if nukesMat(nukeLogi,4) == 0
            %get the no spot data:
            noSpotLogi = noSpotMat(:,2) == frame & noSpotMat(:,3) == nuke;
            noSpotRow = noSpotMat(noSpotLogi,:);
            %assign the no spot values to the nukesMat at the correct row:
            nukesMat(nukeLogi,7:10) = noSpotRow(7:10);
        end
    end
end

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

%set plot colors
noSignalC = [0.3 0.3 0.3];
signalC = [0.23 0.44 0.33];

for i = 0:figTotal - 1
    figure(figNum+i);
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
        % try to determine the number of lines needed for each trace
        %find the number of events:
        binaryV = mat(:,4);
        d = [true, diff(binaryV') ~= 0, true];  % TRUE if values change
        n = diff(find(d));  %length of this is number of events
        eventV = cumsum(n);
        eventNum = length(eventV);
        for k = 1:eventNum
            maxInd = eventV(k);
            try
                minInd = eventV(k - 1) + 1;
                flag = 1;
            catch
                minInd = 1;
                flag = 0;
            end
            h = plot(mat(minInd:maxInd,1),mat(minInd:maxInd,10),'LineWidth',1);
            if mat(minInd,4) == 0
                set(h,'Color',noSignalC);
            else
                set(h,'Color',signalC);
            end
            hold on
            % fill in the gap:
            if flag == 1
                h2 = plot(mat(minInd - 1:minInd,1),mat(minInd - 1:minInd,10),'LineWidth',1); % this makes the time step curve just prior to the event green
                set(h2,'Color',signalC);
            end
        end 
        %plot the entire curve one color:
%         %grab the time values
%         ts = mat(:,1);
%         %grab intensity values
%         ys = mat(:,10);
%         plot(ts,ys,'Color',[0.23 0.44 0.33]);shg
        % plot the signal/no signal different colors:
%         signal = mat(mat(:,4) == 1,:);
%         noSignal = mat(mat(:,4) == 0,:);
%         plot(signal(:,1),signal(:,10),'.','Color',[0.23 0.44 0.33]);
%         hold on;plot(noSignal(:,1),noSignal(:,10),'.','Color',[0.3 0.3 0.3]);shg
        title(sprintf('%d',nuke));
        if j == axesLabel
            xlabel('Time (s)','FontSize',16)    %Puts one x axis label per figure
            ylabel('integrated intensity (au)','FontSize',16)
        end
    hold on    
    end
end

