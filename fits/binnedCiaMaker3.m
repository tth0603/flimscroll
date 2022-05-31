function mih = binnedCiaMaker3(varargin)
% mih = binnedCiaMaker(analData1,analData2,analData3,...,bins,binsOfInterest)
% 
% will make a cia from any number of analData files with only events within
% binsOfInterest. as for binnedCiaMaker & binnedCiaMaker2, but takes any
% number of input analData structures.
%
% INPUTS
%   bins == a vector covering the AP length of interest: eg bins = [-0.20:0.02:0.20];
%   binsOfInterest == the indicies of the bins that you want to select events from: eg boi = 5:8;
%
% OUTPUS
%   MIH = a cia
%
% Timothy Harden 20210818

%initialize some variables:
numberOfReps = length(varargin) - 2;
for i = 1:numberOfReps
    tempStruc = varargin{i};
    ciaCell{i} = tempStruc.cia;
end
bins = varargin{numberOfReps + 1};
binsOfInterest = varargin{numberOfReps + 2};

% place the nukes into bins:
for i = 1:numberOfReps
    out{i} = positionDistribution(varargin{i},bins);
end
% total number of nukes in each bin:
nBinTotals = zeros(1,size(out{1}.totalBin,2))';
for i = 1:numberOfReps
    nBinTotals = nBinTotals + out{i}.totalBin;
end

%for each rep, pull out the nukes that we want:
for i = 1:numberOfReps
    eventPosCell{i} = out{i}.eventPos; %[1.nukeNumber 2.nukePostion 3.binAssignment]
end

comboCia = [];
for i = 1:numberOfReps % one for each analData
    ciaMat = ciaCell{i};
    eventPosMat = eventPosCell{i};
    for j = 1:length(binsOfInterest)
        % get the nukes that belong in bin j
        logi1 = eventPosMat(:,3) == binsOfInterest(j);
        % if there are any nukes:
        if sum(logi1) ~= 0
            nukeNums = eventPosMat(logi1,1);
        % add these to a new cia, 
            for k = 1:length(nukeNums)
                logi2 = ciaMat(:,1) == nukeNums(k);
                addOn = ciaMat(logi2,:);
                if i > 1
                    addOn(:,1) = addOn(:,1) + 1000*(i - 1); % tack on a number so we dont confuse nukes from two different experiments
                end
                comboCia = [comboCia;addOn];
            end
        end
    end
end
mih = comboCia;