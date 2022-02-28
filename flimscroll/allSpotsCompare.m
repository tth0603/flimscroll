function mib = allSpotsCompare(AllSpots1,AllSpots2)
%   mib == allSpotsCompare(AllSpots1,AllSpots2)
%
%   For use within flimscroll2. Will take two allSpots structures and
%   return a single structure containg info for spots that reside in both
%   input structure. For two condition spot selection.
%
%   Timothy Harden 2019

cell1 = AllSpots1.AllSpotsCells; %this is where all interesting info lives
cell2 = AllSpots2.AllSpotsCells;
%the spot comparing:
len = size(cell1,1); %the number of frames we picked spots for
for i = 1:len
    locs1 = cell1{i,1}; %lists of locations
    locs2 = cell2{i,1};
    len2 = size(cell1{i,1},1); %number of spots in frame i of spot set 1
    logiV = false(1,len2)'; %allot space
    for j = 1:len2
        tstAoi = locs1(j,:); %test spot to compare against all in locs2
        dtest = sqrt( (locs2(:,1) - tstAoi(1,1)).^2 + (locs2(:,2) - tstAoi(1,2)).^2); %distances btwn test spot and all others
        logi = any(dtest < 1.5); %check if test spot is within 1.5 pixels of any spot in cell2 (if so keep that spot)
        logiV(j) = logi; %list of spots to keep
    end
    AllSpotsCells{i,1} = cell1{i,1}(logiV,:); %assign the matching spot locations to a new cell array
    AllSpotsCells{i,2} = sum(logiV);  %number of spots
    AllSpotsCells{i,3} = cell1{i,3};  %frame number
    AllSpotsCells{i,4} = cell1{i,4}(logiV,:); %[amp sigma offset integratedInt] list for spots in frame i
end
%make the struc, add the spots cell array:
AllSpots.AllSpotsCells = AllSpotsCells;
%save all bits of AllSpots that do not change btwn the two cndns:
AllSpots.AllSpotsCellsDescription = AllSpots1.AllSpotsCellsDescription;
AllSpots.FrameVector = AllSpots1.FrameVector;
AllSpots.Parameters = [AllSpots1.Parameters AllSpots2.Parameters]; %include the parameters from each condition
AllSpots.ParametersDescripton = AllSpots1.ParametersDescripton;
AllSpots.aoiinfo2 = AllSpots1.aoiinfo2; %this doesn't really matter
AllSpots.aoiinfo2Description = AllSpots1.aoiinfo2Description;

mib = AllSpots;