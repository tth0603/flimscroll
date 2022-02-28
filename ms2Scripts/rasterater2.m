function rasterater2(cia,nukeMat,minLength,figNum,varargin) 
% mib = rasterater(cia,minLength,figNum,unsort)
%
% This is for use with analData.nukesMat straight outta flimscroll. like
% the fn rasterater, but using cia from makeInts. should provide more
% flexibility
%
% unfinished
%
% INPUTS:
%   cia == a nuke mat from the analData struc from flimscroll.
%       format: '1.nucleus 2.low or high (=-2,0,2 or -3,1,3) 3.frameStart 4.frameEnd 5.deltaFrames 6.deltaTime(sec)';
%   nukeMat == a nuke mat from the analData struc from flimscroll.
%       format: %[1.time(s) 2.frameNum 3.nucleus 4.binarySpot 5.nukeXpos 6.nukeYpos 7.spotGaussAmp 8.spotSigma 9.spotOffset 10.spotIntegratedIntensity 11.relativeFrameNumber 12.NC]
%   minLength == a scaler indicating the maximal length a trace must be to
%   be included
%   figNum == where you gunna put it? (a scaler)

%
% OUTPUTS
%   a raster plot of events (with heat map).
%
% USAGE
%   rasterater2(analData.cia,1,30)
%   rasterater2(cia,1,30)
%
% Harden 2018

%get the min/max frame. we get this from nukeMat. 
mx = max(nukeMat(:,2));
mn = min(nukeMat(:,2));
%for each nuke, make a row:
nukes = unique(cia(:,1));
rasterMat = [];
for i = 1:length(nukes)
    %get all ints for this nuke
    mat = cia(cia(:,1) == nukes(i),:);
    noe = size(mat,1);
    for j = 1:noe
    end
    
end

%plot
figure(figNum+1);
imagesc(rasterMatSort,[0 300]);shg
j = parula; 
%make bkgrnd white
j(1,:) = [1 1 1];
colormap(j)
colorbar;

end

%notes for posterity
% % get the NaN's outta here 
% mat(isnan(mat)) = 0;
% % get those nukes:
% nukes = unique(mat(:,3));
% rasterMat = [];
% for i = 1:length(nukes)
%     mat2 = mat(mat(:,3) == nukes(i),:); %get all rose and cols for a given nuke
%     if sum(mat2(:,4)) > 0 %any spots here? (of proper length)
%         rasterMat = [rasterMat;mat2(:,7)']; %cols become rose
%     end
% end
% if minLength > 0
%     logiMat = rasterMat > 0 ;
% %sort rasterMat (there are a million ways to do this. this shitty soln prolly isnt the best):
% indV = [];
% for j = 1:size(rasterMat,1)
%     [~,col] = (find(rasterMat(j,:),1)); %find first col with a non zero value in each row
%     indV = [indV;col]; %place these column indices in a vector the length of rasterMat
% end
% %append these col indices to the first column of rasterMat
% rasterMatInd = [indV rasterMat];
% % sort this based on teh first col:
% sortMat = sortrows(rasterMatInd,1);
% %get rid of the first col (the indicies)
% rasterMatSort = sortMat(:,2:end);
% %if you want unsorted output:
% if nargin == 3 && ischar(varargin{1})
%     rasterMatSort = rasterMat;
% end