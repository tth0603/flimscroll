function out = matchMakerV0(aoifits, allspots)
% work on MS2/his analysis;



%working with MS2 spot, expl from flimscroll
frms=handles.AllSpots.FrameVector;      % List of frames over which we found spots
hold on
[rose col]=size(frms);

for indx=1:max(rose,col)                           % Cycle through frame range
    
    spotnum=handles.AllSpots.AllSpotsCells{indx,2}; % number of spots found in the current frame
     xy=handles.AllSpots.AllSpotsCells{indx,1}(1:spotnum,:);    % xy pairs of spots in current frame
     plot(xy(:,1),xy(:,2),'y.');                % Plot the spots for current frame
end
hold off