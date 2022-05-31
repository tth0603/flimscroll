
                surf4 = surfbf{:,4};
            %meta data stuff:
            zNSeries = size(DSPINMeta,1);
            zNPlanes = size(meat{1,1},1);
            %         Number of z stacks
            zNSlices=str2double(DSPINMeta.getPixelsSizeZ(0)); %in this data, the z stacks of one channel are all listed together first (ie index 1-4, for a set with 4 z stacks), and the other channel second (index 5-8)
            %         Number of channels
            %NChannels=surf4.getChannelCount(0);
            zImSeries = meat{1,1};
            zTestImage = zImSeries{zNPlanes,1}; 
            %position info:
             zPositionX = str2double(DSPINMeta.getPlanePositionX(0,0));      %TH may need to switch X and Y
            zPositionY = str2double(DSPINMeta.getPlanePositionY(0,0));   
             zx = str2double(DSPINMeta.getPixelsSizeX(0));
             zy = str2double(DSPINMeta.getPixelsSizeY(0));
            zPicSizeY = str2double(DSPINMeta.getPixelsPhysicalSizeY(0));
            %mx project the his channel:   
            hChan = 2;
            zImMat = [];
            %if hChan == 2    
                for i = 1:zNSlices %NSlices+1:NPlanes   
                zPlane = zImSeries{i,1};%+zNSlices,1};    %index for the last half of the images, which correspond to the second channel
                zImMat(:,:,i) = zPlane; %put all the planes in a 3D mat
                end
            %else %for his channel == channel 1
%                 for i = 1:NSlices %NSlices+1:NPlanes   
%                 pPlane = pImSeries{i,1};    %index for the last half of the images, which correspond to the second channel
%                 pImMat(:,:,i) = pPlane;
%                 end
%             end
            zMaxProj = max(zImMat,[],3);
            figure(71);imshow(imadjust(mat2gray(zMaxProj,[0 65535])),'DisplayRange',[],'InitialMagnification',100);shg
           
            ZoomRatio = zPicSizeY/PixelSizeX;
            zProjRe=imresize(zMaxProj, ZoomRatio);
            figure(72);imshow(imadjust(mat2gray(zProjRe,[0 65535])),'DisplayRange',[],'InitialMagnification',100);shg
                        
            figure(82);imshowpair(zProjRe,sMaxProj,'montage')
            C = normxcorr2(zProjRe, sMaxProj);
            %for troubleshooting, comment out when done:           
            figure(81); surf(C), shading flat
            %find the peak in the cross correlation:
            [ypeak, xpeak] = find(C==max(C(:)));
            %the 'normxcorr2' documentation explains this is to "Account for the padding that normxcorr2 adds":
            yoffSet = ypeak-size(zProjRe,1);
            xoffSet = xpeak-size(zProjRe,2);
            
            figure(90);
            hFig = figure;
            hAx  = axes;
        
            imshow(imadjust(mat2gray(sMaxProj,[0 65535])),'Parent', hAx);
            %a better way to plot a rectangle, 
            imrect(hAx, [xoffSet+1, yoffSet+1, size(zProjRe,2), size(zProjRe,1)]);
            
            
            
%%%%%%%%%%%%%%%%%%%%
%deleted from v5 of addParticle position
                p1 = zPositionX - coordA(1);
                q1 = zPositionY - coordA(2);
                %pos of acq image rel to ANT in pixels of surf image:
                p2 = p1/PixelSizeX;
                q2 = q1/PixelSizeX;
                %try somting wild:
                p3 = p2 + coordAx(1);
                q3 = coordAy(1) - q2;% + surfPicSizeY/2;
                
                %pos of acq image rel to ant surf im center in um:
                p1 = abs(sPositionX-zPositionX);
                q1 = abs(zPositionY-sPositionY);
                %pos of acq image rel to ant surf im center in pixels of surf image:
                p2 = p1/PixelSizeX;
                q2 = q1/PixelSizeX;
                %try somting wild:
                p3 = p2 + surfPicSizeX/2;
                q3 = q2 + surfPicSizeY/2;% + surfPicSizeY/2;

            hold on;plot(p3,q3,'r.','MarkerSize',10);
            hold on;plot(coordAx(1)+p2,coordAy(1)-q2,'y.','MarkerSize',10);shg
            figure(63);plot([zPositionX ],[zPositionY ],'r.','MarkerSize',10);shg
            hold on;plot([coordA(1)],[coordA(2)],'y.','MarkerSize',10);shg
            %hold on;plot(100,512,'b.','MarkerSize',10);shg
            figure(64);plot([zPositionX ],[zPositionY ],'r.','MarkerSize',10);shg
            hold on;plot([sPositionX ],[sPositionY ],'b.','MarkerSize',10);shg
            hold on;plot([tPositionX ],[tPositionY ],'c.','MarkerSize',10);shg
            pPositionX = str2double(pmeta4.getPlanePositionX(0,0));      %th these are in um, are are reversed for some reason
            pPositionY = str2double(pmeta4.getPlanePositionY(0,0));
            hold on;plot([-pPositionY ],[pPositionX],'c.','MarkerSize',10);shg
            
            hold on; figure(61);hold on
            plot(p3,q3,'y.','MarkerSize',10);shg
            hold on;plot(surfPicSizeX/2,surfPicSizeY/2,'c.','MarkerSize',10);shg
                %             %a series of transformations and scalings:
    %             surfx = sPositionX-(sx*PixelSizeX)/2;   %use LL corner (0,0) location of the surf image
    %             surfy = sPositionY-(sy*PixelSizeX)/2;
    %             %center of acq area realative to ANT position in um
    %             zPx2 = zPositionX - surfx;
    %             zPy2 = zPositionY - surfy;
    %             %center of acq area realative to ANT position in pix
    %             zPx3 = zPx2/(PixelSizeX);   %need to use the pixel size of the full embryo images, cuz in this fov, that's the conversion from pix to um
    %             zPy3 = zPy2/(PixelSizeX);
    %             %center of acq area realative to coord syst of fig 34 in pix
    %             zPx4 = zPx3;% 
    %             zPy4 = zPy3;% 
                %pos of acq image rel to ANT in um:
            
            %plot acq area and ANT pos
            hold on;plot(zPx4,zPy4,'c.','MarkerSize',10);
            hold on;plot(coordAx(1),coordAy(1),'r+','MarkerSize',10);
            acqArea =[zPx4-zImageSizeXpix*scl/2 zPy4-zImageSizeYpix*scl/2 zImageSizeXpix*scl zImageSizeYpix*scl]; %[LLxPos LL yPos width height] I'm sorry this illegible
            hold on;rectangle('Position',acqArea,'EdgeColor','c','LineStyle','-','LineWidth',5);shg
            end
            %geometry:      
        APAngle=atan2((coordP(2)-coordA(2)),(coordP(1)-coordA(1)));  %this assumes no funny distortions between the ANT image and the zoomed image, ie scaling in the x and y direction are the same. I thing this is safe to assume
        APlength=sqrt((coordP(1)-coordA(1))^2+(coordP(2)-coordA(2))^2); %units of um
        APlength2 = APlength/PixelSizeZoom; %units of pixels in the zoomed (ie acq'n) image
        %get coords of A & P in ref frame of zoomed image in units of pix
        aZM1 = [coordA(1)-zPositionX coordA(2)-zPositionY];%units of um
        pZM1 = [coordP(1)-zPositionX coordP(2)-zPositionY];
        aZM2 = aZM1/PixelSizeZoom; %units of zm image pix
        pZM2 = pZM1/PixelSizeZoom;
        
        
%%%%%%%%%%%%%%%%%%aligning mid sag images for findAP axis:


            %         Number of z stacks
            aNSlices=str2double(ameta4.getPixelsSizeZ(0)); %in this data, the z stacks of one channel are all listed together first (ie index 1-4, for a set with 4 z stacks), and the other channel second (index 5-8)
            %         Number of channels
            %NChannels=surf4.getChannelCount(0);
            aNPlanes = size(ameat{1,1},1);
            aImSeries = ameat{1,1};
            %position info:
            %mx project the his channel:   
            aImMat = [];
            %if hChan == 2    
                for i = aNSlices+1:aNPlanes 
                aPlane = aImSeries{i,1};%+zNSlices,1};    %index for the last half of the images, which correspond to the second channel
                aImMat(:,:,i) = aPlane; %put all the planes in a 3D mat
                end
            aMaxProj = max(aImMat,[],3);
            figure(71);imshow(imadjust(mat2gray(aMaxProj,[0 65535])),'DisplayRange',[],'InitialMagnification',100);shg
            
            pNSlices=str2double(pmeta4.getPixelsSizeZ(0)); %in this data, the z stacks of one channel are all listed together first (ie index 1-4, for a set with 4 z stacks), and the other channel second (index 5-8)
            %         Number of channels
            %NChannels=surf4.getChannelCount(0);
            pNPlanes = size(peta{1,1},1);
            pImSeries = peta{1,1};
            %position info:
            %mx project the his channel:   
            pImMat = [];
            %if hChan == 2    
                for i = pNSlices+1:pNPlanes 
                pPlane = pImSeries{i,1};%+zNSlices,1};    %index for the last half of the images, which correspond to the second channel
                pImMat(:,:,i) = pPlane; %put all the planes in a 3D mat
                end
            pMaxProj = max(pImMat,[],3);
            figure(72);imshow(imadjust(mat2gray(pMaxProj,[0 65535])),'DisplayRange',[],'InitialMagnification',100);shg
            figure(82);imshowpair(pMaxProj,aMaxProj,'montage')
            
            C = normxcorr2(pMaxProj, aMaxProj);
            %for troubleshooting, comment out when done:             
            figure(81); surf(C), shading flat
            %find the peak in the cross correlation:
            [ypeak, xpeak] = find(C==max(C(:)));
            %the 'normxcorr2' documentation explains this is to "Account for the padding that normxcorr2 adds":
            yoffSet = ypeak-size(pMaxProj,1);
            xoffSet = xpeak-size(pMaxProj,2);
            
            figure(90);
            %hFig = figure;     %may need this, 'imshow' below may plot to a weird fig, and this would solve that?
            hAx  = axes;
            imshow(imadjust(mat2gray(aMaxProj,[0 65535])),'Parent', hAx);
            %a better way to plot a rectangle, I don't understand why we
            %add '+1' to the x and y position, a fencepost correction? anyways the 
            %image position specifies the UL corner of the image:
            imrect(hAx, [xoffSet+1, 1, size(pMaxProj,2), size(pMaxProj,1)]); 
            %add A position to figure:
            hold on;plot(coordAx(1),coordAy(1),'r.','MarkerSize',10);
            
            
            %%%try something different:
            figure(50);
            imshowpair(aMaxProj, pMaxProj,'Scaling','joint')
            [optimizer, metric] = imregconfig('multimodal');
            movingRegistered = imregister(pMaxProj, aMaxProj, 'translation', optimizer, metric);
            %View the registered images.

            figure(51); 
            xLimits = [0 672*2];
            yLimits = [0 512*2];
            panoramaView = imref2d([512*2 672*2], xLimits, yLimits);
            imshowpair(aMaxProj,movingRegistered, 'OutputView', panoramaView);
            figure(89);
            imshow(warpedImage)
            
            figure(51);
            axh = axes([xLimits yLimits]);
            imshowpair(aMaxProj, movingRegistered,'Scaling','joint')
            tform = imregtform(pMaxProj, aMaxProj, 'translation', optimizer, metric);
            
            %put POST position in ref frame of ANT image:
            [tPostX,tPostY] = transformPointsForward(tform,postX,postY);
            
            figure(60);rectangle('Position',[0 0 lengthX lengthY]);
            hold on;plot(coordPpre(1), coordPpre(2),'r+','MarkerSize',10);
            hold on;plot(coordP(1), coordP(2),'c+','MarkerSize',10);
            hold on;plot(coordA(1), coordA(2),'b+','MarkerSize',10);
            
       