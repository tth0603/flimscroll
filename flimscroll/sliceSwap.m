function sliceSwap(framesToBeSwapped,varargin)
% sliceSwap(framesToBeSwapped)
%
% This script addresses when the NIC diSPIM acts up and switches channels midway
% through an acquisition. Use maxProjLS to create the maximum projection
% tiffs. Then feed this script the frames you would like to swap between
% two tiff files. As of now these files must have the same number of
% frames. The script will ask the user to select the files by hand.
%
% INPUTS:
%   framesToBeSwapped == a vector containing the frames to be swapped. Fiji
%   is useful to determine this range. Ex: [1:5,115:121]
%
% OUPUTS:
%   It will write two tiff files containing the corrected time series to
%   the same directory the original tiff files are contained in. The files
%   names will be appended with '2' at the end. 
%
% USAGE:
%   sliceSwap([1:5,115:121]);
%
% Timothy Harden 2019

tic;
f = framesToBeSwapped;
%UI get the files:
[fn1 fp1] = uigetfile('*.tiff','pick a max proj tiff');
currDir = cd;
cd(fp1);
[fn2 fp2] = uigetfile('*.tiff','pick the other max proj tiff');
cd(currDir);
%load the meta data
r1 = bfGetReader([fp1 fn1]); %low level retireival fn to get file reader without loading the data
meta1 = r1.getMetadataStore(); 
r2 = bfGetReader([fp2 fn2]);
meta2 = r2.getMetadataStore(); 
%find the number of total frames and such:
nT = meta1.getPixelsSizeZ(0).getValue(); %bc the way I save max projs with bf save, the time frames are interpreted as Z values
nX = meta1.getPixelsSizeX(0).getValue();
nY = meta1.getPixelsSizeY(0).getValue();
%first get the original (fucked up) image stack from file 1:
newStack1 = ones(nY,nX,nT);
for iPlane = 1:nT
    newStack1(:,:,iPlane) = bfGetPlane(r1, iPlane);
end
%get the misplaced frames in file 2 to place into file 1 im stack:
for iPlane = f
    %correct frames from 2:
    newStack1(:,:,iPlane) = bfGetPlane(r2, iPlane);
end
%do the same for file 2:
newStack2 = ones(nY,nX,nT);
for iPlane = 1:nT
    newStack2(:,:,iPlane) = bfGetPlane(r2, iPlane);
end
%get the misplaced frames in file 1 to place into file 2 im stack:
for iPlane = f
    %correct frames from 1:
    newStack2(:,:,iPlane) = bfGetPlane(r1, iPlane);
end
%save the files (with fn appended):
%file 1:
fn1prefix = strtok(fn1,'.');
bfsave(newStack1,[fp1 fn1prefix '2.tiff'], 'BigTiff', true);
%file 2:
fn2prefix = strtok(fn2,'.');
bfsave(newStack2,[fp2 fn2prefix '2.tiff'], 'BigTiff', true);
toc