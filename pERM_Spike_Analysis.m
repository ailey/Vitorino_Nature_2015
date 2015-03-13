function pERM_Spike_Analysis(varargin)
%--------------------------------------------------------------------------
% pERM_Spike_Analysis.m
% This function reads in an RGB 20x image (DAPI in blue, actin in
% red, & pERM in green)
% and computes morphological parameters of the whole cell and spikes
% specifically.
%
%   Author: Ailey Crow
%   Date created: March 2012
%   Most recent revision: March 13, 2013
%
%   Inputs:     - datadir: directory where raw data, results folder, etc. are
%                  located
%               - samplename: name of sample to be analyzed (as entered in
%                  datadir)
%
%   Outputs:    None
%
%   Exports:    - _raw.tif: original image
%               - _mask.tif: cell and fibrillar adhesion mask overlayed on
%                   adhesion channel
%               - _periph.tif: peripheral region mask & actin or 2ary
%                   adhesion mask overlayed on 2ary adhesion channel
%               - Adhesion_Detection.xls: quantification results
%
% Notes:
% WC = whole cell, NP = No Protrusions
%--------------------------------------------------------------------------

%checks for formatting of input arguments
 if(nargin<1)
 	fprintf('Usage:\n');
 	fprintf('Adhesion_Detection -datadir DIRECTORY -samplename SAMPLENAME\n');
 	fprintf('  where DIRECTORY is directory containing image data\n');
        fprintf('  where SAMPLENAME is the name of the image to be analyzed\n');
 
 	fprintf('Example:\n');
 	fprintf('Adhesion_Detection -datadir datadir -samplename char(List2Analyze(1))\n');
 	return;
 end
 
 datadir = '';
 
 %assigns working directory datadir & samplename according to VARARGIN
 for k=1:nargin
 	switch varargin{k}
 		case '-datadir'
 			if(nargin>=(k+1)) %makes sure there is an input after "-datadir"
 				datadir = varargin{k+1}; %assigns input after "-datadir" to datadir
 			end
 		case '-samplename'
            if(nargin>=(k+1)) %makes sure there is an input after "-samplename"
                samplename = varargin{k+1}; %assigns input after "-samplename" to samplename
            end
        otherwise
 	end
 end


%% Setup
%Select relevant directories
pdir = '~/Adhesion_Detection/';        % directory where program is saved
tdir = '~/Toolbox/';                           % directory of useful generic m-files
resdir = strcat(datadir, 'Results');        % directory where output & parameters.mat are saved. 

% load parameters
cd(resdir)
load('parameters.mat')


%% Load & assemble image
% This section reads in the image of interest using subroutines specific to the 
% imaging systems used for this project. The image is saved as a 3D matrix of 
% intensity values called im. The scale (pixel-to-micron ratio) is also loaded.
% To run the code, replace these subroutines with code to load images for your 
% specific system. (Contact the vendor of your imaging system if you need help 
% creating or acquiring such code.)
cd(tdir)
switch acquisition
    case 'Ariol'
        [im, scale] = Load_ARIOL_im(datadir, samplename);
        savename = samplename;
    case 'TissueGnostics'
        [im, scale, datadir, samplename] = Load_TG_im(datadir, samplename, 0);
        savename = samplename;
    case 'NanoZoomer'
        lowmagpower = 10;
        [im, scale, Xpix, Ypix] = Load_gSlide_im(samplename, lowmagpower);
        tmp = regexp(samplename, '/', 'split');
        tmp = regexp(tmp{numel(tmp)}, ' ', 'split');
        savename = tmp{1};
    case 'ImageXpress'
        cd(tdir)
        [im, scale] = Load_ImageXpress_im(datadir, samplename, mag);
        savename = samplename;
    otherwise
end

imRed = uint8(im(:,:,actin_ch));
imGreen = uint8(im(:,:,Ad_ch));
imBlue = uint8(im(:,:,DAPI_ch));

RGBim = cat(3, imRed, imGreen, imBlue);


%% Significant parameters
Disp = 0;                        % Disp = 1 displays a number of figures. Set to 0 for batch analysis.
Pstel = strel('disk', ceil(15*scale));        % Structuring element for distinguishing cell body from protrusions during morphologial opening. Default disk, 4
minP_size = ceil(100*scale);                   % minimum protrusion size (pixels). Default: 7
imG_thresh = 40;                 % Threshold for counting pERM+ staining in protrusion. Default: 40.

close all


%% Create whole-cell mask
imB_mask = imBlue > 1.2* imfilter(imBlue, fspecial('disk', ceil(500*scale)), 'symmetric', 'conv');    % Segment DAPI image to create mask of nuclei
imB_mask = imB_mask & ~bwareaopen(imB_mask, ceil(scale*1e5), 4);
imB_mask = bwareaopen(imB_mask, ceil(1000*scale), 4);                        % remove small junk in DAPI channel
imB_mask = imfill(imB_mask, 'holes');

Rmask1 = imRed > 0.7*mean(imRed(imB_mask));
Rmask2 = imRed > 1.1*imfilter(imRed, fspecial('disk', 50), 'symmetric', 'conv');
Rmask = imclose(Rmask1 | Rmask2, strel('disk', 2));

Gmask1 = imGreen > 0.7*mean(imGreen(imB_mask));
Gmask2 = imGreen > 1.1*imfilter(imGreen, fspecial('disk', 50), 'symmetric', 'conv');
Gmask3 = edge(imGreen, 'Canny', [0.05,0.1]);
Gmask = Gmask1 | Gmask2 | Gmask3;

im_WCmask = Gmask | Rmask;

im_WCmask = ~bwareaopen(~im_WCmask, ceil(1000*scale));                     % fill in holes (low-intensity areas) within cell body
im_WCmask = bwareaopen(im_WCmask, ceil(10*scale), 8);
im_WCmask = bwareaopen(im_WCmask, 1000, 8);                 % Remove all non-cell (smaller) objects


%% Create mask of cell body without protrusions
im_NPmask = imopen(im_WCmask, Pstel);                       % Morphologically open image to remove protrusions
im_NPmask = bwareaopen(im_NPmask, 800, 4);                 % Remove smaller objects


%% Isolate single cells & remove balled-up (unhealthy) cells
[Cells, numCells] = bwlabel(im_WCmask);                         % Find number of cells of interest

for c = 1:numCells
    tempR_mask = ismember(Cells, c);                         % Select cell mass (connected area). Note: bwlabel creates a mask where each connected area is given a distinct integer value.
    tempB_mask = imB_mask & tempR_mask;                    % Select nuclei in cell mass
    [~, numNuclei] = bwlabel(tempB_mask);                   % Count nuclei in cell mass
    if (numNuclei<1)                              % If cell mass has no nucleus
        im_NPmask(tempR_mask) = 0;                          % delete from no-protrusion image mask
        im_WCmask (tempR_mask) = 0;                         % delete from whole-cell image mask
    end
end

im_WCmask = bwareaopen(im_WCmask, 1000, 4);                  % Remove all non-cell (smaller) objects

% Highlight the edge of the NP mask for display purposes
if(Disp==1)
    figure; imshow(imoverlay(imGreen, bwperim(im_NPmask) | bwperim(im_WCmask), [0.3 1 0.3]))
end


%% Create mask of protrusions only
prot_mask = im_WCmask & ~im_NPmask;                 % Subtract cell body (No protrusion) mask from whole cell mask
prot_mask = bwareaopen(prot_mask, minP_size);      % Remove tiny "protrusions"
prot_maskBIG = bwareaopen(prot_mask, ceil(2000*scale));        % Identify huge "protrusions"
prot_mask = prot_mask & ~prot_maskBIG;              % Remove huge protrusions

% % check intensity of protrusion
PL = bwlabel(prot_mask);
Pprop = regionprops(logical(prot_mask), imGreen, 'MeanIntensity', 'Solidity');
Int_index = find([Pprop.MeanIntensity] > 0.75*mean(imGreen(imB_mask)));                % Select Cells objects with nuclei
Pprot_mask1 = ismember(PL, Int_index);
Int_index = find([Pprop.Solidity] > 0.65);                % Select Cells objects with nuclei
Pprot_mask2 = ismember(PL, Int_index);
Pprot_mask = Pprot_mask1 & Pprot_mask2;

% Calculate protrusion properties before isolating single cells
[~, numNuc] = bwlabel(imB_mask);                     % Find number of nuclei
[~, numProt] = bwlabel(Pprot_mask);                     % Find number of nuclei
perim = sum(sum(bwperim(im_NPmask)));
ProtperNuc = numProt / numNuc;
ProtperPerim = numProt / perim;

% Save image for all cells
Poutline = imdilate(bwperim(Pprot_mask), strel('disk', 2));
NPoutline = bwperim(im_NPmask);
outfilename = sprintf('%s_mask_AllCells.tif', savename);
temp_img = cat(3, uint8(imRed), uint8(imGreen), uint8(Poutline | NPoutline).*220);
Save_Image(temp_img, resdir, outfilename, 5e8);


%% Second approach to identifying protrusions
Prot_Area = imdilate(im_NPmask, strel('disk', ceil(30*scale))) & ~im_NPmask;   % Define area just outside cell as potential protrusion territory
prot_mask_alt = (imGreen > 0.9*mean(imGreen(im_NPmask))) & Prot_Area;
prot_mask_alt = bwareaopen(prot_mask_alt, ceil(minP_size/3));      % Remove tiny "protrusions"
prot_maskBIG = bwareaopen(prot_mask_alt, ceil(2000*scale));        % Identify huge "protrusions"
prot_mask_alt = prot_mask_alt & ~prot_maskBIG;              % Remove huge protrusions

% Calculate protrusion properties before isolating single cells
[~, numProt] = bwlabel(prot_mask_alt);                     % Find number of nuclei
ProtperNuc_alt = numProt / numNuc;
ProtperPerim_alt = numProt / perim;

% Save image for all cells
Poutline = imdilate(bwperim(prot_mask_alt), strel('disk', 2));
NPoutline = bwperim(im_NPmask);
outfilename = sprintf('%s_mask_alt_AllCells.tif', savename);
temp_img = cat(3, uint8(imRed), uint8(imGreen), uint8(Poutline | NPoutline).*220);
Save_Image(temp_img, resdir, outfilename, 5e8);


%% Isolate single cells and clear border for rest of code
[Cells, numCells] = bwlabel(im_WCmask);                         % Find number of cells of interest

for c = 1:numCells
    tempR_mask = ismember(Cells, c);                         % Select cell mass (connected area). Note: bwlabel creates a mask where each connected area is given a distinct integer value.
    tempB_mask = imB_mask & tempR_mask;                    % Select nuclei in cell mass
    [~, numNuclei] = bwlabel(tempB_mask);                   % Count nuclei in cell mass
    if (~isequal(numNuclei,1))                              % If cell
        im_NPmask(tempR_mask) = 0;                          % delete from no-protrusion image mask
        im_WCmask (tempR_mask) = 0;                         % delete from whole-cell image mask
    end
end

im_WCmask = bwareaopen(im_WCmask, 1000, 4);                  % Remove all non-cell (smaller) objects
im_WCmask = imclearborder(im_WCmask);
prot_mask = prot_mask & im_WCmask;


%% Determine cell & protrusion properties for each cell of interest & record
% Display cell objects and protrusion objects
[~, numNuc] = bwlabel(imB_mask);                     % Find number of nuclei
[CellsWC, numCells] = bwlabel(im_WCmask);
Cell_WC = regionprops(im_WCmask, 'BoundingBox', 'Perimeter', 'Area');
if(Disp==1);
    CellDisp = CellsWC;
    CellDisp(prot_mask) = 10;
    figure; imagesc(CellDisp);
end

for c = 1:numCells
    
    % Cell stats
    tempWC_mask = imcrop(ismember(CellsWC, c), Cell_WC(c,1).BoundingBox);
    tempNP_mask = imcrop(im_NPmask, Cell_WC(c,1).BoundingBox);                        % Select single-cell no protrusion mask
    tempNP_mask = tempNP_mask & tempWC_mask;
    if(max(max(tempNP_mask)) == 0), break, end
    Cell_stats = regionprops(tempNP_mask, 'Area', 'Orientation', 'Perimeter', 'Centroid', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength');    % Determine cell properties
    O = round(Cell_stats(1).Orientation);
    Cell_AR = Cell_stats(1).MajorAxisLength / Cell_stats(1).MinorAxisLength;
    
    % Nucleus stats
    tempB_mask = imcrop(imB_mask, Cell_WC(c,1).BoundingBox);
    tempB_mask = tempNP_mask & tempB_mask;                           % Select cell's nucleus
    [~, numNuc] = bwlabel(tempB_mask);
    Nuc_stats = regionprops(tempB_mask, 'Area', 'Perimeter', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength');    % Determine nucleus properties
    Nuc_AR = max(Nuc_stats.MajorAxisLength) / max(Nuc_stats.MinorAxisLength);
    
    % Protrusion stats
    tempWC_mask = imcrop(im_WCmask, Cell_WC(c,1).BoundingBox);
    tempWC_mask = bwareaopen(tempWC_mask, 1000, 4);         % Remove all objects not connected to cell of interest
    tempP_mask = tempWC_mask & ~tempNP_mask;
    tempP_mask = bwareaopen(tempP_mask, minP_size);      % Remove tiny "protrusions"
    
    [ProtL, numProts] = bwlabel(tempP_mask);
    Prots = regionprops(ProtL, 'Area', 'MajorAxisLength');  % Determine parameters for all protrusions
    Length_ave = mean([Prots(:).MajorAxisLength]);
    Length_std = std([Prots(:).MajorAxisLength]);
    Area_ave = mean([Prots(:).Area]);
    Area_std = std([Prots(:).Area]);
    
    % Look for a cluster of protrusions & determine location of cluster
    % with respect to long & shortaxes of cell body
    C = nan;
    if(numProts >= 3)                                           % if less than 3 protrusions, don't look for a cluster
        tempNP_bound = bwboundaries(tempNP_mask, 8, 'noholes'); % Find coordinates of no protrusion boundary
        tempP_bound = bwboundaries(tempWC_mask, 8, 'noholes');  % Find coordinates of whole cell boundary
        tNP = logical(tempNP_mask);
        CBox = regionprops(tNP, 'Centroid');
        IN1 = inpolygon(CBox(1).Centroid(2),  CBox(1).Centroid(1), tempNP_bound{1}(:,1), tempNP_bound{1}(:,2));
        IN2 = inpolygon(CBox(1).Centroid(2),  CBox(1).Centroid(1), tempP_bound{1}(:,1), tempP_bound{1}(:,2));
        if IN1 && IN2                                                  % if centroid is not within cell boundaries (unusual) skip this cell
            [NPdist, NPangle] = signature(tempNP_bound{1}, CBox(1).Centroid(2), CBox(1).Centroid(1));    % Determine amplitude vs angle along perimenter wrt centroid
            [WCdist, WCangle] = signature(tempP_bound{1}, CBox(1).Centroid(2), CBox(1).Centroid(1));     % Determine amplitude vs angle along perimenter wrt centroid
            
            % Unfortunately signature does not provide an amplitude for
            % every angle and to subtract NPdist from WCdist we need
            % the matrices to be equivalent so we now embark on
            % interpolating to fill in missing angles
            for a = 0:359
                b=a+1;                                             % b marks row position in matrix
                if ~any(NPangle == a)                      % if angle is not in NPangle list
                    if(a==0);
                        x1 = a-1; x2 = NPangle(b);
                        y1 = NPdist(numel(NPdist)); y2 = NPdist(b);
                    elseif (a==numel(NPangle))
                        x1 = a-1; x2 = 360;
                        y1 = NPdist(a); y2 = NPdist(1);
                    else
                        x1 = a-1; x2 = NPangle(b);
                        y1 = NPdist(a); y2 = NPdist(b);
                    end
                    [tmpAngle, tmpDist] = interpPts(x1, x2, y1, y2);    % interpolate for every integer value between x1 & x2
                    NPangle = [NPangle(1:a); tmpAngle; NPangle(b:numel(NPangle))];  % Add interpolated points
                    NPdist = [NPdist(1:a); tmpDist; NPdist(b:numel(NPdist))];
                end
                
                if ~any(WCangle == a)                      % if angle is not in WCangle list
                    if(a==0);
                        x1 = a-1; x2 = WCangle(b);
                        y1 = WCdist(numel(WCdist)); y2 = WCdist(b);
                    elseif (a==numel(WCangle))
                        x1 = a-1; x2 = 360;
                        y1 = WCdist(a); y2 = WCdist(1);
                    else
                        x1 = a-1; x2 = WCangle(b);
                        y1 = WCdist(a); y2 = WCdist(b);
                    end
                    [tmpAngle, tmpDist] = interpPts(x1, x2, y1, y2);    % interpolate for every integer value between x1 & x2
                    WCangle = [WCangle(1:a); tmpAngle; WCangle(b:numel(WCangle))]; % Add interpolated points
                    WCdist = [WCdist(1:a); tmpDist; WCdist(b:numel(WCdist))];
                end
            end
            
            % only include angles 0-359 since 0 & 360 are equivalent
            if any(NPangle == 360)
                NPangle = NPangle(1:360);
                NPdist = NPdist(1:360);
            end
            if any(WCangle == 360)
                WCangle = WCangle(1:360);
                WCdist = WCdist(1:360);
            end
            
            %Pangle = [0:360]';
            %                 if(Disp==1), figure; plot(NPangle, NPdist, 'r'); hold on; plot(WCangle, WCdist, 'b'); title('Signatures w/ & w/out protrusions'); end
            Pdist = WCdist - NPdist;
            %                 if(Disp==1), figure; plot(NPangle, Pdist); title('Protrusion signature'); end
            
            Pdist = [Pdist; Pdist; Pdist];                          % Create periodic boundaries
            [~, locs] = findpeaks(Pdist, 'minpeakdistance', 3, 'minpeakheight', 2);   % identify protrusions (peaks)
            P = zeros(size(Pdist));
            P(locs) = 1;
            Pl = smooth(P, 45, 'loess');                % Smooth using a 45 degree window & local regression weighted linear least squares method to find clusters of peaks
            Clust = Pl( O+360 : O+720);
            %                 if(Disp==1), figure; plot(Clust), end
            if ~any(Clust > 0.17)                        % Threshold for cluster is smoothed height > 0.17
                C = nan;
            else
                C = find(Clust == max(Clust));
            end
        end
    end

    % find pERM+ protrusions and subtract pERM- protrusions from protrusion mask
    tempG = imcrop(imGreen, Cell_WC(c,1).BoundingBox) > imG_thresh;
    tempG(~ProtL) = 0;
    Gprots = regionprops(ProtL, tempG, 'MaxIntensity');
    idx = find([Gprots.MaxIntensity]);
    ProtL = ismember(ProtL, idx);
    
    
    %pERM+ Protrusion stats
    [~, numProtsG] = bwlabel(ProtL);
    tempG = imcrop(imGreen, Cell_WC(c,1).BoundingBox);
    ProtsG = regionprops(ProtL,tempG, 'Area', 'MajorAxisLength', 'MaxIntensity', 'MeanIntensity');
    LengthG_ave = mean([ProtsG(:).MajorAxisLength]);
    LengthG_std = std([ProtsG(:).MajorAxisLength]);
    AreaG_ave = mean([ProtsG(:).Area]);
    AreaG_std = std([ProtsG(:).Area]);
    MaxIG_ave = mean([ProtsG(:).MaxIntensity]);
    MaxIG_std = std(double([ProtsG(:).MaxIntensity]));
    MeanIG_ave = mean([ProtsG(:).MeanIntensity]);
    MeanIG_std = std([ProtsG(:).MeanIntensity]);
    
    
    %% Record results in excel-friendly format
    cd(resdir)
    [cvf] = fopen(savefile,'a+');                                           % Open or create file for recording analysis results
    logger = sprintf('\n%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%fs', ...
        samplename , ProtperNuc, ProtperPerim/scale, ProtperNuc_alt, ProtperPerim_alt/scale, Cell_WC(c).Area , Cell_WC(c).Perimeter, ...
        Cell_AR, Nuc_stats.Area , Nuc_stats.Perimeter, Nuc_AR, numNuc, numProts, numProtsG, ...
        Area_ave, Area_std, AreaG_ave, AreaG_std, Length_ave, Length_std, LengthG_ave, LengthG_std, ...
        MaxIG_ave, MaxIG_std, MeanIG_ave, MeanIG_std, O, C);
    fwrite(cvf, logger);
    fclose(cvf);
  
end

    %% Record summary results in excel-friendly format
    cd(resdir)
    [cvf] = fopen(savefile_summary,'a+');                                           % Open or create file for recording analysis results
    logger = sprintf('\n%s\t%f\t%f\t%f\t%f', ...
        samplename , ProtperNuc, ProtperPerim/scale, ProtperNuc_alt, ProtperPerim_alt/scale);
    fwrite(cvf, logger);
    fclose(cvf);

%% Save image for QC afterwards
cd(tdir)

outfilename = sprintf('%s_mask.tif', savename);
temp_img = cat(3, im2uint8(prot_mask)*200, uint8(im_WCmask).*150, uint8(imB_mask).*100);
Save_Image(temp_img, resdir, outfilename, 5e8);

Poutline = bwperim(prot_mask);
NPoutline = bwperim(im_NPmask);
outfilename = sprintf('%s_mask_outline.tif', savename);
temp_img = cat(3, uint8(imRed), uint8(imGreen), uint8(Poutline | NPoutline).*100);
Save_Image(temp_img, resdir, outfilename, 5e8);

% save original as well
outfilename = sprintf('%s_raw.tif', savename);
Save_Image(RGBim, resdir, outfilename, 5e8);

% Return to original directory
cd(pdir)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Save_Image(img, resdir, outfilename, imsize)
% This m-files saves an RGB image as a tif to the current directory. If 
% the image is larger than imsize, it will reduce the image size.
% imsize = 5e8
pdir = cd;

cd(resdir)
info = whos('img');
if info.bytes>imsize
    tifscale = (imsize)/info.bytes;
else
    tifscale = 1;
end
imwrite(imresize(img,tifscale), outfilename, 'tif')
cd(pdir)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xPts, yPts] = interpPts(x1, x2, y1, y2)
% interpPts: Interpolate between points for integer x values
% This m-file is based on the Digitial Image Processing using MATLAB
% additional custom M-Functions in Appendix C of the textbook by Gonzalez,
% Woods, & Eddins.
% [xpts, ypts] = interpPts(x1, x2, y1, y2) computes an approximation to the
% line segment joining (x1, y1) & (x2, y2) with integer x-values. x1 & x2
% should be integers.

dx = abs(x2-x1);
dy = abs(y2-y1);

% Check for degenerate case.
if ((dx == 0) && (dy == 0))
    x=x1;
    y=y1;
    return;
end

flip = 0;
if (x1 > x2)
    % Always "draw" from left to right.
    t = x1; x1 = x2; x2 = t;
    t = y1; y1 = y2; y2 = t;
    flip = 1;
end

m = (y2 - y1)/(x2 - x1);
xPts = ((x1+1) : (x2-1)).';
yPts = y1 + m*(xPts - x1);

if (flip)
    xPts = flipud(xPts);
    yPts = flipup(yPts);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


