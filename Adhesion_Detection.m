function Adhesion_Detection(varargin) 
%--------------------------------------------------------------------------
% Adhesion_Detection.m
% This m-file reads in a 3-channel image: DAPI, adhesion stain, and actin
% stain (or paxillin stain). First cells are segmented, then adhesions
% are segmented, and finally elongated focal adhesions are segmented.
%
%   Author: Ailey Crow
%   Most recent revision: Dec 5, 2012
%
%   Inputs:     - datadir: directory where raw data, results folder, etc. are
%                  located
%               - samplename: name of sample to be analyzed (as entered in
%                  datadir)
%   
%   Outputs:    None
%
%   Exports:    - _raw.tif: original image (low-res)
%               - _mask.tif: cell and elongated adhesion mask overlayed on
%                   adhesion channel
%               - _periph.tif: peripheral region mask & actin or 2ary
%                   adhesion mask overlayed on 2ary adhesion channel
%               - Adhesion_Detection.xls: quantification results
%
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
 
%% Setup directories & parameters
pdir = '~/Adhesion_Detection/';        % directory where program is saved
tdir = '~/Toolbox/';                           % directory of useful generic m-files
resdir = strcat(datadir, 'Results_Super');      % directory where output & parameters.mat are saved. 
savefile_AreaHist = 'Adhesion_Area_Hist.xls';   % additional output files
savefile_lengthHist = 'Ahdesion_Length_Hist.xls';

% load parameters
cd(resdir)
load('parameters.mat')


%% structuring elements
d100 = strel('disk', 50);
d10 = strel('disk', 5);
d3 = strel('square', 3);
d2 = strel('square', 2);


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


%% load ROI tiff image of tumor area, and resize to fit im
% This step is run only if a region of interest (ROI) has been previously created
if (uroi == 1)
    switch illumination
        case 'IF'   % For immunofluorescence acquisition
            roidir = strcat(datadir, '/rois/');
            if ~isempty(dir(char(strcat(roidir, samplename , '_roi.tif' ))))
                sroi = imread(char(strcat(roidir, samplename , '_roi.tif' )));
            else
                sroi = true(size(im,1), size(im,2));
            end
        case 'BF'   % for Brightfield acquisition
            ROI_regions = VectorAPI(samplename, 'GETALL');   % Load (any) regions from gSlide
            if ~isempty(ROI_regions)
                sroi = zeros(size(im,1), size(im,2)); % If ROI regions have been specified
                clear opts
                opts = struct();
                opts.MAGPOWER = lowmagpower;
                opts.LEFT = 0;           % Set boundaries for import to entire image
                opts.RIGHT = Xpix+1;
                opts.TOP = 0;
                opts.BOTTOM = Ypix+1;
                
                for j=1:numel(ROI_regions)
                    opts.REGIONID = ROI_regions(j).ID;
                    opts.MASK = 1;
                    [status, roi_mask] = CaptureAPI(samplename,opts);
                    sroi = (imfill(~im2bw(roi_mask), 'holes') & ~sroi ) | sroi;
                end
            else
                sroi = true(size(im,1), size(im,2));
            end
        otherwise
    end
    crop_mask = logical(im2bw(sroi));
    [r1, c1] = find(crop_mask);
    CropBox.BoundingBox(1) = min(c1);
    CropBox.BoundingBox(2) = min(r1);
    CropBox.BoundingBox(3) = max(c1) - min(c1);
    CropBox.BoundingBox(4) = max(r1) - min(r1);
    sroi = imcrop(sroi, CropBox.BoundingBox);
    switch illumination
        case 'IF'
            im = cat(3, imcrop(im(:,:,1), CropBox.BoundingBox), imcrop(im(:,:,2), CropBox.BoundingBox), imcrop(im(:,:,3), CropBox.BoundingBox), imcrop(im(:,:,4), CropBox.BoundingBox));
        case 'BF'
            im = imcrop(im, CropBox.BoundingBox);
    end
    clear crop_mask
    
else
    sroi = true(size(im,1), size(im,2));
end

%% Count nuclei
NucIm = mat2gray(im(:,:,DAPI_ch));                   % Select nuclear/DAPI channel
NucIm = imtophat(NucIm, d100);                  % top-hat filter to smooth intensities
dNuc = fspecial('disk', ceil(50*scale));           % create filter roughly the size of a nucleus
Nuc_mask = NucIm > 2*mean(NucIm(:));
Nuc_mask = imopen(Nuc_mask, d10);               % separate nuclei close together (replace with watershed if unsuccessful)
Nuc_mask = bwareaopen(imfill(Nuc_mask, 'holes'), ceil(500*scale^2));

% remove low intensity nucs
Nuc_mask = Nuc_mask & (NucIm > mean(NucIm(Nuc_mask))/2);
Nuc_mask = bwareaopen(imfill(Nuc_mask, 'holes'), ceil(5000*scale^2));

[~, Num_Nuc] = bwlabel(Nuc_mask);
%imshow(imoverlay(NucIm, bwperim(Nuc_mask), [0.3 1 0.3]));  % Preview mask - uncomment for troubleshooting

%% Select cells using actin channel
dcell = strel('disk', ceil(20*scale));          % scale of intracellular interest

actim = mat2gray(im(:,:,actin_ch));                  % Select actin channel
actim = imtophat(actim, d100);                  % top-halt filter to smooth intensities
Adim = mat2gray(im(:,:,Ad_ch));                      % Select adhesion channel
Adim = imtophat(Adim, d100);                    % top-halt filter to smooth intensities
Filt = fspecial('disk', ceil(10*scale));                     % filter to be used for segmentation

act_mask2 = actim > 1.2*mean2(actim);
act_mask = actim > imfilter(actim ,Filt, 'symmetric', 'conv');
act_mask = (actim > mean(actim(Nuc_mask))/5) & act_mask & act_mask2;
act_mask = imclose(act_mask, d10);

Ad_mask = Adim > imfilter(Adim, Filt, 'symmetric', 'conv');         % blur image & find peaks above blur
Testarea = imdilate(Nuc_mask, dcell);
Ad_mask = Ad_mask & (Adim > mean(Adim(Testarea)));      % /5

Cells = Nuc_mask | act_mask | Ad_mask; %Adim>0.9*mean2(Adim);               % Add adhesion channel to help fill in inside of cells

Cells = imclose(Cells, dcell);
Cells = ~bwareaopen(~Cells, 1000);
Cells = imopen(Cells, d3);
Cells = bwareaopen(Cells, ceil(10000*scale^2));
Cperim = bwperim(Cells);

% ensure cell obects have nuclei
CL = bwlabel(Cells, 4);
Nuc_test = regionprops(CL, Nuc_mask, 'MaxIntensity');
Cell_Int = [Nuc_test.MaxIntensity];
Cell_index = find(Cell_Int > 0);                % Select Cells objects with nuclei
Cells = ismember(CL, Cell_index);

Nuc_mask = Nuc_mask & Cells;

%imshow(imoverlay(imadjust(actim), bwperim(Cells), [0.3 1 0.3]))    % Preview mask - uncomment for troubleshooting

%% Select adhesions
Ad_mask = Adim > imfilter(Adim, Filt, 'symmetric', 'conv');         % blur image & find peaks above blur
Ad_mask = Ad_mask & (Adim > 0.75*mean(Adim(Cells)));

Adim(~Cells) = 0;
Ad_mask(~Cells) = 0;

Ad_mask = bwareaopen(Ad_mask, mag, 4);              % remove puncta with # pixels < amgnification, based on 4-connectivity
Ad_mask = bwmorph(Ad_mask, 'hbreak');

% figure; imshow(cat(3, im2uint8(imadjust(Adim)), uint8(Ad_mask)*100, Adim*0))  % Preview mask - uncomment for troubleshooting


%% Calculate adhesion characteristics & classify elongated adhesions
[AdL, numAd] = bwlabel(Ad_mask, 4);              % identify 4-connected objects (diagonal connections do not count)
Ad_prop = regionprops(AdL, Adim, 'Area', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Solidity', 'MeanIntensity');
Ad_Int = [Ad_prop.MeanIntensity];
Ad_Area = [Ad_prop.Area];
Ad_Sol = [Ad_prop.Solidity];
Ad_MaL = [Ad_prop.MajorAxisLength];
Ad_MiL = [Ad_prop.MinorAxisLength];  
Ad_Ecc = Ad_MaL./Ad_MiL;

Area_index = find(Ad_Area>mag);                 % Select larger adhesions
Ecc_index = find(Ad_Ecc>5);                   % Select long skinny adhesions
Int_index = find(Ad_Int > 1.5*mean2(Adim));     % Select adhsions of decent intensity
Sol_index = find(Ad_Sol > 0.3);                % Select adhesions that are not highly branched

Elong_maska = ismember(AdL, Area_index);
Elong_maske = ismember(AdL, Ecc_index);
Elong_masks = ismember(AdL, Sol_index);
Elong_maski = ismember(AdL, Int_index);

Elong_mask = Elong_maska & Elong_maske & Elong_maski & Elong_masks; 


[~, numElong] = bwlabel(Elong_mask);


%% Create histogram of # elongated adhesions per cell
% First separate cell clumps using watershed
W = watershed(bwdist(Nuc_mask));
WL = W==0;
CellObjects = Cells & ~WL;

% remove cell objects without nuclei
COL = bwlabel(CellObjects);
props = regionprops(COL, Nuc_mask, 'MeanIntensity');
NoNuc_idx = find([props.MeanIntensity] == 0);
NoNuc = ismember(COL, NoNuc_idx);
CellObjects = CellObjects & ~ NoNuc;

% Create mask of just centroids of elongated adhesions
Elong_Cent = false(size(CellObjects));
props = regionprops(Elong_mask, 'Centroid');
rc = [props.Centroid];
r = round(rc(1:2:numel(rc)));
c = round(rc(2:2:numel(rc)));
linear_indices = sub2ind(size(Elong_Cent), c, r);
Elong_Cent(linear_indices) = 1;

% Count elongated adhesions per cell object
COL = bwlabel(CellObjects);
props = regionprops(COL, _Cent, 'Area', 'MeanIntensity');
NumElong = [props.Area].*[props.MeanIntensity];
Hist = hist(NumElong, [0:20]);


%% Cell periphery calculations
Cperim = imdilate(Cperim, strel('disk', ceil(scale*5)));       % thicken border to include periphery of cell. Should be 10.
Cperim = Cperim & Cells;        % Restrict perimeter mask to intracellular area

% Now repeat for non-adhesion, non-nuclear channel (paxillin sometimes)
Pax_mask = actim > imfilter(actim, Filt, 'symmetric', 'conv');         % blur image & find peaks above blur
Pax_mask = Pax_mask & Cells;
Pax_mask = imopen(Pax_mask, d2);                                  % remove tiny peaks
CPPax = Cperim & Pax_mask;
[CPPaxL, numPPax] = bwlabel(CPPax);
numPPax = numPPax / sum(CPPax(:));     % normalize by peripheral area
CPPax_props = regionprops(CPPaxL, actim, 'Area', 'MeanIntensity');
PPax_Area = mean([CPPax_props.Area]);
PPax_Int = mean(actim(Cperim));


%% log stats
Cell_Area = sum(Cells(:));
Elong_Area = sum(Elong_mask(:));
MeanAdInt = mean(Adim(Cells));
MeanAdInt_mask = mean(Adim(Ad_mask));
Ad_Area = sum(Ad_mask(:));
PixA = scale^2;

Histout = num2str(Hist(1));
for h=2:numel(Hist)
    Histout = strcat(Histout, '^', num2str(Hist(h)));
end

cd(resdir)
[cvf1] = fopen(savefile,'a+');
fwrite(cvf1,sprintf('%s^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%s\n', char(savename), ...
    numElong, numAd, Cell_Area*PixA, Num_Nuc, MeanAdInt, MeanAdInt_mask, 100*Ad_Area / Cell_Area, 100*Elong_Area / Cell_Area, ...
    nuElong / Num_Nuc, numElong / (Cell_Area*PixA), mean(Ad_MaL), ... %PLave, PLstd, ...
    numPPax/PixA, PPax_Area, PPax_Int, Histout));
fclose(cvf1);

% Histogram of Adhesion Areas
[cvf1] = fopen(savefile_AreaHist,'a+');
fwrite(cvf1,sprintf('%s^%s\n', char(savename), num2str(Ad_Area)));
fclose(cvf1);

% Histogram of Adhesion lengths
[cvf1] = fopen(savefile_lengthHist,'a+');
fwrite(cvf1,sprintf('%s^%s\n', char(savename), num2str(Ad_MaL)));
fclose(cvf1);


%% Save image for QC afterwards
cd(tdir)
%Ad_perim = bwperim(Ad_mask);
Elong_perim = bwperim(Elong_mask);
Nuc_perim = bwperim(Nuc_mask);
Cperim = bwperim(Cells);
%Adim = mat2gray(Adim);              % convert double values to double values in range [0 1]
NucIm = imadjust(NucIm);
Adim = imadjust(Adim);              % adjust contrast
actim = imadjust(actim);

outfilename = sprintf('%s_Elong_mask.tif', savename);
temp_img = cat(3, im2uint8(Adim), uint8(Elong_perim).*100, uint8(Cperim | Nuc_perim).*100);
Save_Image(temp_img, resdir, outfilename, 5e8);

outfilename = sprintf('%s_Adhesion_mask.tif', savename);
temp_img = cat(3, im2uint8(Adim), uint8(Elong_mask).*100, uint8(Ad_mask).*100);
Save_Image(temp_img, resdir, outfilename, 5e8);

% Save peripheral adhesion immask
outfilename = sprintf('%s_periph.tif', savename);
temp_img = cat(3, im2uint8(imadjust(actim)), uint8(Pax_mask).*100, uint8(Cperim).*100);
Save_Image(temp_img, resdir, outfilename, 5e8);

% save original as well
% actim = imadjust(mat2gray(actim));
% NucIm = imadjust(mat2gray(NucIm));
outfilename = sprintf('%s_raw.tif', savename);
temp_img = cat(3, im2uint8(actim), im2uint8(Adim), im2uint8(NucIm));
Save_Image(temp_img, resdir, outfilename, 5e8);


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
