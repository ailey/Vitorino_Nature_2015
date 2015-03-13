function Pathological_Vasculature_Analysis(varargin) 


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
 
%% directories & parameters
pdir = '~/Adhesion_Detection/';        % directory where program is saved
tdir = '~/Toolbox/';                           % directory of useful generic m-files
Imdir = strcat(datadir, 'Images/');
resdir = strcat(datadir, 'Results/');
um2pix_XY = 2.486;          % um/pixel for 10x Confocal in xy

% load parameters
cd(resdir)
load('parameters.mat')

%% structuring elements
d100 = strel('disk', 50);
d40 = strel('disk', 20);
d30 = strel('disk', 14);
d20 = strel('disk', 10);
d14 = strel('disk', 7);
d10 = strel('disk', 5);
d7 = strel('disk', 4);
d3 = strel('square', 3);
d2 = strel('square', 2);

d5 = [ 0 1 1 1 0;
    1 1 1 1 1;
    1 1 1 1 1;
    1 1 1 1 1;
    0 1 1 1 0];
d5 = logical(d5);

%% Manage z-stack
% Cycle through each plane and create binary segmented image. Add binary slice to z stack
cd(Imdir), 
savename = regexp(samplename, '01c1.tif', 'split');     % create root of filename
savename = savename{1};
search4 =  strcat(savename, '*c', num2str(Vess_ch), '.tif');   % search for all z planes of vessel channel for given file
ImList = dir(search4);
ImList = [{ImList.name}];
zsteps = numel(ImList);
im = imread(char(ImList(1)));
Ret_Mask3D = false(size(im,1), size(im,2), zsteps);         % initialize 3D binary matrix for segmented z stack
Vess_Mask3D = false(size(im,1), size(im,2), zsteps);         % initialize 3D binary matrix for segmented z stack

for z = 1:zsteps
    im = tofloat(imread(char(ImList(z))));           % read-in z slice & make floating point
    im = imtophat(im, d100);
    Ret_mask = imfill(im > 0.001, 'holes');                     % global threshold to segment for retina mask
    
    Filt = fspecial('disk', 10);                     % filter to be used for segmentation
    V_mask1 = im > imfilter(im, Filt, 'symmetric', 'conv');         % blur image & find peaks above blur 
    V_mask2 = im > mean2(im(Ret_mask))/3;           % identify regions 3x brighter than average signal in retina
    V_mask = V_mask1 & V_mask2 & Ret_mask;    
    im(~V_mask) = 0;
    
    Ret_Mask3D(:,:,z) = Ret_mask;
    Vess_Mask3D(:,:,z) = im;
end

%% Create z projections
Ret_ZP = sum(Ret_Mask3D,3)/zsteps;
Vess_ZP = sum(Vess_Mask3D,3)/zsteps;


%% Find retinal area
Ret_mask = Vess_ZP>mean2(Vess_ZP);      % global threshold for all significant signal (threshold level set to mean of image)
Ret_mask = imclose(Ret_mask, d20);
Ret_mask = imfill(Ret_mask, 'holes');   % fill in holes in retina mask
Ret_mask = imopen(Ret_mask, d10);
Ret_mask = bwareaopen(Ret_mask, 1e5);   % remove small external objects

% remove optic nerve
Cent = regionprops(Ret_mask, 'Centroid');
CentX = Cent.Centroid(1);
CentY = Cent.Centroid(2);
Cent_mask = false(size(Vess_ZP));
Cent_mask(CentY-199:CentY+199,CentX-199:CentX+199) = getnhood(strel('disk', 200));
ON_mask = Vess_ZP;
ON_mask(~Cent_mask) = 0;
ON_mask = ON_mask > 1.5*mean2(Vess_ZP(Cent_mask));
ON_mask = imopen(ON_mask, d5);
ON_mask = imfill(ON_mask, 'holes');
ON_mask = imopen(ON_mask, d10);
ON_mask = bwareaopen(imclose(ON_mask, d10), 500);
Ret_mask(ON_mask) = 0;

rp = imdilate(bwperim(Ret_mask), d3);
%imshow(imoverlay(Vess_ZP, rp, [0.3 1 0.3]))    %Uncomment to preview for troubleshooting


%% Select features (all vasculature)
Vim = imtophat(Vess_ZP, d100);          % large-scale smoothing to dampen uneven illumination (due to hemorraging etc.)

Filt = fspecial('disk', 100);                     % filter to be used for segmentation
V_mask2 = Vim > imfilter(Vim, Filt, 'symmetric', 'conv');         % blur image & find peaks above blur 
V_mask2 = bwareaopen(V_mask2, 50,8);

Filt = fspecial('disk', 3);                     % filter to be used for segmentation
V_mask3 = Vim > imfilter(Vim, Filt, 'symmetric', 'conv');         % blur image & find peaks above blur 
V_mask3 = bwareaopen(V_mask3, 50,8);
V_mask = (V_mask2 | V_mask3) & Ret_mask;


%% Remove normal vessels to select aneurysms only
An_mask = imopen(V_mask, d5);
An_mask = bwareaopen(An_mask, 50); 

% Find aneurysm areas by another method
V_mask1 = Vim > mean2(Vim);
Filt = fspecial('disk', 25);
Holes = Vim < (imfilter(Vim, Filt, 'symmetric', 'conv')/2);
V_mask1 = V_mask1 & ~Holes & Ret_mask;
V_mask1 = imopen(V_mask1, d10);
V_mask1 = ~bwareaopen(~V_mask1, 10);

An_mask = An_mask | V_mask1;

% remove low intensity areas from aneurysm mask (to isolate vessels vs
% aneurysms)
VimA = Vim;
VimA(~An_mask) = 0;
An_mask = VimA > mean(Vim(V_mask));

% Select aneurysms based on texture
std = stdfilt(Vim, getnhood(d3));
std(~An_mask) = 0;
An_mask2 = (std < 12) & (std > 0) & Ret_mask;
An_mask2 = imopen(An_mask2, d2);
An_mask2 = bwareaopen(An_mask2, 50);

% Remove large normal vessels based on aspect ratio
Al = bwlabel(An_mask);
Aprops = regionprops(An_mask, 'MajorAxisLength', 'MinorAxisLength');
AR = [Aprops.MajorAxisLength] ./ [Aprops.MinorAxisLength];
idx = find(AR>4);       
An_normal = ismember(Al, idx);


% % Identify large pathological areas (likely to be real)
Al = bwlabel(An_mask);
Aprops = regionprops(An_mask, 'Area');
idx = find([Aprops.Area] > 5e4);
An_large = ismember(Al, idx);

% Find areas that are barely separated & join
An_bits = imopen(V_mask, d40);
An_bits = bwareaopen(An_bits, 7000);
An_bits = imclose(An_bits, d3);
An_bits = bwareaopen(An_bits, 3e4);

% Identify smaller aneurysms that are either circular or solid (fewer
% protrusions like a vessel). 
Al = bwlabel(An_mask);
Aprops = regionprops(An_mask, 'Area', 'Solidity', 'Extent', 'MajorAxisLength', 'MinorAxisLength');
idx = find([Aprops.Solidity] > 0.5);
An_mask1 = ismember(Al, idx);
idx = find([Aprops.Extent] > 0.3);
An_mask2 = ismember(Al, idx);
AR = [Aprops.MajorAxisLength] ./ [Aprops.MinorAxisLength];
idx = find(AR<3);
An_mask3 = ismember(Al, idx);
An_mask = An_mask2 & An_mask3;

% Select aneurysms based on texture
Al = bwlabel(An_mask);
E1 = edge(Vess_ZP, 'Sobel', 0.02);
Aprops = regionprops(An_mask, E1, 'MeanIntensity');
idx = find([Aprops.MeanIntensity] > 0.01);
An_mask = ismember(Al, idx);

An_mask = bwareaopen(((An_mask | An_bits | An_large) & ~An_normal), 25);        


%% Calculate vascular Coverage
Vasc = imclose(V_mask, d100);
Vasc_Cov = sum(Vasc(:)) / sum(Ret_mask(:));


%% Save masks for troubleshooting
cd(tdir)
outfilename = sprintf('%s_mask.tif', savename);
temp_img = cat(3, im2uint8(bwperim(An_mask)).*200, im2uint8(imadjust(Vess_ZP)), im2uint8(bwperim(V_mask)).*200);
Save_Image(temp_img, resdir, outfilename, 5e8);


%% Save aneurysm ROIs for future manipulation (i.e. manual correction)
cd(roidir)
roiname = sprintf('%s_roi.tif', savename);
imwrite(An_mask, roiname, 'tif', 'Compression', 'lzw');
save(strcat(savename, '.mat'), 'Vess_ZP', 'An_mask', 'V_mask', 'Ret_mask')


%% log stats
Ret_Area = sum(Ret_mask(:));
Vess_Area = sum(V_mask(:));
An_Area = sum(An_mask(:));
AS = um2pix_XY^2;

cd(resdir)
[cvf1] = fopen(savefile,'a+');
fwrite(cvf1,sprintf('%s^%f^%f^%f^%f^%f^%f^%f\n', char(savename), ...
    Ret_Area*AS, Vess_Area*AS, An_Area*AS, 100*An_Area/Vess_Area, 100*An_Area/Ret_Area, 100*Vess_Area/Ret_Area, Vasc_Cov*100));
fclose(cvf1);

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

