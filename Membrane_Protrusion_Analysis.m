
function Membrane_Protrusion_Analysis(varargin)
%--------------------------------------------------------------------------
% Membrane_Protrusion_Analysis
% This m-file reads in a confocal image of the retinal front with vascular stain. The
% vasculature is segmented & then the front is examined for 2-classes
% of small protrusions: veil (~width of a single endothelial nucleus) &
% tip-cells: slightly larger. These protrusions are counted & normalized by
% vascular front perimeter.
%
%
%   Author: Ailey Crow
%   Date created: April 1 2013
%   Most recent revision: June 12 2013 (adjusted for confocal images)
%
%   Inputs:     - datadir: directory where raw data, results folder, etc. are
%                  located
%               - samplename: name of sample to be analyzed (as entered in
%                  datadir)
%               - Param: structure contatining analysis parameters
%
%   Outputs:    None
%
%   Exports:    - _raw.tif: original image
%               - _mask.tif: cell and fibrillar adhesion mask overlayed on
%                   adhesion channel
%               - ERG_Headskin_quant.xls: quantification results
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


%% Setup
%Select relevant directories
pdir = '~/Adhesion_Detection/';        % directory where program is saved
tdir = '~/Toolbox/';                           % directory of useful generic m-files
resdir = strcat(datadir, 'Results');        % directory where output & parameters.mat are saved. 

% load parameters
cd(resdir)
load('parameters.mat')

%% loop through all images
for i = 1:numel(List2Anal)
    samplename = char(List2Anal(i));
    
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
    
    V_im = im(:,:,Vess_ch);
    
    %% load ROI tiff image of tumor area, and resize to fit im
    % This step is run only if a region of interest (ROI) has been previously created
    if (uroi == 1)
        switch illumination
            case 'IF'   % immunofluorescence acquisition
                roidir = strcat(datadir, '/rois/');
                if ~isempty(dir(char(strcat(roidir, samplename , '_roi.tif' ))))
                    sroi = imread(char(strcat(roidir, samplename , '_roi.tif' )));
                    %sroi = imresize(sroi, [ size(im, 1) size(im,2)], 'nearest' );
                else
                    sroi = true(size(im,1), size(im,2));
                end
            case 'BF'   % Brightfield acquisition
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
    
    
    %% Segment whole retina
    V_imMF = medfilt2(V_im, [4 4], 'symmetric');        % median filter to remove speckle noise
    Ret_mask1 = V_imMF >  imfilter(V_imMF, fspecial('Gaussian', 15, 1), 'symmetric');        % local thresholding
    Ret_mask1 = bwareaopen(Ret_mask1, 10, 4);
    Ret_mask2 = V_imMF > 1*mean(V_imMF(Ret_mask1));      % factor for inducible KO rescue study: 0.7. For rescue repeat: 0.5, For smi: 0.9 or 1
    Ret_mask_noFill = bwareaopen(Ret_mask1 | Ret_mask2, 5000, 4);
    Ret_mask = imfill(Ret_mask_noFill, 'holes');
    Ret_mask = imopen(Ret_mask, strel('disk', 1));
    Ret_mask = bwareaopen(Ret_mask, 5000);
    
    %% Border operations
    Ret_mask = padarray(Ret_mask, [1 1], 1);        % Pad with border of true pixels
    Ret_mask = ~bwareaopen(~Ret_mask, 100000);  % Fill in holes touching border
    Ret_mask = Ret_mask(2:(size(Ret_mask,1)-1), 2:(size(Ret_mask,2)-1));        % remove padding
    
    TouchBorder = Ret_mask & bwperim(ones(size(Ret_mask,1),size(Ret_mask,2)));  % Create mask of pixels where Ret_mask touches border
    
    %% Compute perimeter parameters
    Rperim = bwareaopen(bwperim(Ret_mask) & ~TouchBorder, 50);              % create perimeter mask without border pixels
    Vperim = sum(sum(Rperim));            % Perimeter of vasculature
    Rperim_flat = imclose(imopen(Ret_mask, strel('disk', ceil(50*scale))), strel('disk', ceil(100*scale)));
    Rperim_flat = imclose(Rperim_flat, strel('square', ceil(300*scale)));
    TouchBorder = Rperim_flat & bwperim(ones(size(Ret_mask,1),size(Ret_mask,2)));  % Create mask of pixels where Ret_mask touches border
    Rperim_flat = bwareaopen(bwperim(Rperim_flat) & ~TouchBorder, 100); % create rough perimeter mask without border pixels
    Vperim_flat = sum(sum(Rperim_flat));       % Equivalent to hand-drawn general perimeter calculation
    
    %% Segment veils & tip cells
    
    % Identify veils
    Veil = Ret_mask & ~imopen(Ret_mask, strel('disk', ceil(7.5*scale)));         % Classify veil objects as thickness < 15um (15um is size of single nucleus)
    Veil = imopen(Veil, strel('disk', 1));      % Remove accessory junk
    Veil = bwareaopen(Veil, 50);                % size exclusion
    Veil = imclearborder(Veil);
    
    % Check intensity
    VL = bwlabel(Veil);
    V_prop = regionprops(VL, V_im, 'MajorAxisLength', 'MinorAxisLength', 'MeanIntensity');
    V_I = [V_prop.MeanIntensity];
    V_MaL = [V_prop.MajorAxisLength];
    V_MiL = [V_prop.MinorAxisLength];
    V_Ecc = V_MaL./V_MiL;
    Ecc_index = find(V_Ecc>2.0);                   % Select long skinny protrusions. 2.5 good cutoff
    Veil_E = ismember(VL, Ecc_index);
    I_index = find(V_I > 0.95*mean(V_im(:)));   % changed from 1.1
    Veil_I = ismember(VL, I_index);
    MaL_index = find(V_MaL>40);                   % Select long protrusions. 40pixel good cutoff
    Veil_L = ismember(VL, MaL_index);
    Veil = Veil_E & Veil_I & Veil_L;
    
    % Check endpoints
    EP = bwmorph(bwmorph(Ret_mask, 'thin', Inf), 'endpoints');
    Veil = imreconstruct(EP, Veil);
    
    % Count
    [~,nVeil] = bwlabel(Veil);
    Perim_NoVeil = sum(sum(Rperim & ~Veil));
    
    % Identify tip cell
    TipCell = Ret_mask & ~imopen(Ret_mask, strel('disk', ceil(15*scale)));         % Classify tip cell candidates as protruding 30um (15um is size of single nucleus)
    TipCell = imopen(Ret_mask_noFill& TipCell & ~Veil, strel('disk', ceil(7.5*scale)));
    TipCell = imclearborder(TipCell);
    
    % Check intensity
    TL = bwlabel(TipCell);
    T_prop = regionprops(TL, V_im, 'MeanIntensity');
    T_I = [T_prop.MeanIntensity];
    I_index = find(T_I > 1.1*mean(V_im(:)));
    TipCell = ismember(TL, I_index);
    
    [~,nTipCell] = bwlabel(TipCell);
    
    
    %% Record results in excel format
    cd(resdir)
    [cvf1] = fopen(savefile,'a+');
    fwrite(cvf1,sprintf('%s^%f^%f^%f^%f^%f\n', char(savename), ...
        scale*Vperim, scale*Vperim_flat, scale*Perim_NoVeil, nVeil, nTipCell));
    fclose(cvf1);

    % Veil metrics
    V_prop = regionprops(Veil, 'MajorAxisLength', 'MinorAxisLength');
    V_MaL = [V_prop.MajorAxisLength];
    V_MiL = [V_prop.MinorAxisLength];
    V_Ecc = V_MaL./V_MiL;
    
    [cvf1] = fopen('Veil_Metrics.xls','a+');
    for o = 1:numel(V_prop)
        fwrite(cvf1,sprintf('%s^%f^%f^%f^%f^%f\n', char(savename), ...
            scale*Vperim, scale*Vperim_flat, scale*Perim_NoVeil, V_Ecc(o), V_MaL(o)));
    end
    fclose(cvf1);
   
    
    %% Save image for QC afterwards
    cd(tdir)
    
    outfilename = sprintf('%s_mask.tif', savename);
    tmp = imoverlay(V_im, bwperim(TipCell), [0.3 1 0.3]);
    temp_img = imoverlay(tmp, bwperim(Veil), [1 0.3 0.3]);
    %temp_img = cat(3, im2uint8(Veil)*150, im2uint8(TipCell).*100, im2uint8(Ret_mask_noFill).*100);
    Save_Image(temp_img, resdir, outfilename, 5e8);
    
    outfilename = sprintf('%s_perim.tif', savename);
    temp_img = cat(3, im2uint8(Rperim_flat)*100, im2uint8(Rperim)*100, V_im);
    Save_Image(temp_img, resdir, outfilename, 5e8);
    
end

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

