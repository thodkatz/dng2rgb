% SOURCE:
% RobSumner,Processing RAW Images in MATLAB,https://rcsumner.net/raw_guide/RAWguide.pdf,May19,2014

% setup the path to include the 'utils' directory
directory = pwd;
addpath(genpath(directory))

root = pwd + "/../";
filename = root + "assignment1/assets/RawImage.DNG";
[rawim, XYZ2Cam, wbcoeffs] = readdng(filename);
% bayertype = "gbrg";
bayertype = "rggb";
% bayertype = "bggr";
% bayertype = "grbg";
method = "nearest";
m = 4000;
n = 6000;
[Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, bayertype, method, m, n);

% figure
% title("Clinear")
% imshow(Clinear)
% 
% figure
% title("Cxyz")
% imshow(Cxyz)
% 
% figure
% title("Ccam")
% imshow(Ccam)

figure
title("Csrgb")
imshow(Csrgb)


function [Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, bayertype, method, M, N)
    % preprocessed raw image (normalized croped raw image)

    % white balancing
    [m, n] = size(rawim);
    mask = wbmask(m,n, wbcoeffs,bayertype);
    balanced_bayer = rawim .* mask;

    % interpolate
    switch method
        case "linear"
            Ccam = bilinear_interpolation(balanced_bayer, bayertype, M, N);
        case "nearest"
            Ccam = nearest_interpolation(balanced_bayer, bayertype, M, N);
    end

    % color space conversion
    xyz2rgb = [3.2406, -1.5372, -0.4986; -0.9689, 1.8758, 0.0415; 0.0557, -0.2040, 1.0570];
    rgb2cam = XYZ2Cam * xyz2rgb^-1;
    rgb2cam = rgb2cam ./ repmat(sum(rgb2cam,2),1,3); % normalize rows to 1, explained by the paper (it is a trick)

    Cxyz = apply_cmatrix(Ccam, inv(XYZ2Cam));
    Clinear = apply_cmatrix(Ccam, inv(rgb2cam));

    % Always keep image clipped b/w 0-1
    Cxyz = max(0,min(Cxyz,1)); 
    Clinear = max(0,min(Clinear,1));

    % gamma correction
    Csrgb = Clinear.^(1/2.2);
end

function corrected = apply_cmatrix(im,cmatrix)
    % CORRECTED = apply_cmatrix(IM,CMATRIX)
    %
    % Applies CMATRIX to RGB input IM. Finds the appropriate weighting of the
    % old color planes to form the new color planes, equivalent to but much
    % more efficient than applying a matrix transformation to each pixel.

    if size(im,3) ~= 3
        error("Apply cmatrix to RGB image only.")
    end
    r = cmatrix(1,1)*im(:,:,1)+cmatrix(1,2)*im(:,:,2)+cmatrix(1,3)*im(:,:,3);
    g = cmatrix(2,1)*im(:,:,1)+cmatrix(2,2)*im(:,:,2)+cmatrix(2,3)*im(:,:,3);
    b = cmatrix(3,1)*im(:,:,1)+cmatrix(3,2)*im(:,:,2)+cmatrix(3,3)*im(:,:,3);
    corrected = cat(3,r,g,b);
end

function Ccam = bilinear_interpolation(input_image, pattern, M, N)
    Ccam = zeros(size(input_image));
    [r_mask,g_mask,b_mask] = cfa_masks(pattern, [M, N]);
    g_filter = [0,1,0;1,4,1;0,1,0]/4;
    rb_filter = [1,2,1;2,4,2;1,2,1]/4;
    r = conv2(input_image .* r_mask, rb_filter, "same");
    b = conv2(input_image .* b_mask, rb_filter, "same");
    g = conv2(input_image .* g_mask, g_filter, "same");

    Ccam(:,:,1) = r;
    Ccam(:,:,2) = g;
    Ccam(:,:,3) = b;

    % Always keep image clipped b/w 0-1
    Ccam = max(0,min(Ccam,1)); 
end

function [r, g, b] = cfa_masks(pattern, cfa_size)
    r = false(cfa_size);
    b = false(cfa_size);
    g = false(cfa_size);
    starting_coordinates = {[0,0], [0,1], [1,0], [1,1]};
    pattern = char(pattern);
    for i = 1:4
        color_pattern = pattern(i);
        y = starting_coordinates{i}(1) + 1;
        x = starting_coordinates{i}(2) + 1;
        switch color_pattern
            case "r"
                r(y:2:end, x:2:end) = true;
            case "g"
                g(y:2:end, x:2:end) = true;
            case "b"
                b(y:2:end, x:2:end) = true;
        end
    end
end

function Ccam = nearest_interpolation(input_image, pattern, M, N)
    Ccam = zeros(size(input_image));
    [r_mask,g_mask,b_mask] = cfa_masks(pattern, [M N]);

    r = fill_nearest(input_image .* r_mask);
    g = fill_nearest(input_image .* g_mask);
    b = fill_nearest(input_image .* b_mask);

    Ccam(:,:,1) = r;
    Ccam(:,:,2) = g;
    Ccam(:,:,3) = b;

    % Always keep image clipped b/w 0-1
    Ccam = max(0,min(Ccam,1)); 
end


function [rawim, XYZ2Cam, wbcoeffs] = readdng(filename)
    % extract the color filter array (CFA)
    warning('off', 'imageio:tiffmexutils:libtiffWarning')
    obj = Tiff(filename, "r");
    offsets = getTag(obj,"SubIFD");
    setSubDirectory(obj, offsets(1));
    rawim = read(obj);
    [m0,n0] = size(rawim);
    close(obj);

    % extract metadata
    meta_info = imfinfo(filename);
    subIFDs = meta_info.SubIFDs{1};

    % Crop to only valid pixels
    % todo: try without cropping
    x_origin = meta_info.SubIFDs{1}.ActiveArea(2)+1; % +1 due to MATLAB indexing
    width = meta_info.SubIFDs{1}.DefaultCropSize(1);
    y_origin = meta_info.SubIFDs{1}.ActiveArea(1)+1;
    height = meta_info.SubIFDs{1}.DefaultCropSize(2);
    fprintf("Size of raw image MXN: [%d,%d]\n", m0,n0);
    rawim = double(rawim(y_origin:y_origin+height-1,x_origin:x_origin+width-1));
    [m0,n0] = size(rawim);
    fprintf("Size of cropped raw image MXN: [%d,%d]\n", m0,n0);

    rawim = linearization(rawim, subIFDs);

    wbcoeffs = (meta_info.AsShotNeutral).^-1;
    wbcoeffs = wbcoeffs/wbcoeffs(2); % green channel will be left unchanged

    XYZ2Cam = meta_info.ColorMatrix2;
    XYZ2Cam = reshape(XYZ2Cam,3,3)';
end

function rawim = linearization(rawim, subIFDs)
    % It is possible that the camera applied a non-linear 
    % transformation to the sensor data for storage purposes (e.g., Nikon cameras)
    % If the values are stored non-linearly, undo that mapping.
    if isfield(subIFDs,"LinearizationTable")
        disp("Non linear transformation detected...")
        ltab=meta_info.SubIFDs{1}.LinearizationTable;
        rawim = ltab(rawim+1);
    end

    % normalization given black and white levels
    blackLevel = subIFDs.BlackLevel;
    whiteLevel = subIFDs.WhiteLevel;

    rawim = (rawim-blackLevel)/(whiteLevel-blackLevel);
    rawim = max(0,min(rawim,1)); % trim noisy values out of the range [blackLeve, whiteLevel]
end

function colormask = wbmask(m,n,wbmults,align)
    % COLORMASK = wbmask(M,N,WBMULTS,ALIGN)
    %
    % Makes a white-balance multiplicative mask for an image of size m-by-n
    % with RGB while balance multipliers WBMULTS = [R_scale G_scale B_scale].
    % ALIGN is string indicating Bayer arrangement: ’rggb’,’gbrg’,’grbg’,’bggr’
    colormask = wbmults(2)*ones(m,n); %Initialize to all green values
    switch align
        case "rggb"
            colormask(1:2:end,1:2:end) = wbmults(1); %r
            colormask(2:2:end,2:2:end) = wbmults(3); %b
        case "bggr"
            colormask(2:2:end,2:2:end) = wbmults(1); %r
            colormask(1:2:end,1:2:end) = wbmults(3); %b
        case "grbg"
            colormask(1:2:end,2:2:end) = wbmults(1); %r
            colormask(1:2:end,2:2:end) = wbmults(3); %b
        case "gbrg"
            colormask(2:2:end,1:2:end) = wbmults(1); %r
            colormask(1:2:end,2:2:end) = wbmults(3); %b
    end
end