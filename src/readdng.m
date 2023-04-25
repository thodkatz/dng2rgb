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