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
