% SOURCE:
% RobSumner,Processing RAW Images in MATLAB,https://rcsumner.net/raw_guide/RAWguide.pdf,May19,2014

% setup the path to include the 'utils' directory
directory = pwd;
addpath(genpath(directory))

% SPECIFY THE FILEPATH OF THE DNG FILE !!!
root = pwd + "/../";
filename = root + "assignment1/assets/RawImage.DNG";

% PREPROCESSING
[rawim, XYZ2Cam, wbcoeffs] = readdng(filename);


% BAYER TYPE
% bayertype = "gbrg";
bayertype = "rggb";
% bayertype = "bggr";
% bayertype = "grbg";

% unused, not implemented
m = 4000;
n = 6000;

% INTERPOLATION METHOD
interpolation_methods = ["linear", "nearest"];
% bayertypes = ["rggb", "gbrg", "bggr", "grbg"];
bayertypes = ["rggb"];
for method = interpolation_methods
    for bayertype = bayertypes
        [Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, bayertype, method, m, n);

        % PLOT
        file_to_save = "Clinear" + "_" + method + "_" + bayertype + ".jpg";
        imwrite(Clinear, file_to_save)

        f = figure;
        hold on
        r = Clinear(:,:,1);
        g = Clinear(:,:,2);
        b = Clinear(:,:,3);
        histogram(r(:), 100, "FaceColor","r")
        histogram(g(:), 100, "FaceColor","g")
        histogram(b(:), 100, "FaceColor","b")
        xlabel("Intensity")
        ylabel("Frequency")
        file_to_save = "Clinear" + "_" + method + "_" + bayertype + "_" + "hist" + ".jpg";
        saveas(f, file_to_save)

        file_to_save = "Csrgb" + "_" + method + "_" + bayertype + ".jpg";
        imwrite(Csrgb, file_to_save)

        f = figure;
        hold on
        r = Csrgb(:,:,1);
        g = Csrgb(:,:,2);
        b = Csrgb(:,:,3);
        histogram(r(:), 100, "FaceColor","r")
        histogram(g(:), 100, "FaceColor","g")
        histogram(b(:), 100, "FaceColor","b")
        xlabel("Intensity")
        ylabel("Frequency")
        file_to_save = "Csrgb" + "_" + method + "_" + bayertype + "_" + "hist" + ".jpg";
        saveas(f, file_to_save)
    end
end