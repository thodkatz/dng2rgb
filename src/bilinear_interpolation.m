function Ccam = bilinear_interpolation(input_image, pattern, M, N)
    Ccam = zeros(size(input_image));
    [r_mask,g_mask,b_mask] = cfa_masks(pattern, size(input_image));
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


