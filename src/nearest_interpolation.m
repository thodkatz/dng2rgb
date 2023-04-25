function Ccam = nearest_interpolation(input_image, pattern, M, N)
    Ccam = zeros(size(input_image));
    [r_mask,g_mask,b_mask] = cfa_masks(pattern, size(input_image));

    r = fill_nearest(input_image .* r_mask);
    g = fill_nearest(input_image .* g_mask);
    b = fill_nearest(input_image .* b_mask);

    Ccam(:,:,1) = r;
    Ccam(:,:,2) = g;
    Ccam(:,:,3) = b;

    % Always keep image clipped b/w 0-1
    Ccam = max(0,min(Ccam,1)); 
end

function array = fill_nearest(array)
    % assume that array comes from a filter channel of a cfa array
    [rows, cols] = size(array);
    if array(1,1) == 0 && array(2,2) ~= 0
        for row=2:2:rows
            for col=2:2:cols
                color_to_fill  = array(row,col);
                if col+1 <= cols
                    array(row,col+1) = color_to_fill;
                end
                if row+1 <= rows
                    array(row+1,col) = color_to_fill;
                end
                if row+1 <=  rows && col+1 <= cols
                    array(row+1,col+1) = color_to_fill;
                end
            end
            array(:,1) = array(:,2);
            array(1,:) = array(2,:);
            array(1,1) = array(2,2);
        end
    elseif array(1,1) ~= 0 && array(2,2) == 0
        for row=1:2:rows
            for col = 1:2:cols
                color_to_fill  = array(row,col);
                if col+1 <= cols
                    array(row,col+1) = color_to_fill;
                end
                if row+1 <= rows
                    array(row+1,col) = color_to_fill;
                end
                if row+1 <=  rows && col+1 <= cols
                    array(row+1,col+1) = color_to_fill;
                end
            end
        end
    elseif array(1,2) ~=0 && array(2,1) == 0
       for row=1:2:rows
            for col = 2:2:cols
                color_to_fill  = array(row,col);
                if col+1 <= cols
                    array(row,col+1) = color_to_fill;
                end
                if row+1 <= rows
                    array(row+1,col) = color_to_fill;
                end
                if row+1 <=  rows && col+1 <= cols
                    array(row+1,col+1) = color_to_fill;
                end
            end
       end
       array(:,1) = array(:,2);
    elseif array(2,1) ~= 0 && array(1,2) == 0
       for row=2:2:rows
            for col = 1:2:cols
                color_to_fill  = array(row,col);
                if col+1 <= cols
                    array(row,col+1) = color_to_fill;
                end
                if row+1 <= rows
                    array(row+1,col) = color_to_fill;
                end
                if row+1 <=  rows && col+1 <= cols
                    array(row+1,col+1) = color_to_fill;
                end
            end
       end
       array(1,:) = array(2,:);
    elseif (array(1,2) ~= 0 && array(2,1) ~= 0) || (array(1,1) ~=0 && array(2,2) ~= 0)
        for row=1:rows
            if array(1,2) ~= 0 && array(2,1) ~= 0
                if mod(row,2) == 0
                    start_col = 1;
                else
                    start_col = 2;
                end
            else
                if mod(row,2) ~= 0
                    start_col = 1;
                else
                    start_col = 2;
                end 
            end
            for col = start_col:2:cols
                color_to_fill  = array(row,col);
                if col+1<=cols
                    array(row,col+1) = color_to_fill;
                end
            end
        end
        array(:,1) = array(:,2);
    else
        error("Not implemented")
    end
    [new_rows, new_cols] = size(array);
    assert(new_rows == rows)
    assert(new_cols == cols)
end