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