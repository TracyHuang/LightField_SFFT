function [image_result basis] = RgbYuvConversion(rgb_image, sources, targets, source_reference, target_reference)

    W_R = 0.299; 
    W_B = 0.114;
    W_G = 1 - W_R - W_B;
    U_m = 0.436;
    V_m = 0.615;
    
    b_r = [1 0 0];
    b_g = [0 1 0];
    b_b = [0 0 1];
    
    b_y = [W_R, W_G, W_B];
    b_u = U_m / (1 - W_B) * (b_b - b_y);
    b_v = V_m / (1 - W_R) * (b_r - b_y);

    num_targets = length(targets);
    num_results = num_targets;
    if num_results == 2
        num_results = 3;
    end
    basis = zeros(num_results, 3);
    
    num_sources = length(sources);
    if (num_sources < 3)
        fprintf(2, 'num_sources = %d, should be at least 3\n', num_sources);
    end
    
    if ~exist('target_reference', 'var')
        target_reference = [];
    end
    
    if ~exist('source_reference', 'var')
        source_reference = [];
    end
    
    for i = 1 : num_targets
        [basis(i, :) target_reference] = getBasis(targets(i), b_r, b_g, b_b, b_y, b_u, b_v, target_reference);
    end
    
    source_basis = zeros(num_sources, 3);
    for i = 1 : num_sources
        [source_basis(i, :) source_reference] = getBasis(sources(i), b_r, b_g, b_b, b_y, b_u, b_v, source_reference);
    end
    
    x = size(rgb_image, 1);
    y = size(rgb_image, 2);
        
    image_result = zeros(x, y, num_results);
    if (num_targets == 2)
        A = [basis(1 : 2, :); 1 1 1];
        third_one = A^(-1)*[0; 0; 1];
        basis(3, :) = third_one(:).';
    end

    coeff = (basis / source_basis);
    for i = 1 : x
        for j = 1 : y
            vec = coeff * double(reshape(rgb_image(i, j, :), [3 1]));
            image_result(i, j, :) = reshape(vec, [1 1 num_results]);
        end
    end    

end


function [vec reference] = getBasis(basis, b_r, b_g, b_b, b_y, b_u, b_v, reference) 

    switch lower(basis)
        case 'r'
            vec = b_r;
        case 'g'
            vec = b_g;
        case 'b'
            vec = b_b;
        case 'y'
            vec = b_y;
        case 'u'
            vec = b_u;
        case 'v'
            vec = b_v;
        otherwise
            if size(reference, 1) > 0
                vec = reference(1, :);
                reference = reference(2 : end, :);
            else
                fprintf(2, 'Unknown basis %s\n', basis);
            end
    end

end
