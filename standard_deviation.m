function [sigma_n] = standard_deviation(distorted_block_2d, error_mask_2d) 
    
    % compute variance
    undistorted_pixel = distorted_block_2d(error_mask_2d ~= 0);
    var_w = var(undistorted_pixel(:));

    % compute normalized standard deviation
    sigma_n = sqrt(var_w)/255;
    if (sigma_n < 0)
        sigma_n = 0;
    elseif (sigma_n > 1)
        sigma_n = 1;
    end
end
