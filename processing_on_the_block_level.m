function extrapolated_block = processing_on_the_block_level(non_reg_distorted_block, block_error_mask, fse_parameters, rho, normalised_standard_div)
    
    % parameters
    fft_size = fse_parameters.fft_size;
    orthogonality_correction = fse_parameters.orthogonality_correction;
    [M, N] = size(non_reg_distorted_block);
    fft_x_offset = floor((fft_size-N)/2);
    fft_y_offset = floor((fft_size-M)/2);

    % weighting function
    w = zeros(fft_size);
    w(fft_y_offset+(1:M), fft_x_offset+(1:N)) = block_error_mask;
    for u = 0:fft_size-1
        for v = 0:fft_size-1
            w(u+1,v+1) = w(u+1,v+1) * rho^(sqrt((u+0.5-(fft_y_offset+M/2))^2 + (v+0.5-(fft_x_offset+N/2))^2));
        end
    end
   
    W = fft2(w);
    %to ensure that we can obtain negative frequencies
    W_overal = [W, W; W, W];
    
    % frequency weighting
    decreasing_propability = ones(fft_size, fft_size/2+1);
    for y=0:fft_size-1
        for x=0:fft_size/2
            y2 = fft_size/2 - abs(y - fft_size/2);
            x2 = fft_size/2 - abs(x - fft_size/2);
            decreasing_propability(y+1, x+1) = 1 - sqrt(x2*x2 + y2*y2)*sqrt(2)/fft_size;
        end
    end

    % pad image to fft window size
    f = zeros(fft_size);
    f(fft_y_offset+(1:M), fft_x_offset+(1:N)) = non_reg_distorted_block;

    % create initial model
    G = zeros(fft_size);

    % calculate initial residual
    Rw = fft2(f.*w);
    % to make sure that we search for the basis functions in the right area
    Rw = Rw(1:fft_size, 1:fft_size/2+1);
    
    % estimating the ideal number of iterrations depending on the standard
    % div value
    if (normalised_standard_div == 0)
        normalised_standard_div = standard_deviation(non_reg_distorted_block, block_error_mask) ;
    end
    num_iters = round(fse_parameters.iter_const * normalised_standard_div);
    if (num_iters < fse_parameters.min_iter)
        num_iters = fse_parameters.min_iter;
    elseif (num_iters > fse_parameters.max_iter)
		num_iters = fse_parameters.max_iter;
    end
    % begin extrapolation
    iter_counter = 0;
    while (iter_counter < num_iters) 
        projection_distances = abs(Rw(:)) .* decreasing_propability(:);
        [~, selected_Basis_function_location] = max(projection_distances);% select the basis function which maximize projection distance
        selected_Basis_function_location = selected_Basis_function_location(1)-1;
        % the indexes of the selected basis function
        v = floor(selected_Basis_function_location/fft_size);
        u = mod(selected_Basis_function_location,fft_size);
        
        % exclude second half of first and middle col
        if (v == 0 && u > fft_size/2 || v == fft_size/2 && u > fft_size/2)
            u_prev = u;
            u = fft_size-u;
            Rw(u+1,v+1) = conj(Rw(u_prev+1,v+1));
        end

        % calculate complex conjugate solution
        u_conj = -1; v_conj = -1;
        % fill first lower col (copy from first upper col)
        if (u >= 1 && u < fft_size/2 && v == 0)
            u_conj = fft_size-u;
            v_conj = v;
        end
        % fill middle lower col (copy from first middle col)
        if (u >= 1 && u < fft_size/2 && v == fft_size/2)
            u_conj = fft_size-u;
            v_conj = v;
        end
        % fill first row right (copy from first row left)
        if (u == 0 && v >= 1 && v < fft_size/2)
            u_conj = u;
            v_conj = fft_size-v;
        end
        % fill middle row right (copy from middle row left)
        if (u == fft_size/2 && v >= 1 && v < fft_size/2)
            u_conj = u;
            v_conj = fft_size-v;
        end
        % fill cell upper right (copy from lower cell left)
        if (u >= fft_size/2+1 && v >= 1 && v < fft_size/2)
            u_conj = fft_size-u;
            v_conj = fft_size-v;
        end
        % fill cell lower right (copy from upper cell left)
        if (u >= 1 && u < fft_size/2 && v >= 1 && v < fft_size/2)
            u_conj = fft_size-u;
            v_conj = fft_size-v;
        end
        
        % add coef to model and update residual
        if (u_conj ~= -1 && v_conj ~= -1)
            expansion_coefficient = orthogonality_correction * Rw(u+1, v+1) / W(1);
            G(u+1, v+1) = G(u+1, v+1) + fft_size^2 * expansion_coefficient;
            G(u_conj+1, v_conj+1) = conj(G(u+1, v+1));
            Rw = Rw -  expansion_coefficient * W_overal(fft_size-u+1:2*fft_size-u, fft_size-v+1:fft_size-v+1+fft_size/2) ...
                    -  conj(expansion_coefficient) * W_overal(fft_size-u_conj+1:2*fft_size-u_conj, fft_size-v_conj+1:fft_size-v_conj+1+fft_size/2);
            iter_counter = iter_counter + 1; % ... as two basis functions were added
        else
            expansion_coefficient = orthogonality_correction * Rw(u+1, v+1) / W(1);
            G(u+1, v+1) = G(u+1, v+1) + fft_size^2 * expansion_coefficient;
            Rw = Rw -  expansion_coefficient * W_overal(fft_size-u+1:2*fft_size-u, fft_size-v+1:fft_size-v+1+fft_size/2);
        end
        
        iter_counter = iter_counter + 1;
    end

    % get pixels from model
    g = ifft2(G);

    % extract reconstructed pixels
    extrapolated_block = real(g(fft_y_offset+(1:M), fft_x_offset+(1:N)));
    orig_samples = find(block_error_mask~=0);
    extrapolated_block(orig_samples) = non_reg_distorted_block(orig_samples);
end
