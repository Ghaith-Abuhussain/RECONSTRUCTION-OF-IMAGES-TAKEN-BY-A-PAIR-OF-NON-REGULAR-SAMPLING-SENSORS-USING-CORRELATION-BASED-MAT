function [processed_img] = processing_whole_image(non_reg_sampled_img, non_reg_sampling_mask, fse_parameters)

        % required parameters
        block_size = fse_parameters.block_size;
        block_size_min = fse_parameters.block_size_min;
        conc_weighting = fse_parameters.conc_weighting;
        fft_size = fse_parameters.fft_size;
        rho = fse_parameters.rhos(1);
        non_reg_sampled_img = double(non_reg_sampled_img);
        processed_img = double(non_reg_sampled_img);
        non_reg_sampling_mask =  double(max(0,sign(non_reg_sampling_mask)));
        border_width = floor(fft_size-block_size)/2;
        img_height = size(non_reg_sampled_img,1);
        img_width  = size(non_reg_sampled_img,2);
        
        % choosing the block size and process
        while (block_size >= block_size_min)
            
            num_block_in_column = ceil(size(non_reg_sampled_img,1)/block_size);
            num_block_in_line = ceil(size(non_reg_sampled_img,2)/block_size);
            % do actual processing per block
           for YBLOCK_INDEX = 0:num_block_in_column-1
                for XBLOCK_INDEX = 0:num_block_in_line-1

                    % calculation of the extrapolation area's borders
                    L_B = min(XBLOCK_INDEX*block_size, border_width);
                    T_B = min(YBLOCK_INDEX*block_size, border_width);
                    R_B = max(0, min(img_width-(XBLOCK_INDEX+1)*block_size, border_width));
                    B_B = max(0, min(img_height-(YBLOCK_INDEX+1)*block_size, border_width));

                    % extract blocks from images
                    EXTRACTED_BLOCK = processed_img((YBLOCK_INDEX*block_size-T_B+1):min(img_height,(YBLOCK_INDEX*block_size+block_size+B_B)), (XBLOCK_INDEX*block_size-L_B+1):min(img_width,(XBLOCK_INDEX*block_size+block_size+R_B)));
                    ERR_MASK = non_reg_sampling_mask((YBLOCK_INDEX*block_size-T_B+1):min(img_height,(YBLOCK_INDEX*block_size+block_size+B_B)), (XBLOCK_INDEX*block_size-L_B+1):min(img_width,(XBLOCK_INDEX*block_size+block_size+R_B)));
                    % get actual stddev value as it is needed to estimate the
                    % best number of iterations
                    S_D = standard_deviation(EXTRACTED_BLOCK, ERR_MASK);
          
                    % actual extrapolation
                    if(prod(prod(ERR_MASK))==0)
                        extrapolated_block_2d = processing_on_the_block_level(EXTRACTED_BLOCK, ERR_MASK, fse_parameters, rho, S_D);
                    else
                        extrapolated_block_2d=EXTRACTED_BLOCK;
                    end
                    % update image and mask
                    processed_img((YBLOCK_INDEX*block_size+1):min(img_height,(YBLOCK_INDEX+1)*block_size), (XBLOCK_INDEX*block_size+1):min(img_width,(XBLOCK_INDEX+1)*block_size)) = extrapolated_block_2d(T_B+1:end-B_B, L_B+1:end-R_B);
                    non_reg_sampling_mask((YBLOCK_INDEX*block_size+1):min(img_height,(YBLOCK_INDEX+1)*block_size), (XBLOCK_INDEX*block_size+1):min(img_width,(XBLOCK_INDEX+1)*block_size)) = ERR_MASK(T_B+1:end-B_B, L_B+1:end-R_B) + (1-sign(ERR_MASK(T_B+1:end-B_B, L_B+1:end-R_B)))*conc_weighting;

                    
                end
           end            
                       
            % set parameters for next extrapolation tasks (higher texture)
            block_size = block_size/2;
            border_width = (fft_size-block_size)/2;
            if block_size == 8
                rho = fse_parameters.rhos(2);
            end
            if block_size == 4
                rho = fse_parameters.rhos(3);
            end
            if block_size == 2
                rho = fse_parameters.rhos(4);
            end        
        end




end