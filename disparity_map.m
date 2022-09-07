function y = disparity_map (IMAGE_LEFT,IMAGE_RIGHT,w_size,t,gamma)
I_left=im2double(rgb2gray(IMAGE_LEFT));%TRANSFER TO GRAY SCALE
I_right=im2double(rgb2gray(IMAGE_RIGHT));
[M,N]=size(I_left);
disp_map = zeros(M,N);
window_per_culumn=floor(M/w_size);% THE NUMER OF BLOCKS IN EACH CULUMN
for i=1:window_per_culumn
    for j=1:N-w_size % MOVE ON EACH PIXLE IN THE ROW OF THE LEFT_IMAGE
        bottom_border_left=min(M,i*w_size);
        Left_Window = I_left((i-1)*w_size+1:bottom_border_left,j:j+w_size-1);
        for k=1:N-w_size% MOVE ON EACH PIXLE IN THE ROW OF RIGHT_IMAGE 
            bottom_border_right=min(M,i*w_size);
            right_Window = I_right((i-1)*w_size+1:bottom_border_right,k:k+w_size-1);
            MSE(1,k)=sqrt(mean2((Left_Window-right_Window).^2));% THE MATCHING STANDARD
        end
        [~,minMSE_index]=min(MSE);
        if min(MSE)<=gamma % WE DETERMINE GAMMA IF IT BELOW OR ABOVE A LIMIT 
            if(t==1) % IF THE DISPARITY IS FROM LEFT TO RIGHT
                disp_map((i-1)*w_size+1:bottom_border_left,j:j+w_size-1)=j-minMSE_index;
            else % IF THE DISPARITY IS FROM RIGHT TO LEFT
                disp_map((i-1)*w_size+1:bottom_border_left,j:j+w_size-1)=minMSE_index-j;
            end
        end
    end
end
y=disp_map;
end