function [sampled_image,error_mask] = non_regular_sampling (LR_image,sampling_ratio)
% TAKE THE SIZE OF THE LOW-RESOLUTION IMAGE
[M,N,K]=size(LR_image);
S_R = sqrt(sampling_ratio);% THE SIZE OF THE SAMPLING BLOCK
err_mask=zeros(S_R*M,S_R*N,3);
s1=zeros(S_R,S_R,3);%SAMPLED BLOCK FOR THE IMAGE
s2=zeros(S_R,S_R);%SAMPLED BLOCK FOR THE MASK
for i=1:M
    for j=1:N
        ord1 =randi([1 S_R]); ord2 =randi([1 S_R]);% PUT THE PIXEL RANDOMELY IN THE SAMPLED BLOCK
        s1(ord1,ord2,1:3)=LR_image(i,j,1:3);
        s2(ord1,ord2)=255;
        V1((i-1)*S_R+1:i*S_R,(j-1)*S_R+1:j*S_R,1:3)=s1;
        V2((i-1)*S_R+1:i*S_R,(j-1)*S_R+1:j*S_R)=s2;
        s1=zeros(S_R,S_R);
        s2=zeros(S_R,S_R);
    end
end
V1=uint8(V1);
err_mask(:,:,1)=V2;
err_mask(:,:,2)=V2;
err_mask(:,:,3)=V2;
err_mask=uint8(err_mask);
sampled_image=V1;
error_mask=err_mask;

end