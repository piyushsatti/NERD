clc;
clear all;
close all;

img = imread('lena_gray_512.tif');
% noise_img = imread('.7lena512.tif')
noise_img = imnoise(img,'salt & pepper',0.8);
noise_density(noise_img)

[out_img, s] = omf(noise_img);
out_img(1,1) = TLCF(out_img,s);
out_img(1,1:length(out_img)) = TRF(out_img);
out_img(1:length(out_img),1) = LCF(out_img);
out_img = RA(out_img);

noise_density(out_img)
imshow(uint8(out_img))

ssim(img,out_img)
psnr(img,out_img)

function [out1, out2] = omf(noise_img)

for i=1:length(noise_img)
    for j=1:length(noise_img)
        if noise_img(i,j)==255 ||noise_img(i,j)==0
            h(i,j) = 0;
        else
            h(i,j) = 1;
        end
    end
end
temp_img = noise_img.*uint8(h);
h=~h;
d = noise_density(temp_img);

med_filt_img=[];
w=1;
while d>0.45
    
    w=w+2;
    d_prev = noise_density(temp_img);
    med_filt_img(:,:,1+(w-3)/2) = medfilt2(noise_img,[w w]);
    for i=1:length(med_filt_img)
        for j=1:length(med_filt_img)
            if med_filt_img(i,j,1+(w-3)/2)==255 ||med_filt_img(i,j,1+(w-3)/2)==0
                m(i,j,1+(w-3)/2) = 0;
            else
                m(i,j,1+(w-3)/2) = 1;
            end
        end
    end
    temp_img = temp_img + uint8(med_filt_img(:,:,1+(w-3)/2)).*uint8(bitand(m(:,:,1+(w-3)/2),h));
    h=bitand(~m(:,:,1+(w-3)/2),h);
    
    d = noise_density(temp_img);
    if (d_prev>0.45)&&(d>0.45)
        out_img=temp_img;
        s=w;
    end
end

out1 = out_img;
out2=s;

end

function out = TLCF(temp_img,s)

mid=[];
if temp_img(1,1)==0
    for i=1:floor(s/2)+2
        for j=1:floor(s/2)+2
            if ~(temp_img(i,j)==0)
                mid = [mid temp_img(i,j)];
            end
        end
    end
    out=median(mid);
else
    out=temp_img(1,1);
end
end

function out = TRF(temp_img)

for i=2:length(temp_img)-1
    if temp_img(1,i)==0
        mid = nonzeros([temp_img(1,i-1) temp_img(2,i-1:i+1) temp_img(1,i+1)]);
        if ~isempty(mid)
            temp_img(1,i) = median(mid);
        else
            temp_img(1,i) = temp_img(1,length(temp_img)-1);
        end
    end
end
mid = nonzeros([temp_img(1,i-1) temp_img(2,i-1:i)]);
if ~isempty(mid)
    temp_img(1,length(temp_img)) = median(mid);
else
    temp_img(1,length(temp_img)) = temp_img(1,length(temp_img)-1);
end
out=temp_img(1,:);

end

function out = LCF(temp_img)

for i=2:length(temp_img)-1
    if temp_img(i,1)==0
        mid = nonzeros([temp_img(i-1,1); temp_img(i-1:i+1,2); temp_img(i+1,1)]);
        if ~isempty(mid)
            temp_img(i,1) = median(mid);
        else
            temp_img(i,1) = temp_img(length(temp_img)-1,1);
        end
    end
end
mid = nonzeros([temp_img(i-1,1) temp_img(i-1:i),2]);
if ~isempty(mid)
    temp_img(length(temp_img),1) = median(mid);
else
    temp_img(1,length(temp_img)) = temp_img(1,length(temp_img)-1);
end
out=temp_img(1,:);

end

function out = RA(t_img)

temp_img=double(t_img);

for i=2:length(temp_img)
    for j=2:length(temp_img)
        if (t_img(i,j)==0 && j~=length(temp_img))
            temp_img(i,j) = (temp_img(i-1,j-1)+temp_img(i-1,j)+temp_img(i,j-1)+temp_img(i-1,j+1))/4;
        elseif (t_img(i,j)==0 && j==length(temp_img))
            temp_img(i,j) = (temp_img(i-1,j-1)+temp_img(i-1,j)+temp_img(i,j-1))/3;
        end
    end
end

temp_img(length(temp_img),length(temp_img)) = (temp_img(i-1,j-1)+temp_img(i-1,j)+temp_img(i,j-1))/3;
out = uint8(temp_img);
end

function out = noise_density(out_img)

count=0;
for i=1:length(out_img)
    for j=1:length(out_img)
        if out_img(i,j)==255||out_img(i,j)==0
            count=count+1;
        end
    end
end

out=count/(length(out_img)^2);

end