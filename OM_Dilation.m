clc
clear all
close all

img = imread('lena_gray_512.tif');
% noise_img = imnoise(img,'salt & pepper',0.70);
noise_img = imread('.7lena512.tif')
[out_img, s] = omf(noise_img);

se = strel('line',11,90);
out_img = imerode(out_img,se);

for i=1:length(out_img)
    for j=1:length(out_img)
        if out_img(i,j)==255 ||out_img(i,j)==0
            out_img(i,j) = img(i,j);
        end
    end
end

% out_img(1,1) = TLCF(out_img,s);
% out_img(1,1:length(out_img)) = TRF(out_img);
% out_img(1:length(out_img),1) = LCF(out_img);
% out_img = RA(out_img);

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
t_img=noise_img;
med_filt_img=[];
w=3;
    
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

out1=temp_img;
out2=w;

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








OutImg2 = out_img;
[row col] = size(noise_img);
% imgZP = zeros(row + 3 * pad, col + 3 * pad); % Zero padding
 
    k = 0; p = 0; % variable
 Kgauss3 = uint8([1 2 1; 2 4 2; 1 2 1]);
    Kmean3 = uint8([1 1 1; 1 1 1; 1 1 1]);
imgZP = padarray(OutImg2,[2 2],'symmetric'); 
 pad =2;
 edged3 = edge(out_img, 'canny', .15);
    % edged3  =  edge(i,'canny',.12);
    % % figure(2)
    % % subplot(1,3,1);imshow(edged1);
    % % subplot(1,3,2);imshow(edged2);
    % % subplot(1,3,3);imshow(edged3);
    %
    % figure(3)
    % subplot(1,3,1);imshow(edged3);
    windowSize = 4;
    % kernel = ones(windowSize) / windowSize ^ 2;
    kernel = [1 2 2 1; 2 4 4 2; 2 4 4 2; 1 2 2 1] ./ 36;
    blurryImage = conv2(single(edged3), kernel, 'same');
    edged3 = blurryImage > .2; % Rethreshold
 
    % subplot(1,3,2);imshow(edged3);
    %
    count999 = 0;
 
        for i = 1 + pad:row + pad
            for j = 1 + pad:col + pad
                tmp = imgZP(i - 1:i + 1, j - 1:j + 1);
                count999 = count999 + 1;
                if edged3(i - pad, j - pad) == 1
                    OutImg2(i - pad, j - pad) = sum(sum(tmp .* Kmean3)) ./ 9;
                else
                    OutImg2(i - pad, j - pad) = sum(sum(tmp .* Kgauss3)) ./ 16;
                end
            end
        end
    
 
OutImg2 = uint8(OutImg2);
