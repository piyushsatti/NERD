
% function OutImg = damf(nImg)
clc;
clear all;
Img = imread('lena_gray_512.tif'); % Reading input image
% % %Img=[162,162,162,163,165,162,155;162,162,162,163,165,162,155;162,162,162,163,165,162,155;160,163,160,159,159,156,155;155,158,159,157,163,158,159;156,156,156,155,158,157,159;158,157,157,159,160,158,156];
d= 0.7;   % Noise density
%  nImg= imnoise(Img,'salt & pepper',d);   % Introducing noise
nImg = imread('.7lena512.tif');
% % nImg= [255,0,0,0,0,255,155;0,255,0,0,255,255,0;255,162,255,255,0,255,0;0,0,0,159,0,0,0;155,0,159,0,0,255,255;0,156,255,255,255,0,159;255,157,0,159,0,158,0];
%nImg=[0,255,255,163,165,162,155;0,255,255,163,165,162,155;162,162,162,163,165,162,155;160,163,160,159,159,255,155;155,158,159,157,163,158,159;156,156,156,155,158,157,159;158,157,157,159,160,158,156]
%         255 0 255 0 119
%         0 120 255 255 118
%         255 255 255 0 0
%         255 0 0 0 255];
pixel=nImg;
pad = 3;
mfw = 3;
[row col] = size(nImg); % Size calculation

k = 0; p = 0; % variable
pix = zeros(mfw * mfw, 1); % Pixels without noise
pixN = zeros(mfw * mfw, 1); % Noisy pixel
% imgZP(pad + 1:row + pad, pad + 1:col + pad) = nImg; % Zero padded image
imgZP = padarray(nImg,[4 4],'symmetric');


rng = (mfw - 1) / 2;
flg = zeros(row, col);
count = 0;
count1 = 0;
for i = 1: row+2*pad 
    for j = 1:col+2*pad
        if ((imgZP(i, j) == 0) || (imgZP(i, j) == 255))
            flg(i, j) = 0;
            count = count + 1;
        else
            flg(i, j) = 1;
        end
    end
end
for i = 1+pad: row + pad
    for j = pad + 1: col + pad
           if (flg(i,j)==0)
                tmp = imgZP(i - rng:i + rng, j - rng:j + rng);
                tmp(tmp == 0) = [];
                tmp(tmp == 255) = [];
                if(length(tmp)>0)
                pixel(i-pad,j-pad) = median(tmp);
                elseif(flg(i,j)==0||flg(i,j)==255)
                tmp = imgZP(i - 2:i + 2, j - 2:j + 2);
                tmp(tmp == 0) = [];
                tmp(tmp == 255) = [];
                  if(length(tmp)>0)
                  pixel(i-pad,j-pad) = median(tmp);
                elseif (flg(i,j)==0||flg(i,j)==255)
                tmp = imgZP(i - 3:i + 3, j - 3:j + 3);
                tmp(tmp == 0) = [];
                tmp(tmp == 255) = [];
                if(length(tmp)>0)
                pixel(i-pad,j-pad) = median(tmp);
                end
                  end
                end
           end       
    end
end
% imshow(pixel);
 
for i= 1: row+2
    for j= 1: col+2
        if((i==1)&&(j==1))
            imgZP(i,j)= pixel(i,j);
        elseif((i==1)&&(j==col+2))
            imgZP(i,j)= pixel(i,j-2);
        elseif((i==row+2)&&(j==1))
            imgZP(i,j)= pixel(i-2,j);
        elseif((i==row+2)&&(j==col+2))
            imgZP(i,j)= pixel(i-2,j-2);
        elseif(i==1)
            imgZP(i,j)= pixel(i,j-1);
        elseif(i==row+2)
            imgZP(i,j)= pixel(i-2,j-1);
        elseif(j==1)
            imgZP(i,j)= pixel(i-1,j);
        elseif(j==col+2)
            imgZP(i,j)= pixel(i-1,j-2);
        else
            imgZP(i,j)= pixel(i-1,j-1);
        end
    end
end   
for i = 1: row+2
    for j = 1:col+2
        if ((imgZP(i, j) == 0) || (imgZP(i, j) == 255))
            flg(i, j) = 0;
            count = count + 1;
        else
            flg(i, j) = 1;
        end
    end
end
OutImg=pixel;
for i = 2: row+1
    for j = 2: col+1
           if (flg(i,j)==0)
                tmp = imgZP(i - 1:i + 1, j - 1:j + 1);
                tmp(tmp == 0) = [];
                tmp(tmp == 255) = [];
                if(length(tmp)>0)
                OutImg(i-1,j-1) = median(tmp);
                end
                  end
                end
end       
% imshow(OutImg);
PSNR = psnr(Img, OutImg)
% end           