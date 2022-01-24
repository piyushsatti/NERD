% CODE for the other one called morphology_mean_filter:

% 
% clc
% clear all
% close all

 function OutImg = morphology_mean_filter(nImg)

% Img = imread('lena_gray_512.tif');
% d=0.7;
% nImg = imnoise(Img,'salt & pepper',d);
% nImg = imread('.7lena512.tif');
pad = 1;
imgZP = double(padarray(nImg,[pad pad],'symmetric'));
mfw =3;
rng = (mfw - 1)/2;
[row,col] = size(imgZP);

for i=1:row
    for j=1:col
        if imgZP(i,j)==0||imgZP(i,j)==255
            b_f(i,j) = 0;
        else
            b_f(i,j) = 1;
        end
    end
end

for i=2:row-1
    for j=2:col-1
        if b_f(i,j) == 0
            temp = min(nonzeros(imgZP(i-1:i+1,j-1:j+1).*b_f(i-1:i+1,j-1:j+1)));
            if ~isempty(temp)
                imgZP(i,j) = temp;
            end
        end
    end
end

for i=1:row
    for j=1:col
        if imgZP(i,j)==0||imgZP(i,j)==255
            b_g1(i,j) = 0;
        else
            b_g1(i,j) = 1;
        end
    end
end

for i=2:row-1
    for j=2:col-1
        if b_g1(i,j) == 0
            temp = max(nonzeros(imgZP(i-1:i+1,j-1:j+1).*b_g1(i-1:i+1,j-1:j+1)));
            if ~isempty(temp)
                imgZP(i,j) = temp;
            end
        end
    end
end

for i=1:row
    for j=1:col
        if imgZP(i,j)==0||imgZP(i,j)==255
            b_g2(i,j) = 0;
        else
            b_g2(i,j) = 1;
        end
    end
end

for i=2:row-1
    for j=2:col-1
        if b_g2(i,j) == 0
            temp = max(nonzeros(imgZP(i-1:i+1,j-1:j+1).*b_g2(i-1:i+1,j-1:j+1)));
            if ~isempty(temp)
                imgZP(i,j) = temp;
            end
        end
    end
end

for i=1:row
    for j=1:col
        if imgZP(i,j)==0||imgZP(i,j)==255
            b_g3(i,j) = 0;
        else
            b_g3(i,j) = 1;
        end
    end
end

for i=2:row-1
    for j=2:col-1
        if b_g3(i,j) == 0
            temp = min(nonzeros(imgZP(i-1:i+1,j-1:j+1).*b_g3(i-1:i+1,j-1:j+1)));
            if ~isempty(temp)
                imgZP(i,j) = temp;
            end
        end
    end
end

for i=2:row-1
    for j=2:col-1
        if b_f(i,j)==0
            imgZP(i,j) = round(mean(nonzeros(imgZP(i-1:i+1,j-1:j+1)),'all'));
        end
    end
end

OutImg = uint8(imgZP(2:row-1,2:col-1));

% PSNR = psnr(Img,OutImg)
% SSIM = ssim(Img,OutImg);

% imshow(OutImg)

 end