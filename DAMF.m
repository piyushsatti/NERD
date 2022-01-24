clc; clear all; close all;

Img = imread('lena_gray_512.tif'); % Reading input image
d= 0.8;   % Noise density
nImg= imnoise(Img,'salt & pepper',d);   % Introducing noise

pixel=nImg;
pad = 3;
mfw = 3;
[row, col] = size(nImg); % Size calculation
imgZP = zeros(row + 2 * pad, col + 2 * pad); % Zero padding
OutImg1 = imgZP;
k = 0; p = 0; % variable
pix = zeros(mfw * mfw, 1); % Pixels without noise
pixN = zeros(mfw * mfw, 1); % Noisy pixel
imgZP(pad + 1:row + pad, pad + 1:col + pad) = nImg; % Zero padded image
    for i=1+pad:row +pad
        for j=1+pad:3 +pad
            if j==1+pad
                imgZP(i,pad)=imgZP(i,j);
            elseif j==2+pad
                imgZP(i,pad-1)=imgZP(i,j);
            else
                imgZP(i,pad-2)=imgZP(i,j);
            end
        end
    end
    
for i=1+pad:3 +pad
        for j=1+pad:col +pad
            if i==1+pad
                imgZP(pad,j)=imgZP(i,j);
            elseif i==2+pad
                imgZP(pad-1,j)=imgZP(i,j);
            else
                imgZP(pad-2,j)=imgZP(i,j);
            end
        end
end

for i=1+pad:row +pad
        for j=col+pad-2:col+pad
            if j==col+pad-2
                imgZP(i,col+pad+3)=imgZP(i,j);
            elseif j==col+pad-1
                imgZP(i,col+pad+2)=imgZP(i,j);
            else
                imgZP(i,col+pad+1)=imgZP(i,j);
            end
        end
end
    
   
for i=row+pad -2 : row +pad
        for j=1+pad:col+pad
            if i==row+pad-2
                imgZP(row+pad+3,j)=imgZP(i,j);
            elseif i==row+pad-1
                imgZP(row+pad+2,j)=imgZP(i,j);
            else
                imgZP(row+pad+1,j)=imgZP(i,j);
            end
        end
end 
imgZP(1,1)=imgZP(6,1);imgZP(1,2)=imgZP(6,2);imgZP(1,3)=imgZP(6,3);
imgZP(2,1)=imgZP(5,1);imgZP(2,2)=imgZP(5,2);imgZP(2,3)=imgZP(5,3);
imgZP(3,1)=imgZP(4,1);imgZP(3,2)=imgZP(4,2);imgZP(3,3)=imgZP(4,3);
imgZP(1,col+pad+1)=imgZP(1,col+pad);imgZP(1,col+pad+2)=imgZP(1,col+pad-1);imgZP(1,col+pad+3)=imgZP(1,col+pad-2);
imgZP(2,col+pad+1)=imgZP(2,col+pad);imgZP(2,col+pad+2)=imgZP(2,col+pad-1);imgZP(2,col+pad+3)=imgZP(2,col+pad-2);
imgZP(3,col+pad+1)=imgZP(3,col+pad);imgZP(3,col+pad+2)=imgZP(3,col+pad-1);imgZP(3,col+pad+3)=imgZP(3,col+pad-2);
imgZP(row+pad+1,1)=imgZP(row+pad,1);imgZP(row+pad+1,2)=imgZP(row+pad,2);imgZP(row+pad+1,3)=imgZP(row+pad,3);
imgZP(row+pad+2,1)=imgZP(row+pad-1,1);imgZP(row+pad+2,2)=imgZP(row+pad-1,2);imgZP(row+pad+2,3)=imgZP(row+pad-1,3);
imgZP(row+pad+3,1)=imgZP(row+pad-2,1);imgZP(row+pad+3,2)=imgZP(row+pad-2,2);imgZP(row+pad+3,3)=imgZP(row+pad-2,3);
% 
imgZP(row+pad+1,col + pad +1)=imgZP(row+pad,col + pad +1);imgZP(row+pad+1,col + pad +2)=imgZP(row+pad,col + pad +2);imgZP(row+pad+1,col + pad +3)=imgZP(row+pad,col + pad +3);
imgZP(row+pad+2,col + pad +1)=imgZP(row+pad-1,col + pad +1);imgZP(row+pad+2,col + pad +2)=imgZP(row+pad-1,col + pad +2);imgZP(row+pad+2,col + pad +3)=imgZP(row+pad-1,col + pad +3);
imgZP(row+pad+3,col + pad +1)=imgZP(row+pad-2,col + pad +1);imgZP(row+pad+3,col + pad +2)=imgZP(row+pad-2,col + pad +2);imgZP(row+pad+3,col + pad +3)=imgZP(row+pad-2,col + pad +3);

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
% PSNR = psnr(Img, OutImg)


    Kgauss3 = [1 2 1; 2 4 2; 1 2 1];
    Kmean3 = [1 1 1; 1 1 1; 1 1 1];

OutImg2 = OutImg;
%     imgZP = zeros(row + 3 * pad, col + 3 * pad); % Zero padding
 
    k = 0; p = 0; % variable
 
%     imgZP(pad + 1:row + pad, pad + 1:col + pad) = OutImg; %Zero padded image
    imgZP = padarray(OutImg,[2 2],'symmetric');
 
    edged3 = edge(OutImg, 'canny', .15);
    windowSize = 4;
    kernel = [1 2 2 1; 2 4 4 2; 2 4 4 2; 1 2 2 1] ./ 36;
    blurryImage = conv2(single(edged3), kernel, 'same');
    edged3 = blurryImage > .2; 
 
    count999 = 0;
 
        for i = 1 + pad:row + pad
            for j = 1 + pad:col + pad
                tmp = double(imgZP(i - 1:i + 1, j - 1:j + 1));
                count999 = count999 + 1;
                if edged3(i - pad, j - pad) == 1
                    OutImg2(i - pad, j - pad) = sum(sum(tmp .* Kmean3)) ./ 9;
             
                end
            end
        end
OutImg2 = uint8(OutImg2);
OutImg = OutImg2;



end           