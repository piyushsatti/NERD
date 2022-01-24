% 2016 AEU 
% Probabilistic decision based filter to remove impulse noise usingpatch else trimmed median
% function OutImg = pdbm(nImg)
clc;clear;

Img= imread('lena_gray_512.tif');   % Reading input image
d= 0.1;   % Noise density
 nImg= imnoise(Img,'salt & pepper',d);   % Introducing noise
% nImg= imread('0.9lina_256.tif');

% nImg= [ 255 255 118 0 0
%         255 0 255 0 119
%         0 120 255 255 118
%         255 255 255 0 0
%         255 0 0 0 255];
pad= 2;
mfw=3; 
[row col]= size(nImg);       % Size calculation
imgZP= zeros(row+2*pad,col+2*pad);  % Zero padding

k=0; p=0;% variable
pix= zeros(mfw*mfw,1);    % Pixels without noise
pixN= zeros(mfw*mfw,1);   % Noisy pixel

imgZP(pad+1:row+pad,pad+1:col+pad)= nImg; %Zero padded image

rng=(mfw-1)/2;
flg= zeros(row,col);

for i= pad+1: row+pad
    for j= pad+1: col+pad
        if((imgZP(i,j)==0)||(imgZP(i,j)==255))
            flg(i-pad, j-pad)= 1;
        else
            flg(i-pad, j-pad)= 0;
        end      
    end
end

Max=3; pp=128;

for i= 1+pad: row+pad
    for j= 1+pad: col+pad
        if((imgZP(i,j)==0)||(imgZP(i,j)==255)) % Check pixel is noisy
            tmp= sort(sort(imgZP(i-rng:i+rng,j-rng:j+rng)));
            if((tmp(1+rng,1+rng)==0)||(tmp(1+rng,1+rng)==255)) % Check pixel after PM is noisy 
                tmp1= imgZP(i-rng:i+rng,j-rng:j+rng); % Compute mf using tm with 3x3 window
                tmp1(tmp1==0)=[];
                tmp1(tmp1==255)=[];
                pixel1= median(tmp1);
                if((pixel1==0)||(pixel1==255))    % Check pixel after tm is noisy 
                    tmp2= imgZP(i-2*rng:i+2*rng,j-2*rng:j+2*rng); % Compute mf using tm with 5x5 window
                    tmp2(tmp2==0)=[];
                    tmp2(tmp2==255)=[];
                    pixel2= median(tmp2);
                    if((pixel2==0)||(pixel2==255))      % Check noisy after tm
                        pixel=pp;                       % Compute mf using ppr
                    else
                        pixel= pixel2;                  % Assign mfed using tm with 5x5
                    end
                else
                    pixel= pixel1;                      % Assign mfed using tm with 3x3
                end
            else
                pixel= tmp(1+rng,1+rng);                % Assign mfed using PM with 3x3
            end
        else
            pixel= imgZP(i,j);                          % Assign non-noisy orignal pixel
        end
        OutImg(i-pad,j-pad)= uint8(pixel);
        pp= uint8(pixel);
    end
end
% 
error= Img-OutImg;
figure(1); imshow(Img);
figure(2); imshow(nImg);
figure(3); imshow(OutImg);

mse= sum(sum((Img-OutImg).^2))/(row*col);
PSNRformula= 10*log10(255^2/mse)
SSIM= ssim(Img, OutImg);
PSNR= psnr(Img,OutImg)

% end