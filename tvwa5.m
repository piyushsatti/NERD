% 2016 PRL 
% Removal of salt-and-pepper noise in corrupted image using three-values-weighted approach with variable-size window

function OutImg = tvwa5(nImg)
% 
% clc; clear all;
% Img = imread('lena_gray_512.tif'); % Reading input image
% d = 0.95; % Noise density
% nImg = imnoise(Img, 'salt & pepper', d); % Introducing noise
% nImg= imread('nImg.tif');
count1=0; count2=0; count3=0;count4=0;count5=0;count6=0;count7=0;
% nImg=imread('0.9lina_256.tif');
pad = 3; mfw = 3;

[row,col] = size(nImg); % Size calculation
imgZP = zeros(row + 2 * pad, col + 2 * pad); % Zero padding
OutImg = nImg;
k = 0; p = 0; % variable

imgZP(pad + 1:row + pad, pad + 1:col + pad) = nImg; %Zero padded image

rng = (mfw - 1) / 2;
flg = zeros(row, col);
mat5=[];
pixel= 127;
ar_max= size(row,col);
ar_min= size(row,col);
for i = 1 + pad: row + pad
    for j = 1 + pad:col + pad
        flgReg=1;
        if (imgZP(i,j) ~= 0)&&(imgZP(i,j) ~= 255)
            pixel = imgZP(i, j);
        else
            while(flgReg)
                Wn = imgZP(i-flgReg*rng:i+flgReg*rng, j-flgReg*rng:j+flgReg*rng);
                Wnf = Wn; Wnf(Wnf == 0) = []; Wnf(Wnf == 255) = [];
                if (length(Wnf) ~= 0)
                    max_val = max(Wnf); min_val = min(Wnf);
                    Gmax = []; Gmin = [];
                    for value = 1:length(Wnf)
                        a = abs(max_val - Wnf(value));
                        b = abs(min_val - Wnf(value));
                        if a <= b
                            Gmax = [Gmax, Wnf(value)];
                        else
                            Gmin = [Gmin, Wnf(value)];
                        end
                    end
                    pmax = length(Gmax)/length(Wnf);
                    pmin = length(Gmin)/length(Wnf);
                    mid_val = (pmax * max_val) + (pmin * min_val);
%                     if(pmax==0)
%                         mid_val= median(Gmin);
%                     elseif(pmin==0)
%                         mid_val= median(Gmax);
%                     else
%                         mid_val= (pmax*median(Gmax)+pmin*median(Gmin));
%                     end
                    if(pmax==0)
                        mid_val= median(Gmin);
                    elseif(pmin==0)
                        mid_val= median(Gmax);
                    else
                        mid_val= (pmax*median(Gmax)+pmin*median(Gmin));
                    end
%                     mid_val= (pmax*median(Gmax)+pmin*median(Gmin));
                    pixel= mid_val; 
                    break;
                else
                    flgReg= flgReg+1;
                    if(flgReg > 3)
                        count7=count7+1;
%                         i, j
%                         if(i>0&&i<row)&&(j>0&&j<col)
%                             pixel = uint8((double(OutImg(i,j-1))+double(OutImg(i-1,j-1))+double(OutImg(i-1,j))+double(OutImg(i-1,j+1)))/4);
%                             count6=count6+1;
%                         else
%                             pixel= mean(reshape(Wn(),1,[])); 
%                         end
                        break;
                    end
                end
            end
        end
    OutImg(i-pad, j-pad) = uint8(pixel);
    end
end

% 
% error = Img-OutImg;
% figure(1); imshow(Img);
% figure(2); imshow(nImg);
% figure(3); imshow(OutImg);
% SSIM = ssim(Img, OutImg);
% PSNR = psnr(Img, OutImg);

end