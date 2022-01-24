% 2016 PRL
% Removal of salt-and-pepper noise in corrupted image using three-values-weighted approach with variable-size window

% function OutImg = tvwa1(nImg)
% 
clc; clear all;
Img = imread('lena_gray_512.tif'); % Reading input image
d = 0.95; % Noise density
nImg = imnoise(Img, 'salt & pepper', d); % Introducing noise
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
for i = pad + 1: row + pad
    for j = pad + 1: col + pad
        if ((imgZP(i, j) == 0) || (imgZP(i, j) == 255))
            flg(i - pad, j - pad) = 1;
        else
            flg(i - pad, j - pad) = 0;
        end
    end
end

for i = 1 + pad: row + pad
    for j = 1 + pad:col + pad
        if (flg(i - pad, j - pad) == 0)
            pixel = imgZP(i, j);
            count1=count1+1;
        else
            tmp1 = imgZP(i - rng:i + rng, j - rng:j + rng);
            tmp = tmp1;
            tmp(tmp == 0) = [];
            tmp(tmp == 255) = [];
            if (length(tmp) ~= 0)
               count2=count2+1;
               maximum = max(tmp); minimum = min(tmp);
               mat1 = []; mat2 = [];
               for value = 1:length(tmp)
                    a = abs(maximum - tmp(value));
                    b = abs(minimum - tmp(value));
                    if a <= b
                        mat1 = [mat1, tmp(value)];
                    else
                        mat2 = [mat2, tmp(value)];
                    end
               end
               pmax = length(mat1) / length(tmp);
               pmin = length(mat2) / length(tmp);
               middle_value = (pmax * maximum) + (pmin * minimum);
               mat1 = []; mat2 = []; mat3 = [];
               for value = 1:length(tmp)
                    a = abs(maximum - tmp(value));
                    b = abs(minimum - tmp(value));
                    c = abs(middle_value - tmp(value));
                    if (a<b && a<c) %Max
                        mat1 = [mat1, tmp(value)];
                    else
                        if (b<a && b<c) % Min
                            mat2 = [mat2, tmp(value)];
                        else % Mid
                            mat3 = [mat3, tmp(value)];
                        end
                    end
               end
               pmax = length(mat1)/length(tmp);
               pmin = length(mat2)/length(tmp);
               pmid = length(mat3)/length(tmp);
               pixel = (pmax * maximum) + (pmin * minimum) + (pmid * middle_value);
            else
                tmp1 = imgZP(i - 2 * rng:i + 2 * rng, j - 2 * rng:j + 2 * rng);
                tmp = tmp1;
                tmp(tmp == 0) = []; tmp(tmp == 255) = [];
                if (length(tmp) ~= 0)
                    count3=count3+1;
                    maximum = max(tmp); minimum = min(tmp);
                    mat1 = []; mat2 = [];
                    for value = 1:length(tmp)
                        a = abs(maximum - tmp(value));
                        b = abs(minimum - tmp(value));
                        if a <= b
                            mat1 = [mat1, tmp(value)];
                        else
                            mat2 = [mat2, tmp(value)];
                        end
                    end
                    pmax = length(mat1) / length(tmp);
                    pmin = length(mat2) / length(tmp);
                    middle_value = (pmax * maximum) + (pmin * minimum);
                    mat1 = []; mat2 = []; mat3 = [];
                    for value = 1:length(tmp)
                        a = abs(maximum - tmp(value));
                        b = abs(minimum - tmp(value));
                        c = abs(middle_value - tmp(value)); 
                        if (a < b) && (a < c)
                            mat1 = [mat1, tmp(value)];
                        elseif (b < a) && (b < c)
                            mat2 = [mat2, tmp(value)];
                        else
                            mat3 = [mat3,tmp(value)];    
                        end
                    end
                    pmax = length(mat1)/length(tmp);
                    pmin = length(mat2)/length(tmp);
                    pmid = length(mat3)/length(tmp);
                    pixel = (pmax * maximum) + (pmin * minimum) + (pmid * middle_value);
                else
                    tmp1 = imgZP(i - 3*rng:i + 3*rng, j - 3*rng:j + 3*rng);
                    tmp = tmp1;
                    tmp(tmp == 0) = [];  tmp(tmp == 255) = [];
                    if (length(tmp)~=0)
                        count4=count4+1;
                        maximum = max(tmp); minimum = min(tmp);
                        mat1 = []; mat2 = [];
                        for value = 1:length(tmp)
                            a = abs(maximum - tmp(value));
                            b = abs(minimum - tmp(value));
                            if a <= b
                                mat1 = [mat1,tmp(value)];
                            else
                                mat2 = [mat2, tmp(value)];
                            end
                        end
                        pmax = length(mat1) / length(tmp);
                        pmin = length(mat2) / length(tmp);
                        middle_value = (pmax * maximum) + (pmin * minimum);
                        mat1 = []; mat2 = []; mat3 = [];
                        for value = 1:length(tmp)
                            a = abs(maximum - tmp(value));
                            b = abs(minimum - tmp(value));
                            c = abs(middle_value - tmp(value));
                            if (a < b) && (a < c)
                                mat1 = [mat1, tmp(value)];
                            elseif (b < a) && (b < c)
                                mat2 = [mat2, tmp(value)];
                            else
                                mat3 = [mat3, tmp(value)];
                            end
                        end
                        pmax = length(mat1)/length(tmp);
                        pmin = length(mat2)/length(tmp);
                        pmid = length(mat3)/length(tmp);
                        pixel = (pmax*maximum) + (pmin*minimum) + (pmid* middle_value);
                    else
                        if(i>0&i<row)&(j>0&j<col)
                            pixel = uint8((double(OutImg(i,j-1))+double(OutImg(i-1,j-1))+double(OutImg(i-1,j))+double(OutImg(i-1,j+1)))/4);
                        else
                            tmp= imgZP(i-3*rng:i+3*rng,j-3*rng:j+3*rng);
                            pixel= mean(reshape(tmp(),1,[]));
                            mat5=[mat5, pixel];
                        end
                    end  
                end
            end
        end
    OutImg(i - pad, j - pad) = uint8(pixel);    
    end
end


error = Img - OutImg;
figure(1); imshow(Img);
figure(2); imshow(nImg);
figure(3); imshow(OutImg);

mse = sum(sum((Img - OutImg) .^ 2)) / (row * col);
PSNR = 10 * log(255 ^ 2 / mse);
SSIM = ssim(Img, OutImg);
PSNR = psnr(Img, OutImg)

% end