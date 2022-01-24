% 2016 PRL
% Adaptive switching weighted median filter framework for suppressing salt-and-pepper noise
% 
% clc; clear all;
 function OutImg = PAtern(nImg)
count1=0; count2=0; count3=0;count4=0;count5=0;count6=0;count7=0;
% Img = imread('lena_gray_512.tif'); % Reading input image
% d = 0.9; % Noise density
%  nImg = imnoise(Img, 'salt & pepper', d); % Introducing noise
%  nImg=imread('.7lena512.tif');

% % nImg= [ 255 255 118 0 0
% %         255 0 255 0 119
% %         0 120 255 255 118
% %         255 255 255 0 0
% %         255 0 0 0 255];
pad = 3; mfw = 3;
[row,col] = size(nImg); % Size calculation
imgZP = zeros(row + 3 * pad, col +3* pad); % Zero padding
OutImg = nImg;
k = 0; p = 0; % variable

imgZP(pad + 1:row + pad, pad + 1:col + pad) = nImg; %Zero padded image

rng = (mfw - 1) / 2;
flg = zeros(row, col);
OutImg1 = imgZP;
for i = pad + 1: row + pad
    for j = pad + 1: col + pad
        if ((imgZP(i, j) == 0) || (imgZP(i, j) == 255))
            flg(i - pad, j - pad) = 1;
            count7=count7+1;
        else
            flg(i - pad, j - pad) = 0;
        end
    end
end
noise_density=count7./(row*col);
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
            if (length(tmp) > 0)
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
                    tmp1 = imgZP(i - 3*rng:i + 3*rng, j - 3*rng:j + 3*rng);
                    tmp = tmp1;
                    tmp(tmp == 0)=[];  tmp(tmp == 255)=[];
                    if (length(tmp)>0)
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
%                        pixel=0;
% %                         imrovement in pSNR by 1 db;
%                         if (noise_density <= .85)
%                         count6=count6+1;
%                         pixel = (double(OutImg1(i - 2, j - 3)) + double(OutImg1(i - 3, j - 2)) + double(OutImg1(i - 3, j - 3)) + double(OutImg1(i - 3, j - 1))) / 4;
%                         else
%                         count6=count6+1;
                        arr = [ OutImg1(i - 2, j - 2:j + 2) OutImg1(i - 1, j - 2:j + 2) OutImg1(i,j-2) OutImg1(i,j-1)];
                        arr = reshape(arr(), 1, []);
                        pixel = mean(reshape(arr(), 1, []));
%                     end
               
                    end  
                end
            end
            end
        
    OutImg(i - pad, j - pad) = uint8(pixel);
    OutImg1(i , j ) = uint8(pixel);
        
end
end



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
                else
                    OutImg2(i - pad, j - pad) = sum(sum(tmp .* Kgauss3)) ./ 16;
                end
            end
        end

OutImg2 = uint8(OutImg2);
OutImg = OutImg2;





% 
% 
% % 
% error = Img - OutImg;
% figure(1); imshow(Img);
% figure(2); imshow(nImg);
%  figure(3); imshow(OutImg);
% 
% mse = sum(sum((Img - OutImg) .^ 2)) / (row * col);
% PSNR = 10 * log10(255 ^ 2 / mse)
% SSIM = ssim(Img, OutImg);
% PSNR = psnr(Img, OutImg)
 end
%  
