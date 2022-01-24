% 2016 PRL 
% Removal of salt-and-pepper noise in corrupted image using three-values-weighted approach with variable-size window

% function OutImg = tvwa2(nImg)
% 
clc; clear all;
Img = imread('lena_gray_256.tif'); % Reading input image
d = 0.9; % Noise density
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

for i = 1 + pad: row + pad
    for j = 1 + pad:col + pad
        flgReg=1;
        if (imgZP(i,j) ~= 0)&&(imgZP(i,j) ~= 255)
            pixel = imgZP(i, j);
        else
            while(flgReg)
                tmp1 = imgZP(i-flgReg*rng:i+flgReg*rng, j-flgReg*rng:j+flgReg*rng);
                tmp = tmp1; tmp(tmp == 0) = []; tmp(tmp == 255) = [];
                if (length(tmp) ~= 0)
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
%                     middle_value = (pmax * maximum) + (pmin * minimum);
                    if(length(mat1)==0)
                        middle_value= median(mat2);
                    elseif(length(mat2)==0)
                        middle_value= median(mat2);
                    else
                        middle_value= (pmax*median(mat1)+pmin*median(mat2));
                    end
%                     middle_value = pmax*median(mat1)+pmin*median(mat2);
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
%                     pixel = (pmax*maximum) + (pmin*minimum) + (pmid*middle_value);
                    if((length(mat1)==0)&&(length(mat2)==0))
                      pixel = pmid*median(mat3);
                    elseif((length(mat2)==0)&&(length(mat3)==0))
                        pixel = pmax*median(mat1);
                    elseif((length(mat3)==0)&&(length(mat1)==0))
                        pixel = pmin*median(mat2);
                    elseif(length(mat1)==0)
                        pixel = (pmin*median(mat2)) + (pmid*median(mat3));
                    elseif(length(mat2)==0)
                        pixel = (pmax*median(mat1)) + (pmid*median(mat3));
                    elseif(length(mat3)==0)
                        pixel = (pmax*median(mat1))+(pmin*median(mat2)) ;
                    else
                        pixel = (pmax*median(mat1)+pmin*median(mat2)+pmid*median(mat3));
                    end
%                     pixel = pmax*median(mat1)+pmin*median(mat2)+pmid*median(mat3);
                    break;
                else
                    flgReg= flgReg+1;
                    if(flgReg > 3)
                        break;
                    end
                end
            end
        end
    OutImg(i-pad, j-pad) = uint8(pixel);    
    end
end
edge= zeros(row, col);

% Sharpening
OutImgI= double(OutImg);
for i= 2: row-1
    for j= 2: col-1
        ed= abs(OutImgI(i-1,j)-OutImgI(i+1,j));
        if(ed>50)
            edge(i,j)= 255;
        end         
    end
end
incr= 20;
for i= 2: row-1
    for j= 2: col-4
        if(edge(i,j)== 255)&&(edge(i,j+1)== 255)&&(edge(i,j+2)== 255)&&(edge(i,j+3)== 255)
          OutImgI(i,j)=OutImgI(i,j)-incr;
          OutImgI(i,j+3)=OutImgI(i,j+3)+incr;
        elseif(edge(i,j)== 255)&&(edge(i,j+1)== 255)&&(edge(i,j+2)== 255)
          OutImgI(i,j)=OutImgI(i,j)-incr;
          OutImgI(i,j+2)=OutImgI(i,j+2)+incr;
        elseif(edge(i,j)== 255)&&(edge(i,j+1)== 255)
          OutImgI(i,j)=OutImgI(i,j)-incr;
          OutImgI(i,j+1)=OutImgI(i,j+1)+incr;
        end         
    end
end

OutImgII= double(OutImgI);
edgeII= zeros(row, col);
for i= 2: row-1
    for j= 2: col-1
        ed= abs(OutImgII(i,j-1)-OutImgII(i,j+1));
        if(ed>50)
            edgeII(i,j)= 255;
        end         
    end
end

for i= 2: row-4
    for j= 2: col-1
        if(edgeII(i,j)== 255)&&(edgeII(i+1,j)== 255)&&(edgeII(i+2,j)== 255)&&(edgeII(i+3,j)== 255)
          OutImgII(i,j)=OutImgII(i,j)-incr;
          OutImgII(i+3,j)=OutImgII(i+3,j)+incr;
        elseif(edgeII(i,j)== 255)&&(edgeII(i+1,j)== 255)&&(edgeII(i+2,j)== 255)
          OutImgII(i,j)=OutImgII(i,j)-incr;
          OutImgII(i+2,j)=OutImgII(i+2,j)+incr;
        elseif(edgeII(i,j)== 255)&&(edgeII(i+1,j)== 255)
          OutImgII(i,j)=OutImgII(i,j)-incr;
          OutImgII(i+1,j)=OutImgII(i+1,j)+incr;
        end         
    end
end

OutImgI= uint8(OutImgI);
OutImgII= uint8(OutImgII);
error = Img-OutImg;
figure(1); imshow(Img);
figure(2); imshow(nImg);
figure(3); imshow(OutImg);
figure(4); imshow(OutImgI);
figure(5); imshow(OutImgII);
SSIM = ssim(Img, OutImg);
SSIM1 = ssim(Img, OutImgI);
SSIM2 = ssim(Img, OutImgII);
PSNR = psnr(Img, OutImg);
PSNR1 = psnr(Img, OutImgI);
PSNR2 = psnr(Img, OutImgII);
% end