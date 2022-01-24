% 2016 PRL
% Adaptive switching weighted median filter framework for suppressing salt-and-pepper noise
% 
% clc; clear all;
 function OutImg = TVWA2final(nImg)
count1=0; count2=0; count3=0;count4=0;count5=0;count6=0;count7=0;
% Img = imread('lena_gray_512.tif'); % Reading input image
% d = 0.97; % Noise density
%  nImg = imnoise(Img, 'salt & pepper', d);
% Introducing noise
%  nImg=imread('.99lena512.tif');

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
            if (length(tmp) >0)
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
            if (length(tmp) > 2)
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
                    if (length(tmp)>1)
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
                       pixel=0;
% %                         imrovement in pSNR by 1 db;
%                         if (noise_density <= .85)
%                         count6=count6+1;
%                         pixel = (double(OutImg1(i - 2, j - 3)) + double(OutImg1(i - 3, j - 2)) + double(OutImg1(i - 3, j - 3)) + double(OutImg1(i - 3, j - 1))) / 4;
%                         else
%                         count6=count6+1;
%                         arr = [ OutImg1(i - 2, j - 2:j + 2) OutImg1(i - 1, j - 2:j + 2) OutImg1(i,j-2) OutImg1(i,j-1)];
%                         arr = reshape(arr(), 1, []);
%                         pixel = mean(reshape(arr(), 1, []));
%                     end
               
                    end  
                end
            end
            end
        
    OutImg(i - pad, j - pad) = uint8(pixel);
    OutImg1(i , j ) = uint8(pixel);
        
end
end


imgZP(1+pad:row+pad,1+pad:col+pad)=OutImg;
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
noise_density=count7./(row*col)
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
            if (length(tmp) > 4)
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
            if (length(tmp) > 6)
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
                    if (length(tmp)>1)
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
                       
% %                         imrovement in pSNR by 1 db;
                        if (noise_density <= .85)
                        count6=count6+1;
                        pixel = (double(OutImg1(i - 2, j - 3)) + double(OutImg1(i - 3, j - 2)) + double(OutImg1(i - 3, j - 3)) + double(OutImg1(i - 3, j - 1))) / 4;
                        else
                        count6=count6+1;
                        arr = [ OutImg1(i - 2, j - 2:j + 2) OutImg1(i - 1, j - 2:j + 2) OutImg1(i,j-2) OutImg1(i,j-1)];
                        arr = reshape(arr(), 1, []);
                        pixel = mean(reshape(arr(), 1, []));
                    end
               
                    end  
                end
            end
            end
        
    OutImg(i - pad, j - pad) = uint8(pixel);
    OutImg1(i , j ) = uint8(pixel);
        
end
end


% 
% % % 
% % % % 
% % % error = Img - OutImg;
% % % figure(1); imshow(Img);
% % % figure(2); imshow(nImg);
% % %  figure(3); imshow(OutImg);
% % % 
% % % mse = sum(sum((Img - OutImg) .^ 2)) / (row * col);
% % % PSNR = 10 * log10(255 ^ 2 / mse)
% % SSIM = ssim(Img, OutImg);
% PSNR = psnr(Img, OutImg)
 end
% %  
