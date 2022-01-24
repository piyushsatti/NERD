function OutImg= ASWMF(nImg)
% clc; clear;
%  Img = imread('lena_gray_256.tif'); % Reading input image
% 
% Img=imread('Lena_gray_512.tif');
% % Img= uint8([159,156,155,162,157;159,158,155,153,158;156,156,156,162,157;157,156,155,155,152;158,155,153,155,156]);
% % Img=uint8([166,169,169,165,166;166,167,166,166,163;167,170,166,165,167;163,168,165,166,170;168,164,164,162,163]);
% % Img=[127 128 
% %     155  200]
% d = 0.9; % Noise density
%  nImg = imnoise(Img, 'salt & pepper', d); % Introducing noise
% nImg=uint8( [166,0,255,0,255;0,255,0,0,0;0,255,0,255,0;255,0,255,255,0;168,0,0,0,255]);
% nImg=imread('0.9lina_256.tif');
% ND=.9;
% nImg=uint8([255,255,255,0,255;255,0,0,255,255;156,255,0,255,255;0,255,255,155,255;255,0,0,155,0]);
%  nImg= uint8([0,0,255,0,0;0,0,0,0,0;0,0,0,0,255;255,0,0,0,0;157,255,0,0,159]);
pad = 4;
mfw = 3;
count = 0;
count1 = 0;
count = 0;
count2 = 0; count3 = 0;




% nImg= [ 255 255 118 0 0
%         255 0 255 0 119
%         0 120 255 255 118
%         255 255 255 0 0
%         255 0 0 0 255];

pixel=0;

[row col] = size(nImg); % Size calculation
imgZP = zeros(row + 2 * pad, col + 2 * pad); % Zero padding
OutImg = nImg;
k = 0; p = 0; % variable
pix = zeros(mfw * mfw, 1); % Pixels without noise
pixN = zeros(mfw * mfw, 1); % Noisy pixel

imgZP(pad + 1:row + pad, pad + 1:col + pad) = nImg; % Zero padded image

rng = (mfw - 1) / 2;
flg = zeros(row, col);

for i = pad + 1: row + pad
    for j = pad + 1: col + pad
        tmp1 = imgZP(i - rng:i + rng, j - rng:j + rng);
        tmp = mean(reshape(tmp1(), 1, []));
        if (((imgZP(i, j) ~= 0) && (imgZP(i, j) ~= 255)) )
           % if imgZP(i,j)~=0 || imgZP(i,j)~=255 
            flg(i - pad, j - pad) = 1;
        else
            flg(i - pad, j - pad) = 0;
        end
    end
end

for i = 1 + pad: row + pad
    for j = 1 + pad:col + pad
        if (flg(i - pad, j - pad) == 1)
            pixel = imgZP(i, j);
         
        else
         
            tmp1 = imgZP(i - rng:i + rng, j - rng:j + rng);
            tmp1(tmp1 == 0) = [];
            tmp1(tmp1 == 255) = [];
            
            
            if (length(tmp1) ~= 0)
                length(tmp1)  ;  
                count3=count3+1;
                tmp1 = reshape(imgZP(i - rng:i + rng, j - rng:j + rng), 1, []);
                tmp1(5) = [];
                tmp = [tmp1, tmp1];
                tmp(tmp == 0) = [];
                tmp(tmp == 255) = [];
                if mod(length(tmp()), 2) == 1
                    pixel = median(reshape(tmp(), 1, []));
                    pp=pixel;
                 
                else
                    C = unique(tmp());
                    if (length(tmp()) ~= length(C()))
                        a = mode(tmp());
                        tmp = [tmp, a];
                        pixel = median(reshape(tmp(), 1, []));
                        pp=pixel;
                    else
%                         error is comming in this step
                        tmp1 = reshape(imgZP(i - rng:i + rng, j - rng:j + rng), 1, []);
                        tmp1(5) = [];
                        tmp = [tmp1, tmp1];
                        tmp(tmp == 0) = [];
                        tmp(tmp == 255) = [];
                        tmp;
                        pixel = median(reshape(tmp(), 1, []));
                        pp=pixel;
                    end
                end
            else
                tmp1 = imgZP(i - 2 * rng:i + 2 * rng, j - 2 * rng:j + 2 * rng);
                tmp1(tmp1 == 0) = [];
                tmp1(tmp1 == 255) = [];
                tmp1;
                if (length(tmp1()) ~= 0)
                    count=count+1;
                    tmp1 = reshape(imgZP(i - 2 * rng:i + 2 * rng, j - 2 * rng:j + 2 * rng), 1, []);
                    b = [imgZP(i - 1, j), imgZP(i + 1, j), imgZP(i, j + 1), imgZP(i, j - 1)];
                    tmp1 = [tmp1, b, b, b, b];
                    c = [imgZP(i - 2, j), imgZP(i, j - 2), imgZP(i, j + 2), imgZP(i + 2, j)];
                    tmp1 = [tmp1, c, c, c];
                    d = [imgZP(i - 2, j - 1), imgZP(i - 1, j - 1), imgZP(i + 1, j - 1), imgZP(i + 2, j - 1), imgZP(i - 2, j + 1), imgZP(i - 1, j + 1), imgZP(i + 1, j + 1), imgZP(i + 2, j + 1)];
                    tmp1 = [tmp1, d, d];
%                     e = [imgZP(i - 2, j - 2), imgZP(i - 2, j + 1), imgZP(i - 2, j + 2), imgZP(i -1, j - 2), imgZP(i - 1, j + 2), imgZP(i + 1, j + 2), imgZP(i + 2, j -2), imgZP(i + 2, j + 2)];
%                     tmp1 = [tmp1, e];
%                     
                    tmp1(tmp1 == 0) = [];
                    tmp1(tmp1 == 255) = [];
                    tmp = tmp1;
                    if mod(length(tmp()), 2) == 1
                        pixel = median(reshape(tmp(), 1, []));
                        pp=pixel;
                    else
                        C = unique(tmp());
                        if (length(tmp()) ~= length(C()))
                            a = mode(tmp());
                            tmp = [tmp, a];
                            pixel = median(reshape(tmp(), 1, []));
                            pp=pixel;
                        else
                            tmp1 = reshape(imgZP(i - 2*rng:i + 2*rng, j - 2*rng:j + 2*rng), 1, []);
                            tmp1(5) = [];
                            tmp = [tmp1, tmp1];
                            tmp(tmp == 0) = [];
                            tmp(tmp == 255) = [];
                            tmp;
                            i;
                            j;
                            pixel = median(reshape(tmp(), 1, []));
                            pp=pixel;
                        end
                    end
                
            
                else
                  
                    count1 = count1 + 1;
                    tmp1 = imgZP(i - 3 * rng, j - 3 * rng:j +3 * rng);
                    tmp2 = imgZP(i - 3 * rng:i + 3 * rng, j - 3 * rng);
                    tmp3 = imgZP(i + 3 * rng, j - 3 * rng:j + 3 * rng);
                    tmp4 = imgZP(i - 3 * rng:i + 3 * rng, j + 3 * rng);
                    pixel1 = 0;
                    if pixel1 == 0
                        for value = 1:length(tmp1())
                            if (tmp1(value) ~= 0 && tmp1(value) ~= 255)
                                pixel1 = tmp1(value);
                            end
                        end
                    end
                    if pixel1 == 0
                        for value = 1:length(tmp2())
                            if (tmp2(value) ~= 0 && tmp2(value) ~= 255)
                                pixel1 = tmp2(value);
                            end
                        end
                    end
                    if pixel1 == 0
                        for value = 1:length(tmp3())
                            if (tmp3(value) ~= 0 && tmp3(value) ~= 255)
                                pixel1 = tmp3(value);
                            end
                        end
                    end
                    if pixel1 == 0
                        for value = 1:length(tmp4())
                            if (tmp4(value) ~= 0 && tmp4(value) ~= 255)
                                pixel1 = tmp4(value);
                            end
                        end
                        
%                     end
                    new=[];
                    if pixel1==0
                       count2=count2+1;
                       tmp1 = imgZP(i - 2 * rng:i + 2 * rng, j - 2 * rng:j + 2 * rng);
                       pixel1=mean(reshape(tmp1(),1,[]));
%                  
                    end
                    pixel=pixel1;

%                       if (noise_density <= .85)
%                         count1=count1+1;
%                         pixel = (double(OutImg(i - 2, j - 3)) + double(OutImg(i - 3, j - 2)) + double(OutImg(i - 3, j - 3)) + double(OutImg(i - 3, j - 1))) / 4;
%                     else
%                         arr = [OutImg1(i - 2, j - 2:j + 2) OutImg1(i - 1, j - 2:j + 2) OutImg1(i,j-2) OutImg1(i,j-1)];
%                         arr = reshape(arr(), 1, []);
%                         pixel = mean(reshape(arr(), 1, []));


%                     end


                end
                end
            end
        end
     
        OutImg(i - pad, j - pad) = pixel;
        pp=pixel;
    end
end


% 
% error = Img - OutImg;
% figure(1); imshow(Img);
% figure(2); imshow(nImg);
% figure(3); imshow(OutImg);
% 
% mse = sum(sum((Img - OutImg) .^ 2)) / (row * col);
% PSNRformula = 10 * log10(255 ^ 2 / mse)
% PSNR = psnr(Img, OutImg, 255)
% SSIM = ssim(Img, OutImg);
% count;

end
% 
% 
% 
% 
% 
% 
