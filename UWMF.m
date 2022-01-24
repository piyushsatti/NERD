
function OutImg = UWMF(nImg)
% clc; clear all; 
% Img = imread('lena_gray_512.tif');
% d= 0.2;  
% nImg= imnoise(Img,'salt & pepper',d);
% nImg = imread('.7lena512.tif');
 pad = 5; % according to window size (mfw+1)/2;
% mfw =9;
[row,col] = size(nImg);
imgZP = padarray(nImg,[pad pad],'symmetric');
imgZP = double(imgZP);
count =0 ;
for i=1+pad:row+pad
    for j=1+pad:col+pad
        if (imgZP(i,j) ==0 || imgZP(i,j) ==255 )
        count=count+1;
        end
    end
end
noise_density = count./(row*col);


% choose according to noisde density 
% Recommended window sizes for the proposed method.
% Noise density (p) Window size
% p < 20                3 × 3
% 20 ? p < 50       5 × 5
% 50 ? p < 70       7 × 7
% 70 ? p < 85       9 × 9
% 85 ? p < 90       11 × 11
% p ? 90                13 × 13


if noise_density < .18
    mfw =3;
elseif (noise_density>=.18 && noise_density<.48)
    mfw=5;
elseif (noise_density>=.48 && noise_density<.68)
   mfw =7;
elseif (noise_density>=.68 && noise_density<.83)
    mfw=9;
elseif (noise_density>=.83 && noise_density<.88)
    mfw=11;
else
    mfw=13;
end

% mfw = 3
pad =(mfw+1)/2;
imgZP = padarray(nImg,[pad pad],'symmetric');
imgZP = double(imgZP);

rng = (mfw - 1) / 2;

k=6;

% tmp =[192 105 255 81 78;208 163 255 89 255;2524 205 255 255 255;225 0 255 255 89;227 228 222 198 133;];
num =1000000;
count1 =0;
count2=0;
count3=0;
count4=0;
p=1;
for i = 1+pad:row+pad
    for j = 1+pad:col+pad
        if ((imgZP(i,j) == 0) || (imgZP(i,j) == 255))
             count4=count4+1;
            sumi=0;sumwdash=0;
            count255=0;count0=0;
            P=0;Q=0;R=0;S=0;T=0;
            tmp = imgZP(i-rng:i+rng,j-rng:j+rng);
            for x=-rng+i:rng+i
                for y = -rng+j:rng+j
                    if  (imgZP(x,y) ~=0) && (imgZP(x,y) ~=255)
%                         check = imgZP(x,y)  
                        dxyx = x-i;
                        dxyy = y-j;
                        wxy(x-i+pad,y-j+pad) = power(power((power(abs(x-i),p)+power(abs(y-j),p)),1/p),-k);
                        P=P+((wxy(x-i+pad,y-j+pad)*dxyx)*dxyx);
                        Q=Q+((wxy(x-i+pad,y-j+pad)*dxyx)*dxyy);
                        R=R+((wxy(x-i+pad,y-j+pad)*dxyx));
                        S=S+((wxy(x-i+pad,y-j+pad)*dxyy)*dxyy);
                        T=T+((wxy(x-i+pad,y-j+pad)*dxyy));
                        
                    end
                end
            end
            R=-R;
            T=-T;
            gydash= ((P*T)-(Q*R))/((-power(Q,2))+(P*S));
            gxdash= ((R)-(Q*gydash))/P;
            for x=1:length(tmp)
                for y =1:length(tmp)
                    if tmp(x,y)==255
                        count255 =count255+1;
                    elseif tmp(x,y) ==0
                        coutn0=count0+1;
                    else
                        wxydash = wxy(x,y) +(wxy(x,y)*((gxdash*dxyx)+(gydash*dxyy)));
                        sumwdash=sumwdash+wxydash;
                        num = tmp(x,y)*wxydash;
                        sumi = sumi + (tmp(x,y)*wxydash);
                    end
                end
            end
            tmp = imgZP(i-rng:i+rng,j-rng:j+rng);
            if (count255+count0) == power(mfw,2)
                if count255>count0
                    imgZP(i,j)=255;
                    count1=count1+1;
                else
                    imgZP(i,j)=0;
                    count2=count2+1;
                end
            else
%                 imgZP(i,j) = (sumi)/(sumwdash);
                   pix =  (sumi)/(sumwdash);
                count3=count3+1;
            end
        else
%            imgZP(i,j) = imgZP(i,j);
                pix = imgZP(i,j);
%             count1=count1+1;

        end
        OutImg(i-pad,j-pad) = uint8(pix);
    end
end

% OutImg = uint8(imgZP(1+pad:row+pad,1+pad:col+pad));

% subplot(1,3,1)
% imshow(Img)
% subplot(1,3,2)
% imshow(nImg)
% subplot(1,3,3)
% imshow(OutImg)




% 
% count = 0;
% pad =2 ;
% for i = pad + 1: row + pad
%     for j = pad + 1: col + pad
%         if ((imgZP(i, j) == 0) || (imgZP(i, j) == 255))
%             flg(i - pad, j - pad) = 1;
%             count = count + 1;
%         else
%             flg(i - pad, j - pad) = 0;
%         end
%     end
% end
% noise_density = count ./ (row * col);
% 
% th = 0.6;
% OutImg2 = double(OutImg);
% imgZP = padarray(OutImg2,[2 2],'symmetric');
% edged3 = edge(OutImg, 'canny', .15);
% windowSize = 4;
% % kernel = ones(windowSize) / windowSize ^ 2;
% kernel = [1 2 2 1; 2 4 4 2; 2 4 4 2; 1 2 2 1] ./ 36;
% blurryImage = conv2(single(edged3), kernel, 'same');
% edged3 = blurryImage > .2; 
% Kgauss3 = [1 2 1; 2 4 2; 1 2 1];
% Kmean3 = [1 1 1; 1 1 1; 1 1 1];
% 
% if noise_density > th
%         for i = 1 + pad:row + pad
%             for j = 1 + pad:col + pad
%                 tmp = imgZP(i - 1:i + 1, j - 1:j + 1);
%                 if edged3(i - pad, j - pad) == 1
%                     OutImg2(i - pad, j - pad) = sum(sum(tmp .* Kmean3)) ./ 9;
%                 else
%                     OutImg2(i - pad, j - pad) = sum(sum(tmp .* Kgauss3)) ./ 16;
%                 end
%             end
%         end
% end
% 
% OutImg2 = uint8(OutImg2);

% subplot(1,2,1)
% imshow(OutImg)
% subplot(1,2,2)
% imshow(OutImg2)

% PSNR1 = psnr(Img,OutImg)
% PSNR2 =psnr(Img,OutImg2)


end



