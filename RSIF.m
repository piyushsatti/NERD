% 2014 SIViP 
% Recursive cubic spline interpolation filter approach for the removal of high density salt-and-pepper noise
function OutImg1 = RSIF(nImg)

% clc;clear;
% Img= imread('lena_gray_256.tif');   % Reading input image
% d= 0.7;   % Noise density
% nImg= imnoise(Img,'salt & pepper',d);   % Introducing noise

pad= 2; mfw=3; 
[row col]= size(nImg);              % Size calculation
imgZP= zeros(row+2*pad,col+2*pad);  % Zero padding
imgZP(pad+1:row+pad,pad+1:col+pad)= nImg; %Zero padded image

rng=(mfw-1)/2;
flg= zeros(row,col);
OutImg=imgZP;
for i= 1+pad:row+pad
    for j= 1+pad: col+pad
        if (imgZP(i,j)==0)||(imgZP(i,j)==255)
              tmp= OutImg(i-rng:i+rng,j-rng:j+rng);
              tmp(tmp==0)=[];
              tmp(tmp==255)=[];
            if(length(tmp)>1)
                x=tmp;
                y=1:length(tmp);
                n= median(y);
                OutImg(i,j)= spline(double(y),double(x),n);
            else
                OutImg(i,j)= mean(reshape(OutImg(i-rng:i+rng,j-rng:j+rng),1,[]));
            end
        end
    end
end

%   
OutImg1= uint8(OutImg(3:row+2,3:col+2));
% figure(1); imshow(Img);
% figure(2); imshow(uint8(OutImg1));                      
% error = Img-OutImg1;                       
% psnr=psnr(Img,OutImg1)                       
% ssim=ssim(Img,OutImg1)

end
