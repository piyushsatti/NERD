clc;
clear all;
img= imread('lena_gray_256.tif');
d=.95;
nimg=imnoise(img,'salt & pepper',d);
% nimg= imread('0.9lina_256.tif');
d=nimg;
rng=1;
siz=length(nimg);
i=0;
% y=[0 1 2 3 4 5 6 7 8];
for i= 1:siz-2
    for j= 1: siz-2
        if (nimg(i+1,j+1)==0||nimg(i+1,j+1)==255)
%                     if(nimg(i,j)==0||nimg(i,j)==255)
%                         i=i+1;
              tmp= nimg(i:i+2,j:j+2);
              tmp(tmp==0)=[];
              tmp(tmp==255)=[];
            if(length(tmp)>1)
%                 x=(reshape(tmp(1:length(tmp),1:length(tmp)),1,[]));
                x=tmp;
                y=1:length(tmp);
                n= median(y);
                nimg(i+1,j+1)= spline(double(y),double(x),n);
            else
                nimg(i+1,j+1)= mean(reshape(nimg(i:i+2,j:j+2),1,[]));
            end
        end
    end
end
% 
figure
imshow(nimg);
% figure 
% imshow(d);                      
% error = img-nimg;                       
psnr=psnr(img,nimg)                       
   ssim=ssim(img,nimg)