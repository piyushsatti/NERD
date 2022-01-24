% palinisamy 2016
clc;
clear all;
count1=0;count2=0;
Img= imread('lena_gray_256.tif');   % Reading input image
% d= 0.9;   % Noise density
%  nImg= imnoise(Img,'salt & pepper',d);   % Introducing noise
 nImg=imread('0.1lina_256.tif');
pad= 2;
mfw=3; a
% nImg= [ 255 255 118 0 0
%         255 0 255 0 121
%         0 120 255 255 118
%         255 255 255 0 0
%         255 0 0 0 255];


[row col]= size(nImg);       % Size calculation
imgZP= zeros(row+2*pad,col+2*pad);  % Zero padding
OutImg1=imgZP;
k=0; p=0;% variable
pix= zeros(mfw*mfw,1);    % Pixels without noise
pixN= zeros(mfw*mfw,1);   % Noisy pixel

imgZP(pad+1:row+pad,pad+1:col+pad)= nImg; %Zero padded image

rng=(mfw-1)/2;
flg= zeros(row,col);
tol=0;            %tolerence of the pixel%

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
for i= pad+1: row+pad
    for j= pad+1: col+pad
        tmp=imgZP(i-rng:i+rng,j-rng:j+rng);
        tmp(tmp==0)=[];
        tmp(tmp==255)=[];
        if length(tmp)==0
            tmp=imgZP(i-2*rng:i+2*rng,j-2*rng:j+2*rng);
            tmp(tmp==0)=[];
            tmp(tmp==255)=[];
            if length(tmp)==0
                pixel=mean(reshape(imgZP(i-2*rng:i+2*rng,j-2*rng:j+2*rng),1,[]));
            else
                    informationPixel=tmp(1);
                    infor1=0;infor2=0;
                    tmp1=[];tmp2=[];
                    ar=25;
                    for (a=1:length(tmp()))
                        if ( informationPixel -ar >   tmp(a)) && (informationPixel + ar <   tmp(a))
                               infor1= infor1+1;
                               tmp1(infor1+1)=tmp(a);
                        else
                            infor2=infor2+1;
                            tmp2(infor2+1)=tmp(a);
                        end
                    end
                    if  (infor1 >= infor2)
                        pixel=mean(reshape(tmp1(),1,[]));
                    else
                        pixel=mean(reshape(tmp2(),1,[]));
                    end
            end
        else
                pixel=mean(reshape(tmp(),1,[]));
        end
        
        
        if (abs(pixel- imgZP(i,j))>=tol)
              n_pixel=pixel;
        else
            n_pixel=imgZP(i,j)
        end
        OutImg(i-pad,j-pad)=uint8(n_pixel);
    end
end

wavelet_thersholding=OutImg;
sX=size(i)
[ca,chd,cvd,cdd]=swt2(wavelet_thersholding,2,'db1');
A1 = wcodemat(ca(:,:,1),255);
H1 = wcodemat(chd(:,:,1),255);
V1 = wcodemat(cvd(:,:,1),255);
D1 = wcodemat(cdd(:,:,1),255);

A2 = wcodemat(ca(:,:,2),255);
H2 = wcodemat(chd(:,:,2),255);
V2 = wcodemat(cvd(:,:,2),255);
D2 = wcodemat(cdd(:,:,2),255);

thr=.4
ythard1 = wthresh(A1,'h',thr);
ytsoft1 = wthresh(A1,'s',thr);

ythard2 = wthresh(H1,'h',thr);
ytsoft2 = wthresh(H1,'s',thr);

ythard3 = wthresh(V1,'h',thr);
ytsoft3 = wthresh(V1,'s',thr);

ythard4 = wthresh(D1,'h',thr);
ytsoft4 = wthresh(D1,'s',thr);


rec1 = uint8(iswt2(ythard1,ythard2,ythard3,ythard4,'db1'));
rec2 = uint8(iswt2(ytsoft1,ytsoft2,ytsoft3,ytsoft4,'db1'));

error = Img - OutImg;
figure(1); imshow(Img);
figure(2); imshow(nImg);
figure(3); imshow(rec1);

mse = sum(sum((Img - rec1) .^ 2)) / (row * col);
PSNRformula = 10 * log10(255 ^ 2 / mse)
PSNR = psnr(Img, rec1, 255)
SSIM = ssim(Img, rec1);

