%Level-2 ITSAMFT 2010
% Improved Tolerance based Selective Arithmetic Mean Filter for Detection and Removal of Impulse Noise
%  function OutImg = l2ITSAMFT(nImg)
clc;clear;
% 
% 
Img= imread('lena_gray_256.tif');   % Reading input image
d= 0.5;   % Noise density
% nImg= imnoise(Img,'salt & pepper',d);   % Introducing noise
nImg=imread('0.1lina_256.tif');
% nImg= [ 255 255 118 0 0
%         255 0 255 0 119
%         0 120 255 255 118
%         255 255 255 0 0
%         255 0 0 0 255];
pad= 2;
mfw=3; 
[row col]= size(nImg);       % Size calculation
imgZP= zeros(row+2*pad,col+2*pad); % Zero padding

k=0; p=0;% variable
pix= zeros(mfw*mfw,1);    % Pixels without noise
pixN= zeros(mfw*mfw,1);   % Noisy pixel

imgZP(pad+1:row+pad,pad+1:col+pad)= nImg; %Zero padded image
size(imgZP);
rng=(mfw-1)/2;
flg= zeros(row,col);

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
tol=0;




for i= 1+pad: row+pad
    for j= 1+pad: col+pad
            rng=1;       
            tmp= imgZP(i-rng:i+rng,j-rng:j+rng);
            siz1 = size(tmp);
            siz1=siz1(1).*siz1(2);
            
            tmp(tmp==0)=[];tmp(tmp==255)=[];
            tmp=reshape(tmp(),1,[]);
            siz2 = size(tmp);
            siz2=siz2(2);
            if (siz2 > (1/3)*siz1)
                
                tmp=median(reshape(tmp(),1,[]));
                if (abs(tmp- imgZP(i,j))>=tol)
                    pixel=tmp;
                else
                    pixel=imgZP(i,j);
                end
            else
                rng=2;
                tmp= imgZP(i-rng:i+rng,j-rng:j+rng);
                tmp(tmp==0)=[];tmp(tmp==255)=[];
                tmp=reshape(tmp(),1,[]);
                
                if (length(tmp())==0)
                    pixel=mean(reshape(imgZP(i-rng:i+rng,j-rng:j+rng),1,[]));
                    continue
                else
                    informationPixel=tmp(1);
                    infor1=0;infor2=0;
                    tmp1=[];tmp2=[];
                    ar=100;
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
                        tmp=mean(reshape(tmp1(),1,[]));
                    else
                        tmp=mean(reshape(tmp2(),1,[]));
                    end
                    if (abs(tmp- imgZP(i,j))>=tol)
                    pixel=tmp;
                    else
                    pixel=imgZP(i,j);
                    end
                end
            end
           OutImg(i-pad,j-pad)=uint8(pixel); 
            
    end
    
end

% PSNRfunc=psnr(Img,OutImg)
% for i= 1+pad: row+pad
%   for i=1+pad:1+pad+2
%     for j= 1+pad: col+pad
%        
%          tmp= imgZP(i-rng:i+rng,j-rng:j+rng);
%          tmp(tmp==0)=[];
%         tmp(tmp==255)=[];
%          OutImg(i-pad,j-pad)= mean(reshape(tmp(),1,[]));
%        
%        
%     end
%     
% end
% end



error= Img-OutImg;
figure(1); imshow(Img);
figure(2); imshow(nImg);
figure(3); imshow(OutImg);

mse= sum(sum((Img-OutImg).^2))/(row*col);
PSNR= 10*log10(255^2/mse)
PSNRfunc=psnr(Img,OutImg,255)
SSIM= ssim(Img, OutImg)




%  end































