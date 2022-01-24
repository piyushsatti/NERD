function OutImg= mdbutm(nImg);
% clc
% clear all % Noise density
% Img= imread('lena_gray_256.tif');
% nImg= imread('0.2lina_256.tif');
OutImg= nImg;

[row col]= size(nImg);       % Size calculation
imgZP= zeros(row+4,col+4);  % Zero padding

pix= zeros(9,1);    % Pixels without noise
pixN= zeros(9,1);   % Noisy pixel

imgZP(3:row+2,3:col+2)= nImg; %Zero padded image

mfw=3; k=0; p=0;% variable 
rng=(mfw-1)/2;
flg= zeros(row,col);
% figure
% imshow(J);
%[row col]= size();
for x=3:row+2
    for y=3:col+2
       i=0;
        if (imgZP(x,y)==0||imgZP(x,y)==255)
          b= reshape(imgZP(x-1:x+1,y-1:y+1),1,[]);
          tmp=b;
          tmp(tmp==0)=[];
          tmp(tmp==255)=[];
            %b=[ J(x,y) J(x,y+1) J(x,y+2) J(x+1,y) J(x+1,y+1) J(x+1,y+2) J(x+2,y) J(x+2,y+1) J(x+2,y+2) ];
%           for x1=1:9
%               if(b(1,x1)==0||b(1,x1)==225)   
          if(sum(tmp)==0)
           OutImg(x-2,y-2)=uint8(mean(b)); 
           
          else
            b=sort(tmp);
            OutImg(x-2,y-2)= uint8(median(b));
        
         end
        end
    end
end
% 
% figure
% imshow(OutImg);
% title('filtered');
% figure 
% imshow(I);
% title('orignal');
%  
%  mse= sum(sum((Img-OutImg).^2))/(row*col);
%  PSNR= 10*log(255^2/mse)
% PSNR=psnr(OutImg,Img)
% SSIM= ssim(OutImg,Img)
end