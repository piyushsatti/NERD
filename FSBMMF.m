% 2014 AEU 
% Fast switching based median–mean filter for high density saltand pepper noise removal SAP noise removal
function OutImg = fsbmmf(nImg)
% clc;clear;
% % 
%  Img= imread('lena_gray_512.tif');   % Reading input image
%  d= 0.1;   % Noise density
%  nImg=imnoise(Img,'salt & pepper',d);


% nImg= imread('0.7lina_256.tif');   % Introducing noise

% nImg= [ 255 255 118 0 0
%         255 0 255 0 119
%         0 120 255 255 118
%         255 255 255 0 0
%         255 0 0 0 255];

[row col]= size(nImg);       % Size calculation
imgZP= zeros(row+4,col+4);  % Zero padding

pix= zeros(9,1);    % Pixels without noise
pixN= zeros(9,1);   % Noisy pixel

imgZP(3:row+2,3:col+2)= nImg; %Zero padded image

mfw=3; k=0; p=0;% variable 
rng=(mfw-1)/2;
flg= zeros(row,col);

for i= 3: row+2
    for j= 3: col+2
        if((imgZP(i,j)==0)||(imgZP(i,j)==255))
            flg(i-2, j-2)= 1;
        else
            flg(i-2, j-2)= 0;
        end      
    end
end

for i= 3: row+2
    for j= 3: col+2
        if(flg(i-2,j-2)==0)
            pixel= imgZP(i,j);
        else
            pixel= median(reshape(imgZP(i-rng:i+rng,j-rng:j+rng),1,[]));
%                     tmp= imgZP(i-rng:i+rng,j-rng:j+rng);
%                     tmp(tmp==0)=[];
%                     tmp(tmp==255)=[];
%                     pixel= median(tmp);
            if((pixel==0)||(pixel==255))%write code to compute noisy pixel
                if(((i==3)&&(j==3))||((i==3)&&(j==col+2))||((i==row+2)&&(j==3))||((i==row+2)&&(j==col+2)))
                    tmp= imgZP(i-2*rng:i+2*rng,j-2*rng:j+2*rng);
                    tmp(tmp==0)=[];
                    tmp(tmp==255)=[];
                    pixel= median(tmp);
                elseif((i==3)||(i==row+2))
                    pixel= OutImg(i-2,j-3); % previous processed signal 
                elseif((j==3)||(j==col+2))
                    pixel= OutImg(i-3,j-2); % previous processed signal
                else
                    pixel= (double(OutImg(i-2,j-3))+double(OutImg(i-3,j-2))+double(OutImg(i-3,j-3))+double(OutImg(i-3,j-1)))/4;
                end
            end
        end
        OutImg(i-2,j-2)= uint8(pixel);
    end
end
%     PSNR=psnr(Img,OutImg)



 end