% 2018 Turkish 
% A new method based on pixel density in salt and pepper noise removal

function OutImg = BPDM2(nImg)
% 
% clc;clear;
% Img= imread('lena_gray_256.tif');   % Reading input image
% psnr1= zeros(1,9);
% ssim1= zeros(1,9);
% d= 9;   % Noise density
% % nImg= imnoise(Img,'salt & pepper',0.1*d);   % Introducing noise
% nImg= imread('nImg.tif');

pad= 3; mfw=3; 
[row col]= size(nImg);              % Size calculation
imgZP= zeros(row+2*pad,col+2*pad);  % Zero padding
imgZP(pad+1:row+pad,pad+1:col+pad)= nImg; %Zero padded image

rng=(mfw-1)/2;
flg= zeros(row,col);
ThMin=5; ThMax=250;

for i= pad+1: row+pad
    for j= pad+1: col+pad
        if((imgZP(i,j)==0)||(imgZP(i,j)==255))
            flg(i-pad, j-pad)= 1;
        else
            flg(i-pad, j-pad)= 0;
        end      
    end
end

count1=0;count2=0;count3=0;count4=0;

for i= 1+pad: row+pad
    for j= 1+pad: col+pad
        if((imgZP(i,j)==0)||(imgZP(i,j)==255))
            tmp= imgZP(i-rng:i+rng,j-rng:j+rng); % Start algo with 3x3
            tmp_nfp=tmp((tmp>ThMin)&(tmp<ThMax));
            tmp_np =tmp((tmp<ThMin+1)|(tmp>ThMax-1));
            if(length(tmp_nfp)>0)
                nfp= mode(tmp_nfp);
                tmp_nfp1= tmp_nfp(tmp_nfp~=nfp);
                nfp1= mode(tmp_nfp1);
                tmp_nfp2= tmp_nfp1(tmp_nfp1~=nfp1);
                nfp2= mode(tmp_nfp2);
                if(~isnan(nfp))&&(~isnan(nfp1))&&(~isnan(nfp2))
%                     pixel = double((double(nfp)+double(nfp1)+double(nfp2))/3);
                    pixel = median([double(nfp),double(nfp1),double(nfp2)]);
                    count2= count2+1;
                elseif(~isnan(nfp))&&(~isnan(nfp1))
                    pixel = double((double(nfp)+double(nfp1))/2);
                else
                    pixel = double(nfp);
                end
            else
                tmp= imgZP(i-2*rng:i+2*rng,j-2*rng:j+2*rng); % Start algo with 5x5
                tmp_nfp=tmp((tmp>ThMin)&(tmp<ThMax));
                tmp_np =tmp((tmp<ThMin+1)|(tmp>ThMax-1));
%                 if(~isempty(tmp_nfp))&&(~isempty(tmp_np)) 
                if(length(tmp_nfp)>0)
                    nfp= mode(tmp_nfp);
                    tmp_nfp1= tmp_nfp(tmp_nfp~=nfp);
                    nfp1= mode(tmp_nfp1);
                    tmp_nfp2= tmp_nfp1(tmp_nfp1~=nfp1);
                    nfp2= mode(tmp_nfp2);
                    if(~isnan(nfp))&&(~isnan(nfp1))&&(~isnan(nfp2))
%                         pixel = double((double(nfp)+double(nfp1)+double(nfp2))/3);
                        pixel = median([double(nfp),double(nfp1),double(nfp2)]);
                        count2= count2+1;
                    elseif(~isnan(nfp))&&(~isnan(nfp1))
                        pixel = double((double(nfp)+double(nfp1))/2);
                    else
                        pixel = double(nfp);
                    end
                else
                    tmp= imgZP(i-3*rng:i+3*rng,j-3*rng:j+3*rng); % Start algo with 7x7
                    tmp_nfp=tmp((tmp>ThMin)&(tmp<ThMax));
                    tmp_np =tmp((tmp<ThMin+1)|(tmp>ThMax-1));
%                     if(~isempty(tmp_nfp))&&(~isempty(tmp_np)) 
                    if(length(tmp_nfp)>0)
                        nfp= mode(tmp_nfp);
                        tmp_nfp1= tmp_nfp(tmp_nfp~=nfp);
                        nfp1= mode(tmp_nfp1);
                        tmp_nfp2= tmp_nfp1(tmp_nfp1~=nfp1);
                        nfp2= mode(tmp_nfp2);
                        if(~isnan(nfp))&&(~isnan(nfp1))&&(~isnan(nfp2))
%                             pixel = double((double(nfp)+double(nfp1)+double(nfp2))/3);
                            pixel = median([double(nfp),double(nfp1),double(nfp2)]);
                            count2= count2+1;
                        elseif(~isnan(nfp))&&(~isnan(nfp1))
                            pixel = double((double(nfp)+double(nfp1))/2);
                        else
                            pixel = double(nfp);
                        end
                    else
                        if(((i==1+pad)&&(j==1+pad))||((i==row+pad)&&(j==col+pad))) % first pixel mean of current window
                            tmp= imgZP(i-pad:i+pad,j-pad:j+pad);
                            pixel= mean(reshape(tmp,1,[]));
                        else
                            if(((i==1+pad)||(i==row+pad))&&(j~=1+pad))
                                pixel= OutImg(i-pad,j-pad-1); % previous processed signal in same row 
                            else
                                if(((j==1+pad)||(j==col+pad)))
                                    pixel= OutImg(i-pad-1,j-pad); % previous processed signal in same col
                                else % not belongs to boundary are estimated by previously processed four pixels
                                    pixel= (double(OutImg(i-3,j-4))+double(OutImg(i-4,j-3))+double(OutImg(i-4,j-4))+double(OutImg(i-4,j-2)))/4;
                                end
                            end
                        end
                    end
                end
            end
        else
            pixel= imgZP(i,j);
        end
        OutImg(i-pad,j-pad)= uint8(pixel);
    end
end

% figure(1); imshow(nImg);
% figure(2); imshow(OutImg);                      
% error = Img-nImg;  
% 
% psnr1=psnr(Img,OutImg);                      
% ssim=ssim(Img,OutImg);

end