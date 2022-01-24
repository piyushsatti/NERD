clc;
clear all;
 Img= imread('lena_gray_256.tif'); 
figure(1)
imshow(Img)

% nImg= imread('0.7lina_256.tif');
 nImg=imnoise(Img,'salt & pepper',.1);
 
%  nImg=[0 0 255 0 255 ;
%      122 123 255 255 255;
%      0 255 255 255 123;
%      0 122 255 124 255;
%      0 255 255 255 255]
figure(2)
imshow(nImg)
[row col]= size(nImg);  
rep=0;
for (i=1:row -2)
    for (j=1:col-2)
     
        new_mat=[nImg(i,j)  nImg(i,j+1)   nImg(i,j+2);
                 nImg(i+1,j)  nImg(i+1,j+1)   nImg(i+1,j+2);
                 nImg(i+2,j)  nImg(i+2,j+1)   nImg(i+2,j+2)];       
            a1=sort(new_mat(1,1:3));
            a2=sort(new_mat(2,1:3));%sorting row wise
            a3=sort(new_mat(3,1:3));
             new_mat=[a1 ; a2 ; a3];
            
            b1=sort(new_mat(1:3,1));
            b2=sort(new_mat(1:3,2));%sorting coloumn wise
            b3=sort(new_mat(1:3,3));
            new_mat=[b1 b2 b3];
            
            
            val3=new_mat(3,1);
            val2=new_mat(2,2);%sorting right Diagonal wise 
            val1=new_mat(1,3);
            if val1 > val2
                rep=val1;val1=val2;val2=rep;
            end
            if val2 > val3
                rep=val2;val2=val3;val3=rep;
            end
            if val1 > val3
                rep=val1;val1=val3;val3=rep;
            end
            
            
            new_mat(3,1)=val3;new_mat(2,2)=val2;new_mat(1,3)=val1;
            
            
            
%             Decision Based Median Filter
            if (new_mat(1,1)<nImg(i+1,j+1)) && (nImg(i+1,j+1)<new_mat(3,3))  && (0<new_mat(3,3)) &&(new_mat(3,3) <255) 
                nImg(i+1,j+1)=new_mat(2,2);
            else
                pixel=median(reshape(new_mat(),1,[]));
                
                if  (pixel>new_mat(1,1)) && (pixel>0) && (pixel<new_mat(3,3)) && (pixel<255)
                    nImg(i+1,j+1)=pixel;
                else
                    nImg(i+1,j+1)=nImg(i+1,j);
                end
               
              
            end  
      
      
      end
end

 PSNR=psnr(Img,nImg)
% Un-symmetric Trimmed Mean Filtering (UTMF)
for (i=1:row-2)
    for (j=1:col-2)
        new_mat=[nImg(i,j)  nImg(i,j+1)   nImg(i,j+2) nImg(i+1,j)  nImg(i+1,j+1)   nImg(i+1,j+2) nImg(i+2,j)  nImg(i+2,j+1)   nImg(i+2,j+2)];
        
        new_mat(new_mat==0)=[];
        new_mat(new_mat==255)=[];
        nImg(i+1,j+1)=mean(reshape(new_mat(),1,[]));
        
    end
    
end





error= Img-nImg;
figure(3)
imshow(nImg)
 
 mse= sum(sum((Img-nImg).^2))/(row*col);
    PSNRformaula= 10* log10((255.^2)./mse)
  PSNR=psnr(Img,nImg)
 
%  pSnR=(PSNR+PSNRformaula)/2
SSIM= ssim(Img, nImg)

