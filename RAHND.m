%%% Function to run all mf functions on benchmark images with varying noise
%%% density.

clc; clear all;                                 % Clear all existing variables
a=.1;                                          % Default value of Noise density
  Img= imread('lena_gray_512.tif');

for d= 1:8
    nImg = imnoise(Img,'salt & pepper',a);   % Introducing noise
 
    DBAMF_OutImg = DBAMF(nImg);                         % Call DBAMF median filter 2007 SPL
    tmp7(d,1) = psnr(Img,DBAMF_OutImg);                   % Calculate PSNR of DBAMF  
    eval(['psnr_DBAMF' num2str(d) '= tmp7(d,1)']);
    tmp8(d,1)=ssim(Img,DBAMF_OutImg);                   % Calculate SSIM
    eval(['ssim_DBAMF' num2str(d) '= tmp8(d,1)']);
    
    ASWMF_OutImg = ASWMF(nImg);                 % Call ASWMF median filter 2010 ICIIS
    tmp7(d,2)=psnr(Img,ASWMF_OutImg);               % Calculate PSNR of fsbmmf  
    eval(['psnr_ASWMF' num2str(d) '= tmp7(d,2)']);
    tmp8(d,2)=ssim(Img,ASWMF_OutImg);               % Calculate SSIM
    eval(['ssim_ASWMF' num2str(d) '= tmp8(d,2)']);
% %     
    mdbutm_OutImg = mdbutm(nImg);                       % Call mdbutm median filter 2011 SPL 
    tmp7(d,3)=psnr(Img,mdbutm_OutImg);                  % Calculate PSNR of mdbutm  
    eval(['psnr_mdbutm' num2str(d) '= tmp7(d,3)']);
    tmp8(d,3)=ssim(Img,mdbutm_OutImg);                  % Calculate SSIM
    eval(['ssim_mdbutm' num2str(d) '= tmp8(d,3)']);
% %     
    fsbmmf_OutImg = fsbmmf(nImg);                       % Call fsbmmf median filter 2014 AEU
    tmp7(d,4)=psnr(Img,fsbmmf_OutImg);                  % Calculate PSNR of fsbmmf  
    eval(['psnr_fsbmmf' num2str(d) '= tmp7(d,4)']);
    tmp8(d,4)=ssim(Img,fsbmmf_OutImg);                  % Calculate SSIM
    eval(['ssim_fsbmmf' num2str(d) '= tmp8(d,4)']);
%     
    RSIF_OutImg = RSIF(nImg);                           % Call RSIF median filter 2014 SIViP
    tmp7(d,5)=psnr(Img,RSIF_OutImg);                    % Calculate PSNR of fsbmmf  
    eval(['psnr_RSIF' num2str(d) '= tmp7(d,5)']);
    tmp8(d,5)=ssim(Img,RSIF_OutImg);                    % Calculate SSIM
    eval(['ssim_RSIF' num2str(d) '= tmp8(d,5)']);

    DAMF_OutImg = DAMF(nImg);                           % Call DAMF median filter 2018 AEU 
    tmp7(d,7)=psnr(Img,DAMF_OutImg);                    % Calculate PSNR of DAMF  
    eval(['psnr_DAMF' num2str(d) '= tmp7(d,7)']);
    tmp8(d,7)=ssim(Img,DAMF_OutImg);                    % Calculate SSIM
    eval(['ssim_DAMF' num2str(d) '= tmp8(d,7)']);
   
      TVWA_OutImg = PAtern(nImg);               % TVWA
      tmp7(d,13)=psnr(Img,TVWA_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_TVWA' num2str(d) '= tmp7(d,13)']);
      tmp8(d,13)=ssim(Img,TVWA_OutImg);                    % Calculate SSIM
      eval(['ssim_TVWA' num2str(d) '= tmp8(d,13)']);
     
      
      new_output = newAlgo(nImg);
      tmp7(d,17)=psnr(Img,new_output);                    % Calculate PSNR of TVWA  
      eval(['psnr_new_output' num2str(d) '= tmp7(d,17)']);
      tmp8(d,17)=ssim(Img,new_output);                    % Calculate SSIM
      eval(['ssim_new_output' num2str(d) '= tmp8(d,17)']);
      
  

a=a+0.02
end

tmp1=tmp7;
tmp2=tmp8;


% 
% % Code to print PSNR
x= [85 87 89 91 93 95 97 99];
% x= linespace(85,99,2)
figure(1);
 plot(x,tmp1(:,1),':b*'); hold on;     % dotted line
plot(x,tmp1(:,2),'--ko'); hold on;    % dashed line
plot(x,tmp1(:,3),'-.md'); hold on;  % solid line with diamond specifier
 plot(x,tmp1(:,4),'-.b^'); hold on;    % 
 plot(x,tmp1(:,5),'-.r^'); hold on;
% plot(x,tmp1(:,6),'-.gs'); hold on;
  plot(x,tmp1(:,7),'-.g^'); hold on;
% % plot(x,tmp1(:,8),'-.b^'); hold on; 
% % plot(x,tmp1(:,9),'-.r^'); hold on;
% %  plot(x,tmp1(:,10),'-.r*'); hold on;s
% % plot(x,tmp1(:,11),'-.rs'); hold on;
 plot(x,tmp1(:,13),'--k<'); hold on;
%  plot(x,tmp1(:,12),'-.rs'); hold on;
%  plot(x,tmp1(:,14),':b*'); hold on;
%  plot(x,tmp1(:,15),'-.md'); hold on;
  plot(x,tmp1(:,17),'-.cx'); hold on;
%  plot(x,tmp1(:,14),'-.ys');
legend('DBAMF','ASWMF','MDBUTM','FSBMMF','RSIF','DAMF','TVWA','PA', 'Location', 'NorthEast');
% %    'fsbmmf','UTMF','DAMF','BPDM','DBAMF','ASWMF','mdbutm','RSIF','pdbm','DBAMF','UTMP',,'
% title('PSNR');
% 'DBAMF','ASWMF', 'MDBUTM','FSBMMF','RSIF','DAMF',
xlabel('Noise Density (%)','FontSize',15,'FontName','Times New Roman');
ylabel('PSNR (db)','FontSize',15,'FontName','Times New Roman');


% % Code to print PSNR
x= [85 87 89 91 93 95 97 99];
% x= linespace(85,99,2)
figure(2);
 plot(x,tmp2(:,1),':b*'); hold on;     % dotted line
plot(x,tmp2(:,2),'--ko'); hold on;    % dashed line
plot(x,tmp2(:,3),'-.md'); hold on;  % solid line with diamond specifier
 plot(x,tmp2(:,4),'-.b^'); hold on;    % 
 plot(x,tmp2(:,5),'-.r^'); hold on;
% plot(x,tmp1(:,6),'-.gs'); hold on;
  plot(x,tmp2(:,7),'-.g^'); hold on;
% % plot(x,tmp1(:,8),'-.b^'); hold on; 
% % plot(x,tmp1(:,9),'-.r^'); hold on;
% %  plot(x,tmp1(:,10),'-.r*'); hold on;s
% % plot(x,tmp1(:,11),'-.rs'); hold on;
 plot(x,tmp2(:,13),'--k<'); hold on;
%  plot(x,tmp1(:,12),'-.rs'); hold on;
%  plot(x,tmp1(:,14),':b*'); hold on;
%  plot(x,tmp1(:,15),'-.md'); hold on;
  plot(x,tmp2(:,17),'-.cx'); hold on;
%  plot(x,tmp1(:,14),'-.ys');
legend('DBAMF','ASWMF','MDBUTM','FSBMMF','RSIF','DAMF','TVWA','PA', 'Location', 'NorthEast');
% %    'fsbmmf','UTMF','DAMF','BPDM','DBAMF','ASWMF','mdbutm','RSIF','pdbm','DBAMF','UTMP',,'
% title('PSNR');
% 'DBAMF','ASWMF', 'MDBUTM','FSBMMF','RSIF','DAMF',
xlabel('Noise Density (%)','FontSize',15,'FontName','Times New Roman');
ylabel('SSIM','FontSize',15,'FontName','Times New Roman');

 