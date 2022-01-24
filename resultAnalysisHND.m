%%% Function to run all mf functions on benchmark images with varying noise
%%% density.

clc; clear all;                                 % Clear all existing variables

% Reading input image
% IMG=[Img1 Img2 Img3];
% Img = imread('lena_gray_512.tif');
a=0.5;                                          % Default value of Noise density
for i=1:1
%   if i==1
  Img=  imread('zelda.png');
%   elseif i==5
%       a=.75;
%   Img= imread('zelda.png');
%   else
%       a=.75;
%   Img= imread('lena_gray_512.tif');
%   end
%     
%     
for d= 1:10
    
    nImg= imnoise(Img,'salt & pepper',a);   % Introducing noise
   
    if i==1
%      DBAMF_OutImg = DBAMF(nImg);                         % Call DBAMF median filter 2007 SPL
%     tmp1(d,1)=psnr(Img,DBAMF_OutImg);                   % Calculate PSNR of DBAMF  
%     eval(['psnr_DBAMF' num2str(d) '= tmp1(d,1)']);
%     tmp2(d,1)=ssim(Img,DBAMF_OutImg);                   % Calculate SSIM
%     eval(['ssim_DBAMF' num2str(d) '= tmp2(d,1)']);
    
    ASWMF_OutImg = ASWMF(nImg);                 % Call ASWMF median filter 2010 ICIIS
    tmp1(d,2)=psnr(Img,ASWMF_OutImg);               % Calculate PSNR of fsbmmf  
    eval(['psnr_ASWMF' num2str(d) '= tmp1(d,2)']);
    tmp2(d,2)=ssim(Img,ASWMF_OutImg);               % Calculate SSIM
    eval(['ssim_ASWMF' num2str(d) '= tmp2(d,2)']);
% %     
    mdbutm_OutImg = mdbutm(nImg);                       % Call mdbutm median filter 2011 SPL 
    tmp1(d,3)=psnr(Img,mdbutm_OutImg);                  % Calculate PSNR of mdbutm  
    eval(['psnr_mdbutm' num2str(d) '= tmp1(d,3)']);
    tmp2(d,3)=ssim(Img,mdbutm_OutImg);                  % Calculate SSIM
    eval(['ssim_mdbutm' num2str(d) '= tmp2(d,3)']);
%     
    fsbmmf_OutImg = fsbmmf(nImg);                       % Call fsbmmf median filter 2014 AEU
    tmp1(d,4)=psnr(Img,fsbmmf_OutImg);                  % Calculate PSNR of fsbmmf  
    eval(['psnr_fsbmmf' num2str(d) '= tmp1(d,4)']);
    tmp2(d,4)=ssim(Img,fsbmmf_OutImg);                  % Calculate SSIM
    eval(['ssim_fsbmmf' num2str(d) '= tmp2(d,4)']);
    
    RSIF_OutImg = RSIF(nImg);                           % Call RSIF median filter 2014 SIViP
    tmp1(d,5)=psnr(Img,RSIF_OutImg);                    % Calculate PSNR of fsbmmf  
    eval(['psnr_RSIF' num2str(d) '= tmp1(d,5)']);
    tmp2(d,5)=ssim(Img,RSIF_OutImg);                    % Calculate SSIM
    eval(['ssim_RSIF' num2str(d) '= tmp2(d,5)']);
% %     
% %     pdbm_OutImg = pdbm(nImg);                           % Call pdbm median filter 2016 AEU
% %     tmp1(d,6)=psnr(Img,pdbm_OutImg);                    % Calculate PSNR of pdbm  
% %     eval(['psnr_pdbm' num2str(d) '= tmp1(d,6)']);
% %     tmp2(d,6)=ssim(Img,pdbm_OutImg);                    % Calculate SSIM
% %     eval(['ssim_pdbm' num2str(d) '= tmp2(d,6)']);
    
    DAMF_OutImg = DAMF(nImg);                           % Call DAMF median filter 2018 AEU 
    tmp1(d,7)=psnr(Img,DAMF_OutImg);                    % Calculate PSNR of DAMF  
    eval(['psnr_DAMF' num2str(d) '= tmp1(d,7)']);
    tmp2(d,7)=ssim(Img,DAMF_OutImg);                    % Calculate SSIM
    eval(['ssim_DAMF' num2str(d) '= tmp2(d,7)']);

%       modified_final_code1_OutImg = modified_final_code1(nImg);               % Final Modified Code
%       tmp1(d,12)=psnr(Img,modified_final_code1_OutImg);                    % Calculate PSNR of BPDM  
%       eval(['psnr_modified_final_code1' num2str(d) '= tmp1(d,12)']);
%       tmp2(d,12)=ssim(Img,modified_final_code1_OutImg);                    % Calculate SSIM
%       eval(['ssim_modified_final_code1' num2str(d) '= tmp2(d,12)']);
      
      TVWA_OutImg = PAtern(nImg);               % TVWA
      tmp1(d,13)=psnr(Img,TVWA_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_TVWA' num2str(d) '= tmp1(d,13)']);
      tmp2(d,13)=ssim(Img,TVWA_OutImg);                    % Calculate SSIM
      eval(['ssim_TVWA' num2str(d) '= tmp2(d,13)']);
%       TVWA2_OutImg = TVWA2(nImg);               % TVWA
%       tmp1(d,15)=psnr(Img,TVWA2_OutImg);                    % Calculate PSNR of TVWA  
%       eval(['psnr_TVWA2' num2str(d) '= tmp7(d,15)']);
%       tmp2(d,15)=ssim(Img,TVWA2_OutImg);                    % Calculate SSIM
%       eval(['ssim_TVWA2' num2str(d) '= tmp8(d,15)']);
      
             TVWA2final_OutImg = TVWA2final(nImg);
      tmp1(d,17)=psnr(Img,TVWA2final_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_TVWA2final' num2str(d) '= tmp1(d,17)']);
      tmp2(d,17)=ssim(Img,TVWA2final_OutImg);                    % Calculate SSIM
      eval(['ssim_TVWA2final' num2str(d) '= tmp2(d,17)']);
                 
      occo_OutImg = occo(nImg);
      tmp1(d,18)=psnr(Img,occo_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_occo' num2str(d) '= tmp1(d,18)']);
      tmp2(d,18)=ssim(Img,occo_OutImg);                    % Calculate SSIM
      eval(['ssim_occo' num2str(d) '= tmp2(d,18)']);
      
      morphology_mean_filter_OutImg = morphology_mean_filter(nImg);
      tmp1(d,19)=psnr(Img,morphology_mean_filter_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_morphology_mean_filter' num2str(d) '= tmp1(d,19)']);
      tmp2(d,19)=ssim(Img,morphology_mean_filter_OutImg);                    % Calculate SSIM
      eval(['ssim_morphology_mean_filter' num2str(d) '= tmp2(d,19)']);
%       
     UWMF_OutImg = UWMF(nImg);
      tmp1(d,20)=psnr(Img,UWMF_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_UWMF' num2str(d) '= tmp1(d,20)']);
      tmp2(d,20)=ssim(Img,UWMF_OutImg);                    % Calculate SSIM
      eval(['ssim_UWMF' num2str(d) '= tmp2(d,20)']);
%       

    end

a=a+0.05
end
end


tmp1 = tmp7;
tmp2 = tmp8;

figure(1);
x=[50  55 60 65 70 75 80 85 90 95]
%  plot(x,tmp1(:,1),':b*'); hold on;     % dotted line
plot(x,tmp1(:,2),'--ko'); hold on;    % dashed line
plot(x,tmp1(:,3),'-.md'); hold on;  % solid line with diamond specifier
 plot(x,tmp1(:,4),'-.b^'); hold on;    % 
 plot(x,tmp1(:,5),'-.gx'); hold on;
% plot(x,tmp1(:,6),'-.gs'); hold on;
 plot(x,tmp1(:,7),'-.r^'); hold on;
% % plot(x,tmp1(:,8),'-.b^'); hold on; 
% plot(x,tmp1(:,9),'-.r^'); hold on;
%  plot(x,tmp1(:,10),'-.r*'); hold on;s
% plot(x,tmp1(:,11),'-.rs'); hold on;
plot(x,tmp1(:,13),'--k<'); hold on;
%  plot(x,tmp1(:,12),'-.rs'); hold on;
%  plot(x,tmp1(:,14),'-.ys');
% plot(x,tmp1(:,17),'-.rs'); hold on;
plot(x,tmp1(:,18),'-.k*'); hold on;
plot(x,tmp1(:,19),'-.c*'); hold on;
plot(x,tmp1(:,20),'--c^'); hold on;
legend('ASWMF', 'MDBUTM','FSBMMF','RSIF','DAMF','TVWA','occo','morphology_mean_filter','UWMF','Location', 'NorthEast');
% %    'fsbmmf','UTMF','DAMF','BPDM','DMF','ASWMF','mdbutm','RSIF','pdbm','DMF','UTMP',,'
% title('PSNR');
xlabel('Noise Density (%)','FontSize',15,'FontName','Times New Roman');
ylabel('PSNR (db)','FontSize',15,'FontName','Times New Roman');


% x= [75 77 79 81 83 85 87 89 91 93 95 97 99];
% x= linespace(85,99,2)
x= [10 20 30 40 50 60 70 80 90];
figure(2);
 plot(x,tmp2(:,1),':b*'); hold on;     % dotted line
plot(x,tmp2(:,2),'--ko'); hold on;    % dashed line
plot(x,tmp2(:,3),'-.md'); hold on;  % solid line with diamond specifier
 plot(x,tmp2(:,4),'-.b^'); hold on;    % 
 plot(x,tmp2(:,5),'-.cx'); hold on;
% plot(x,tmp1(:,6),'-.gs'); hold on;
  plot(x,tmp2(:,7),'-.g^'); hold on;
% % plot(x,tmp1(:,8),'-.b^'); hold on; 
% plot(x,tmp1(:,9),'-.r^'); hold on;
%  plot(x,tmp1(:,10),'-.r*'); hold on;s
% plot(x,tmp1(:,11),'-.rs'); hold on;
 plot(x,tmp2(:,13),'--k<'); hold on;
%  plot(x,tmp2(:,12),'-.rs'); hold on;
%   plot(x,tmp1(:,14),'-.b^'); hold on;
%    plot(x,tmp2(:,15),'-.rs'); hold on;
%  plot(x,tmp1(:,14),'-.ys');
legend('DBAMF','ASWMF', 'MDBUTM','FSBMMF','RSIF','DAMF','TVWA', 'Location', 'NorthEast');
% %    'fsbmmf','UTMF','DAMF','BPDM','DBAMF','ASWMF','mdbutm','RSIF','pdbm','DBAMF','UTMP',,'
% title('SSIM');
xlabel('Noise Density (%)','FontSize',15,'FontName','Times New Roman');
ylabel('SSIM','FontSize',15,'FontName','Times New Roman');


