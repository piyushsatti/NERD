%%% Function to run all mf functions on benchmark images with varying noise
%%% density.

clc; clear all;                                 % Clear all existing variables

Img= imread('lena_gray_512.tif');               % Reading input image
% d=0.5;                                          % Default value of Noise density
x= [1 2 3 4 5 6 7 8 9];
for d= 1:1   
    nImg= imnoise(Img,'salt & pepper',.1*d);   % Introducing noise
% 
%      DBAMF_OutImg = DBAMF(nImg);                         % Call DBAMF median filter 2007 SPL
%     tmp1(d,1)=psnr(Img,DBAMF_OutImg);                   % Calculate PSNR of DBAMF  
%     eval(['psnr_DBAMF' num2str(d) '= tmp1(d,1)']);
%     tmp2(d,1)=ssim(Img,DBAMF_OutImg);                   % Calculate SSIM
%     eval(['ssim_DBAMF' num2str(d) '= tmp2(d,1)']);
    
tic;
    ASWMF_OutImg = ASWMF(nImg);                 % Call ASWMF median filter 2010 ICIIS
    tmp1(d,2)=psnr(Img,ASWMF_OutImg);               % Calculate PSNR of fsbmmf  
    eval(['psnr_ASWMF' num2str(d) '= tmp1(d,2)']);
    tmp2(d,2)=ssim(Img,ASWMF_OutImg);               % Calculate SSIM
    eval(['ssim_ASWMF' num2str(d) '= tmp2(d,2)']);
% %
time(d,1) = toc;

tic
    mdbutm_OutImg = mdbutm(nImg);                       % Call mdbutm median filter 2011 SPL 
    tmp1(d,3)=psnr(Img,mdbutm_OutImg);                  % Calculate PSNR of mdbutm  
    eval(['psnr_mdbutm' num2str(d) '= tmp1(d,3)']);
    tmp2(d,3)=ssim(Img,mdbutm_OutImg);                  % Calculate SSIM
    eval(['ssim_mdbutm' num2str(d) '= tmp2(d,3)']);
time(d,2) = toc;

tic;
    fsbmmf_OutImg = fsbmmf(nImg);                       % Call fsbmmf median filter 2014 AEU
    tmp1(d,4)=psnr(Img,fsbmmf_OutImg);                  % Calculate PSNR of fsbmmf  
    eval(['psnr_fsbmmf' num2str(d) '= tmp1(d,4)']);
    tmp2(d,4)=ssim(Img,fsbmmf_OutImg);                  % Calculate SSIM
    eval(['ssim_fsbmmf' num2str(d) '= tmp2(d,4)']);
    time(d,3) = toc;
    


    
    tic
    RSIF_OutImg = RSIF(nImg);                           % Call RSIF median filter 2014 SIViP
    tmp1(d,5)=psnr(Img,RSIF_OutImg);                    % Calculate PSNR of fsbmmf  
    eval(['psnr_RSIF' num2str(d) '= tmp1(d,5)']);
    tmp2(d,5)=ssim(Img,RSIF_OutImg);                    % Calculate SSIM
    eval(['ssim_RSIF' num2str(d) '= tmp2(d,5)']);
time(d,4) = toc;
    % %     
% % %     pdbm_OutImg = pdbm(nImg);                           % Call pdbm median filter 2016 AEU
% % %     tmp1(d,6)=psnr(Img,pdbm_OutImg);                    % Calculate PSNR of pdbm  
% % %     eval(['psnr_pdbm' num2str(d) '= tmp1(d,6)']);
% % %     tmp2(d,6)=ssim(Img,pdbm_OutImg);                    % Calculate SSIM
% % %     eval(['ssim_pdbm' num2str(d) '= tmp2(d,6)']);
%     
tic
    DAMF_OutImg = DAMF(nImg);                           % Call DAMF median filter 2018 AEU 
    tmp1(d,7)=psnr(Img,DAMF_OutImg);                    % Calculate PSNR of DAMF  
    eval(['psnr_DAMF' num2str(d) '= tmp1(d,7)']);
    tmp2(d,7)=ssim(Img,DAMF_OutImg);                    % Calculate SSIM
    eval(['ssim_DAMF' num2str(d) '= tmp2(d,7)']);
time(d,5) = toc;
% %       modified_final_code1_OutImg = modified_final_code1(nImg);               % Final Modified Code
% %       tmp1(d,12)=psnr(Img,modified_final_code1_OutImg);                    % Calculate PSNR of BPDM  
% %       eval(['psnr_modified_final_code1' num2str(d) '= tmp1(d,12)']);
% %       tmp2(d,12)=ssim(Img,modified_final_code1_OutImg);                    % Calculate SSIM
% %       eval(['ssim_modified_final_code1' num2str(d) '= tmp2(d,12)']);
%       
tic
      TVWA_OutImg = PAtern(nImg);               % TVWA
      tmp1(d,13)=psnr(Img,TVWA_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_TVWA' num2str(d) '= tmp1(d,13)']);
      tmp2(d,13)=ssim(Img,TVWA_OutImg);                    % Calculate SSIM
      eval(['ssim_TVWA' num2str(d) '= tmp2(d,13)']);
time(d,6) = toc;
      % %       TVWA2_OutImg = TVWA2(nImg);               % TVWA
% %       tmp1(d,15)=psnr(Img,TVWA2_OutImg);                    % Calculate PSNR of TVWA  
% %       eval(['psnr_TVWA2' num2str(d) '= tmp7(d,15)']);
% %       tmp2(d,15)=ssim(Img,TVWA2_OutImg);                    % Calculate SSIM
% %       eval(['ssim_TVWA2' num2str(d) '= tmp8(d,15)']);
%       
%  TVWA2final_OutImg = TVWA2final(nImg);
%       tmp1(d,17)=psnr(Img,TVWA2final_OutImg);                    % Calculate PSNR of TVWA  
%       eval(['psnr_TVWA2final' num2str(d) '= tmp1(d,17)']);
%       tmp2(d,17)=ssim(Img,TVWA2final_OutImg);                    % Calculate SSIM
%       eval(['ssim_TVWA2final' num2str(d) '= tmp2(d,17)']);
%                  
tic     
occo_OutImg = occo(nImg);
      tmp1(d,18)=psnr(Img,occo_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_occo' num2str(d) '= tmp1(d,18)']);
      tmp2(d,18)=ssim(Img,occo_OutImg);                    % Calculate SSIM
      eval(['ssim_occo' num2str(d) '= tmp2(d,18)']);
    time(d,7) = toc;  
    
    tic
      morphology_mean_filter_OutImg = morphology_mean_filter(nImg);
      tmp1(d,19)=psnr(Img,morphology_mean_filter_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_morphology_mean_filter' num2str(d) '= tmp1(d,19)']);
      tmp2(d,19)=ssim(Img,morphology_mean_filter_OutImg);                    % Calculate SSIM
      eval(['ssim_morphology_mean_filter' num2str(d) '= tmp2(d,19)']);
      time(d,8) = toc;
      
      tic
     UWMF_OutImg = UWMF(nImg);
      tmp1(d,20)=psnr(Img,UWMF_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_UWMF' num2str(d) '= tmp1(d,20)']);
      tmp2(d,20)=ssim(Img,UWMF_OutImg);                    % Calculate SSIM
      eval(['ssim_UWMF' num2str(d) '= tmp2(d,20)']);
       time(d,9) = toc;
end

% 
% % Code to print PSNR
figure(1);
x=[10 20 30 40 50 60 70 80 90]
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
plot(x,tmp1(:,17),'-.rs'); hold on;
plot(x,tmp1(:,18),'-.k*'); hold on;
plot(x,tmp1(:,19),'-.c*'); hold on;
plot(x,tmp1(:,20),'--c^'); hold on;
legend('ASWMF', 'MDBUTM','FSBMMF','RSIF','DAMF','TVWA','2ndPaper','occo','morphology_mean_filter','UWMF','Location', 'NorthEast');
% %    'fsbmmf','UTMF','DAMF','BPDM','DMF','ASWMF','mdbutm','RSIF','pdbm','DMF','UTMP',,'
% title('PSNR');
xlabel('Noise Density (%)','FontSize',15,'FontName','Times New Roman');
ylabel('PSNR (db)','FontSize',15,'FontName','Times New Roman');


% % % Code to print PSNR
% figure(2);
%  plot(x,tmp8(:,1),':b*'); hold on;     % dotted line
% plot(x,tmp2(:,2),'--ko'); hold on;    % dashed line
% plot(x,tmp2(:,3),'-.md'); hold on;  % solid line with diamond specifier
%  plot(x,tmp2(:,4),'-.b^'); hold on;    % 
%  plot(x,tmp2(:,5),'-.cx'); hold on;
% % plot(x,tmp2(:,6),'-.gs'); hold on;
% plot(x,tmp2(:,7),'-.r^'); hold on;
% % % plot(x,tmp2(:,8),'-.b^'); hold on; 
% % % plot(x,tmp2(:,9),'-.r^'); hold on;
% %  plot(x,tmp2(:,10),'-.r*'); hold on;s
% % plot(x,tmp2(:,11),'-.rs'); hold on;
% %  plot(x,tmp2(:,13),'-.k<'); hold on;
%  plot(x,tmp2(:,12),'-.rs'); hold on;
% %  plot(x,tmp2(:,14),'-.ys');
%  plot(x,tmp2(:,12),'-.rs'); hold on;
% legend('DBAMF','ASWMF', 'MDBUTM','FSBMMF','RSIF','DAMF','TVWA','PA', 'Location', 'NorthEast');
% % %    'fsbmmf','UTMF','DAMF','BPDM','DMF','ASWMF','mdbutm','RSIF','pdbm','DMF','UTMP',,'
% % title('SSIM');
% xlabel('Noise Density (%)','FontSize',15,'FontName','Times New Roman');
% ylabel('SSIM','FontSize',15,'FontName','Times New Roman');


