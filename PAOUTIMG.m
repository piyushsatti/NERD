%%% Function to run all mf functions on benchmark images with varying noise
%%% density.

clc; clear all;                                 % Clear all existing variables;
a=0.91;                                          % Default value of Noise density
for i=1:1
  if i==1
  Img=  imread('lena_gray_512.tif');
  elseif i==2
      a=.85;
  Img= imread('zelda.png');
  else
      a=.85;
  Img= imread('lena_gray_512.tif');
  end
    
    
    nImg = imnoise(Img,'salt & pepper',a);   % Introducing noise
    if i==1
%     DBAMF_OutImg = DBAMF(nImg);                         % Call DBAMF median filter 2007 SPL
%     tmp7(d,1) = psnr(Img,DBAMF_OutImg);                   % Calculate PSNR of DBAMF  
%     eval(['psnr_DBAMF' num2str(d) '= tmp7(d,1)']);
%     tmp8(d,1)=ssim(Img,DBAMF_OutImg);                   % Calculate SSIM
%     eval(['ssim_DBAMF' num2str(d) '= tmp8(d,1)']);
%     
%     ASWMF_OutImg = ASWMF(nImg);                 % Call ASWMF median filter 2010 ICIIS
%     tmp7(d,2)=psnr(Img,ASWMF_OutImg);               % Calculate PSNR of fsbmmf  
%     eval(['psnr_ASWMF' num2str(d) '= tmp7(d,2)']);
%     tmp8(d,2)=ssim(Img,ASWMF_OutImg);               % Calculate SSIM
%     eval(['ssim_ASWMF' num2str(d) '= tmp8(d,2)']);
% % %     
%     mdbutm_OutImg = mdbutm(nImg);                       % Call mdbutm median filter 2011 SPL 
%     tmp7(d,3)=psnr(Img,mdbutm_OutImg);                  % Calculate PSNR of mdbutm  
%     eval(['psnr_mdbutm' num2str(d) '= tmp7(d,3)']);
%     tmp8(d,3)=ssim(Img,mdbutm_OutImg);                  % Calculate SSIM
%     eval(['ssim_mdbutm' num2str(d) '= tmp8(d,3)']);
% % %     
%     fsbmmf_OutImg = fsbmmf(nImg);                       % Call fsbmmf median filter 2014 AEU
%     tmp7(d,4)=psnr(Img,fsbmmf_OutImg);                  % Calculate PSNR of fsbmmf  
%     eval(['psnr_fsbmmf' num2str(d) '= tmp7(d,4)']);
%     tmp8(d,4)=ssim(Img,fsbmmf_OutImg);                  % Calculate SSIM
%     eval(['ssim_fsbmmf' num2str(d) '= tmp8(d,4)']);
% %     
%     RSIF_OutImg = RSIF(nImg);                           % Call RSIF median filter 2014 SIViP
%     tmp7(d,5)=psnr(Img,RSIF_OutImg);                    % Calculate PSNR of fsbmmf  
%     eval(['psnr_RSIF' num2str(d) '= tmp7(d,5)']);
%     tmp8(d,5)=ssim(Img,RSIF_OutImg);                    % Calculate SSIM
%     eval(['ssim_RSIF' num2str(d) '= tmp8(d,5)']);
% % %     
% % %     pdbm_OutImg = pdbm(nImg);                           % Call pdbm median filter 2016 AEU
% % %     tmp1(d,6)=psnr(Img,pdbm_OutImg);                    % Calculate PSNR of pdbm  
% % %     eval(['psnr_pdbm' num2str(d) '= tmp1(d,6)']);
% % %     tmp2(d,6)=ssim(Img,pdbm_OutImg);                    % Calculate SSIM
% % %     eval(['ssim_pdbm' num2str(d) '= tmp2(d,6)']);
% %     
%     DAMF_OutImg = DAMF(nImg);                           % Call DAMF median filter 2018 AEU 
%     tmp7(d,7)=psnr(Img,DAMF_OutImg);                    % Calculate PSNR of DAMF  
%     eval(['psnr_DAMF' num2str(d) '= tmp7(d,7)']);
%     tmp8(d,7)=ssim(Img,DAMF_OutImg);                    % Calculate SSIM
%     eval(['ssim_DAMF' num2str(d) '= tmp8(d,7)']);
% 
% %       modified_final_code1_OutImg = modified_final_code1(nImg);               % Final Modified Code
% %       tmp7(d,12)=psnr(Img,modified_final_code1_OutImg);                    % Calculate PSNR of BPDM  
% %       eval(['psnr_modified_final_code1' num2str(d) '= tmp7(d,12)']);
% %       tmp8(d,12)=ssim(Img,modified_final_code1_OutImg);                    % Calculate SSIM
% %       eval(['ssim_modified_final_code1' num2str(d) '= tmp8(d,12)']);
%       
%       TVWA_OutImg = PAtern(nImg);               % TVWA
%       tmp7(d,13)=psnr(Img,TVWA_OutImg);                    % Calculate PSNR of TVWA  
%       eval(['psnr_TVWA' num2str(d) '= tmp7(d,13)']);
%       tmp8(d,13)=ssim(Img,TVWA_OutImg);                    % Calculate SSIM
%       eval(['ssim_TVWA' num2str(d) '= tmp8(d,13)']);
% %       
      
      TVWA2final_OutImg = TVWA2final(nImg);
%       tmp7(d,17)=psnr(Img,TVWA2final_OutImg);                    % Calculate PSNR of TVWA  
%       eval(['psnr_TVWA2final' num2str(d) '= tmp7(d,17)']);
%       tmp8(d,17)=ssim(Img,TVWA2final_OutImg);                    % Calculate SSIM
%       eval(['ssim_TVWA2final' num2str(d) '= tmp8(d,17)']);
      
      
      
      
      
%       
%       
    elseif i==2
%     DBAMF_OutImg = DBAMF(nImg);                         % Call DBAMF median filter 2007 SPL
%     tmp3(d,1)=psnr(Img,DBAMF_OutImg);                   % Calculate PSNR of DBAMF  
%     eval(['psnr_DBAMF' num2str(d) '= tmp3(d,1)']);
%     tmp4(d,1)=ssim(Img,DBAMF_OutImg);                   % Calculate SSIM
%     eval(['ssim_DBAMF' num2str(d) '= tmp4(d,1)']);
%     
%     ASWMF_OutImg = ASWMF(nImg);                 % Call ASWMF median filter 2010 ICIIS
%     tmp3(d,2)=psnr(Img,ASWMF_OutImg);               % Calculate PSNR of fsbmmf  
%     eval(['psnr_ASWMF' num2str(d) '= tmp3(d,2)']);
%     tmp4(d,2)=ssim(Img,ASWMF_OutImg);               % Calculate SSIM
%     eval(['ssim_ASWMF' num2str(d) '= tmp4(d,2)']);
% %     
%     mdbutm_OutImg = mdbutm(nImg);                       % Call mdbutm median filter 2011 SPL 
%     tmp3(d,3)=psnr(Img,mdbutm_OutImg);                  % Calculate PSNR of mdbutm  
%     eval(['psnr_mdbutm' num2str(d) '= tmp3(d,3)']);
%     tmp4(d,3)=ssim(Img,mdbutm_OutImg);                  % Calculate SSIM
%     eval(['ssim_mdbutm' num2str(d) '= tmp4(d,3)']);
% %     
%     fsbmmf_OutImg = fsbmmf(nImg);                       % Call fsbmmf median filter 2014 AEU
%     tmp3(d,4)=psnr(Img,fsbmmf_OutImg);                  % Calculate PSNR of fsbmmf  
%     eval(['psnr_fsbmmf' num2str(d) '= tmp3(d,4)']);
%     tmp4(d,4)=ssim(Img,fsbmmf_OutImg);                  % Calculate SSIM
%     eval(['ssim_fsbmmf' num2str(d) '= tmp4(d,4)']);
%     
%     RSIF_OutImg = RSIF(nImg);                           % Call RSIF median filter 2014 SIViP
%     tmp3(d,5)=psnr(Img,RSIF_OutImg);                    % Calculate PSNR of fsbmmf  
%     eval(['psnr_RSIF' num2str(d) '= tmp3(d,5)']);
%     tmp4(d,5)=ssim(Img,RSIF_OutImg);                    % Calculate SSIM
%     eval(['ssim_RSIF' num2str(d) '= tmp4(d,5)']);
% %     
% %     pdbm_OutImg = pdbm(nImg);                           % Call pdbm median filter 2016 AEU
% %     tmp1(d,6)=psnr(Img,pdbm_OutImg);                    % Calculate PSNR of pdbm  
% %     eval(['psnr_pdbm' num2str(d) '= tmp1(d,6)']);
% %     tmp2(d,6)=ssim(Img,pdbm_OutImg);                    % Calculate SSIM
% %     eval(['ssim_pdbm' num2str(d) '= tmp2(d,6)']);
%     
%     DAMF_OutImg = DAMF(nImg);                           % Call DAMF median filter 2018 AEU 
%     tmp3(d,7)=psnr(Img,DAMF_OutImg);                    % Calculate PSNR of DAMF  
%     eval(['psnr_DAMF' num2str(d) '= tmp3(d,7)']);
%     tmp4(d,7)=ssim(Img,DAMF_OutImg);                    % Calculate SSIM
%     eval(['ssim_DAMF' num2str(d) '= tmp4(d,7)']);

      modified_final_code1_OutImg = modified_final_code1(nImg);               % Final Modified Code
      tmp3(d,12)=psnr(Img,modified_final_code1_OutImg);                    % Calculate PSNR of BPDM  
      eval(['psnr_modified_final_code1' num2str(d) '= tmp3(d,12)']);
      tmp4(d,12)=ssim(Img,modified_final_code1_OutImg);                    % Calculate SSIM
      eval(['ssim_modified_final_code1' num2str(d) '= tmp4(d,12)']);
      
      TVWA_OutImg = PAtern(nImg);               % TVWA
      tmp3(d,13)=psnr(Img,TVWA_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_TVWA' num2str(d) '= tmp3(d,13)']);
      tmp4(d,13)=ssim(Img,TVWA_OutImg);                    % Calculate SSIM
      eval(['ssim_TVWA' num2str(d) '= tmp4(d,13)']);
      
%       
%       NAMF_OutImg = NAMF(nImg);               % TVWA
%       tmp1(d,14)=psnr(Img,NAMF_OutImg);                    % Calculate PSNR of TVWA  
%       eval(['psnr_NAMF' num2str(d) '= tmp1(d,14)']);
%       tmp2(d,14)=ssim(Img,NAMF_OutImg);                    % Calculate SSIM
%       eval(['ssim_NAMF' num2str(d) '= tmp2(d,14)']);
%     
    else
%         DBAMF_OutImg = DBAMF(nImg);                         % Call DBAMF median filter 2007 SPL
%     tmp5(d,1)=psnr(Img,DBAMF_OutImg);                   % Calculate PSNR of DBAMF  
%     eval(['psnr_DBAMF' num2str(d) '= tmp5(d,1)']);
%     tmp6(d,1)=ssim(Img,DBAMF_OutImg);                   % Calculate SSIM
%     eval(['ssim_DBAMF' num2str(d) '= tmp6(d,1)']);
%     
%     ASWMF_OutImg = ASWMF(nImg);                 % Call ASWMF median filter 2010 ICIIS
%     tmp5(d,2)=psnr(Img,ASWMF_OutImg);               % Calculate PSNR of fsbmmf  
%     eval(['psnr_ASWMF' num2str(d) '= tmp5(d,2)']);
%     tmp6(d,2)=ssim(Img,ASWMF_OutImg);               % Calculate SSIM
%     eval(['ssim_ASWMF' num2str(d) '= tmp6(d,2)']);
% %     
%     mdbutm_OutImg = mdbutm(nImg);                       % Call mdbutm median filter 2011 SPL 
%     tmp5(d,3)=psnr(Img,mdbutm_OutImg);                  % Calculate PSNR of mdbutm  
%     eval(['psnr_mdbutm' num2str(d) '= tmp5(d,3)']);
%     tmp6(d,3)=ssim(Img,mdbutm_OutImg);                  % Calculate SSIM
%     eval(['ssim_mdbutm' num2str(d) '= tmp6(d,3)']);
% %     
%     fsbmmf_OutImg = fsbmmf(nImg);                       % Call fsbmmf median filter 2014 AEU
%     tmp5(d,4)=psnr(Img,fsbmmf_OutImg);                  % Calculate PSNR of fsbmmf  
%     eval(['psnr_fsbmmf' num2str(d) '= tmp5(d,4)']);
%     tmp6(d,4)=ssim(Img,fsbmmf_OutImg);                  % Calculate SSIM
%     eval(['ssim_fsbmmf' num2str(d) '= tmp6(d,4)']);
%     
%     RSIF_OutImg = RSIF(nImg);                           % Call RSIF median filter 2014 SIViP
%     tmp5(d,5)=psnr(Img,RSIF_OutImg);                    % Calculate PSNR of fsbmmf  
%     eval(['psnr_RSIF' num2str(d) '= tmp5(d,5)']);
%     tmp6(d,5)=ssim(Img,RSIF_OutImg);                    % Calculate SSIM
%     eval(['ssim_RSIF' num2str(d) '= tmp6(d,5)']);
% %     
% %     pdbm_OutImg = pdbm(nImg);                           % Call pdbm median filter 2016 AEU
% %     tmp1(d,6)=psnr(Img,pdbm_OutImg);                    % Calculate PSNR of pdbm  
% %     eval(['psnr_pdbm' num2str(d) '= tmp1(d,6)']);
% %     tmp2(d,6)=ssim(Img,pdbm_OutImg);                    % Calculate SSIM
% %     eval(['ssim_pdbm' num2str(d) '= tmp2(d,6)']);
%     
%     DAMF_OutImg = DAMF(nImg);                           % Call DAMF median filter 2018 AEU 
%     tmp5(d,7)=psnr(Img,DAMF_OutImg);                    % Calculate PSNR of DAMF  
%     eval(['psnr_DAMF' num2str(d) '= tmp5(d,7)']);
%     tmp6(d,7)=ssim(Img,DAMF_OutImg);                    % Calculate SSIM
%     eval(['ssim_DAMF' num2str(d) '= tmp6(d,7)']);
% 
      modified_final_code1_OutImg = modified_final_code1(nImg);               % Final Modified Code
      tmp5(d,12)=psnr(Img,modified_final_code1_OutImg);                    % Calculate PSNR of BPDM  
      eval(['psnr_modified_final_code1' num2str(d) '= tmp5(d,12)']);
      tmp6(d,12)=ssim(Img,modified_final_code1_OutImg);                    % Calculate SSIM
      eval(['ssim_modified_final_code1' num2str(d) '= tmp6(d,12)']);
      
      TVWA_OutImg = PAtern(nImg);               % TVWA
      tmp5(d,13)=psnr(Img,TVWA_OutImg);                    % Calculate PSNR of TVWA  
      eval(['psnr_TVWA' num2str(d) '= tmp5(d,13)']);
      tmp6(d,13)=ssim(Img,TVWA_OutImg);                    % Calculate SSIM
      eval(['ssim_TVWA' num2str(d) '= tmp6(d,13)']);
      
%       
%       NAMF_OutImg = NAMF(nImg);               % TVWA
%       tmp1(d,14)=psnr(Img,NAMF_OutImg);                    % Calculate PSNR of TVWA  
%       eval(['psnr_NAMF' num2str(d) '= tmp1(d,14)']);
%       tmp2(d,14)=ssim(Img,NAMF_OutImg);                    % Calculate SSIM
%       eval(['ssim_NAMF' num2str(d) '= tmp2(d,14)']);
%     
    end



end
% 

% t

imwrite(TVWA2final_OutImg,'PAOUTIMG.tif');