%%%%%%%%% the signal components are not crossing %%%%%%%%%%%%%%
% this example is adopted from paper£ºChen S, Peng Z, Yang Y, et al, Intrinsic chirp component decomposition by using Fourier Series representation, Signal Processing, 2017.
clc
clear
close all
SampFreq = 250;
t = 0:1/SampFreq:6;
amp1 = exp(0.2*t);
amp2 = exp(0.15*t);
amp3 = exp(0.1*t);
Sig1 = amp1.*cos(2*pi*(1.3+6*t + 4*t.^2 - 0.8*t.^3 + 0.07*t.^4));
IF1 = 6 + 8*t - 2.4*t.^2 + 0.28*t.^3;
Sig2 = amp2.*cos(2*pi*(2.6+12*t + 8*t.^2 -1.6*t.^3 + 0.14*t.^4));
IF2 = 12 + 16*t - 4.8*t.^2 + 0.56*t.^3;
Sig3 = amp3.*cos(2*pi*(3.9+18*t + 12*t.^2 - 2.4*t.^3 + 0.21*t.^4));
IF3 = 18 + 24*t - 7.2*t.^2 + 0.84*t.^3;
Sig4 = 1./(1.2+cos(2*pi*0.165*t))+0.5; % a trend component

Sig = Sig1 + Sig2 + Sig3 + Sig4;

% figure
% set(gcf,'Position',[20 100 320 250]);	    
% set(gcf,'Color','w'); 
% plot(t,Sig);
% xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
% ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
% set(gca,'FontSize',12)
% set(gca,'linewidth',1);
% axis([0 6 -2 12])
%% STFT
window = 256;
Nfrebin = 1024;
figure
[Spec,f] = STFT(Sig',SampFreq,Nfrebin,window);
imagesc(t,f,abs(Spec)); 
axis([0 6 0 100]);
set(gcf,'Position',[20 100 320 250]);	 
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',12);
set(gcf,'Color','w');	
%% ridge extraction and smoothing
bw = SampFreq/80;% the bandwidth of the TF filter for ridge extraction
beta1 = 1e-4; % beta1 should be larger than the following beta
num = 4; % the number of the components
delta = 20;
[fidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, beta1,bw,Nfrebin,window);

figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,f(fidexmult),'linewidth',3); % detected ridge curves
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 6 0 100]);

% %% parameter setting
% alpha = 1e-5; 
% beta = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
% iniIF = curvesmooth(f(fidexmult),beta); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges
% var = 0;% noise variance
% tol = 1e-7;%
% tic
% [IFmset IA smset] = VNCMD(Sig,SampFreq,iniIF,alpha,beta,var,tol);
% toc
% 
% figure
% set(gcf,'Position',[20 100 640 500]);	    
% set(gcf,'Color','w'); 
% plot(t,iniIF,'linewidth',3); % initial IFs
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
% set(gca,'FontSize',24)
% set(gca,'linewidth',2);
% axis([0 6 0 100]);
% 
% %% Relative errors of the initial IFs
%    iniRE1 =  norm(iniIF(2,:)-IF1)/norm(IF1)
%    iniRE2 =  norm(iniIF(3,:)-IF2)/norm(IF2)
%    iniRE3 =  norm(iniIF(4,:)-IF3)/norm(IF3)
% %% estimated IF
% figure
% plot(t,[IF1;IF2;IF3],'b','linewidth',3) % true IFs
% hold on
% plot(t,IFmset(:,:,end),'r','linewidth',3) % finally estimated IFs
% set(gcf,'Position',[20 100 640 500]);	 
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
% set(gca,'YDir','normal')
% set(gca,'FontSize',24);
% set(gca,'linewidth',2);
% set(gcf,'Color','w');	
% axis([0 6 0 100])
% %% Relative errors of the finally estimated IFs
%    RE1 =  norm(IFmset(2,:,end)-IF1)/norm(IF1)   % one can find that the accuracy of the finally estimated IFs is significantly improved compared with the initial IFs
%    RE2 =  norm(IFmset(3,:,end)-IF2)/norm(IF2)
%    RE3 =  norm(IFmset(4,:,end)-IF3)/norm(IF3)
% %% Reconstructed modes
% figure
% set(gcf,'Position',[20 100 640 200]);	    
% set(gcf,'Color','w'); 
% plot(t,smset(1,:,end),'linewidth',2)	% estimated mode
% hold on
% plot(t,Sig4 - smset(1,:,end),'k','linewidth',2) % estimation errors
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('m1','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
% set(gca,'FontSize',24);
% set(gca,'linewidth',2);	
% axis([0 6 -0.5 6])
% 
% figure
% set(gcf,'Position',[20 100 640 200]);	    
% set(gcf,'Color','w'); 
% plot(t,smset(2,:,end),'linewidth',2)	% estimated mode
% hold on
% plot(t,Sig1 - smset(2,:,end),'k','linewidth',2)  % estimation errors
% hold on
% plot(t,IA(2,:),'r','linewidth',3) % estimated IAs
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('m2','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
% set(gca,'FontSize',24);
% set(gca,'linewidth',2);	
% axis([0 6 -4 4])
% 
% figure
% set(gcf,'Position',[20 100 640 200]);	    
% set(gcf,'Color','w'); 
% plot(t,smset(3,:,end),'linewidth',2)	% estimated mode
% hold on
% plot(t,Sig2 - smset(3,:,end),'k','linewidth',2)  % estimation errors
% hold on
% plot(t,IA(3,:),'r','linewidth',3) % estimated IAs
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('m3','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
% set(gca,'FontSize',24);
% set(gca,'linewidth',2);	
% axis([0 6 -3 3])
% 
% figure
% set(gcf,'Position',[20 100 640 200]);	    
% set(gcf,'Color','w'); 
% plot(t,smset(4,:,end),'linewidth',2)	% estimated mode
% hold on
% plot(t,Sig3 - smset(4,:,end),'k','linewidth',2)  % estimation errors
% hold on
% plot(t,IA(4,:),'r','linewidth',3) % estimated IAs
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('m3','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
% set(gca,'FontSize',24);
% set(gca,'linewidth',2);	
% axis([0 6 -2 2])
