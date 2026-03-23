%%%%%%%%% the signal components are crossing %%%%%%%%%%%%%%
% this example is adopted from paper£ºChen S, Peng Z, Yang Y, et al, Intrinsic chirp component decomposition by using Fourier Series representation, Signal Processing, 2017.
clc
clear
close all
SampFreq = 100;
t = 0:1/SampFreq:15;
Sig1 = cos(2*pi*(0.23+15*t + 0.2*t.^2));
IF1 = 15 + 0.4*t;
amp2 = 0.5*cos(2*pi*0.3*t)+1;
Sig2 = amp2.*cos(2*pi*(5*sin(pi/4*t) +5*t + 1.2*t.^2 ));
IF2 = 5*pi/4*cos(pi*t/4) + 2.4*t + 5;
Sig3 = cos(2*pi*(0.35+35*t - 0.8*t.^2));
IF3 = 35 - 1.6*t;

Sig = Sig1 + Sig2 + Sig3;

figure
set(gcf,'Position',[20 100 320 250]);	    
set(gcf,'Color','w'); 
plot(t,Sig);
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',12)
set(gca,'linewidth',1);
%% STFT
window = 128;
Nfrebin = 1024;
figure
[Spec,f] = STFT(Sig',SampFreq,Nfrebin,window);
imagesc(t,f,abs(Spec)); 
axis([0 15 0 50]);
set(gcf,'Position',[20 100 320 250]);	 
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',12);
set(gcf,'Color','w');	
%% ridge extraction and smoothing
bw = SampFreq/80;% the bandwidth of the TF filter for ridge extraction
beta1 = 1e-4; % beta1 should be larger than the following beta
num = 3; % the number of the components
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
axis([0 15 0 50]);
%% ridge path regrouping (RPRG)
% the RPRG algorithm is developed to extract ridge curves of crossed signal modes
% more details about the RPRG can be found in paper: Chen S, Dong X, Xing G, et al, Separation of Overlapped Non-Stationary Signals by Ridge Path Regrouping and Intrinsic Chirp Component Decomposition, IEEE Sensors Journal, 2017.
thrf = length(f)/30;

[findex,interset] = RPRG(fidexmult,thrf);

figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,f(findex),'linewidth',3);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 15 0 50]);
%% parameter setting
alpha = 1e-5; 
beta = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
iniIF = curvesmooth(f(findex),beta); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges
var = 0;% noise variance
tol = 1e-8;%
tic
[IFmset IA smset] = VNCMD(Sig,SampFreq,iniIF,alpha,beta,var,tol);
toc

figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,iniIF(:,:),'linewidth',3); % initial IFs
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
axis([0 15 0 50]);

%% Relative errors of the initial IFs
   iniRE1 =  norm(iniIF(1,:)-IF1)/norm(IF1) 
   iniRE2 =  norm(iniIF(2,:)-IF2)/norm(IF2)
   iniRE3 =  norm(iniIF(3,:)-IF3)/norm(IF3)
%% estimated IF
figure
plot(t,[IF1;IF2;IF3],'b','linewidth',3) % true IFs
hold on
plot(t,IFmset(:,:,end),'r','linewidth',3) % finally estimated IFs
set(gcf,'Position',[20 100 640 500]);	 
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');	
axis([0 15 0 50])
%% Relative errors of the finally estimated IFs
   RE1 =  norm(IFmset(1,:,end)-IF1)/norm(IF1) % one can find that the accuracy of the finally estimated IFs is significantly improved compared with the initial IFs
   RE2 =  norm(IFmset(2,:,end)-IF2)/norm(IF2)
   RE3 =  norm(IFmset(3,:,end)-IF3)/norm(IF3)
%% Reconstructed modes

figure
set(gcf,'Position',[20 100 640 200]);	    
set(gcf,'Color','w'); 
plot(t,smset(1,:,end),'linewidth',2)	% estimated mode
hold on
plot(t,Sig1 - smset(1,:,end),'k','linewidth',2)  % estimation errors
hold on
plot(t,IA(1,:),'r','linewidth',3) % estimated IAs
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('m1','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);	
axis([0 15 -1.5 1.5])

figure
set(gcf,'Position',[20 100 640 200]);	    
set(gcf,'Color','w'); 
plot(t,smset(2,:,end),'linewidth',2)	% estimated mode
hold on
plot(t,Sig2 - smset(2,:,end),'k','linewidth',2)  % estimation errors
hold on
plot(t,IA(2,:),'r','linewidth',3) % estimated IAs
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('m2','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);	
axis([0 15 -2 2])

figure
set(gcf,'Position',[20 100 640 200]);	    
set(gcf,'Color','w'); 
plot(t,smset(3,:,end),'linewidth',2)	% estimated mode
hold on
plot(t,Sig3 - smset(3,:,end),'k','linewidth',2)  % estimation errors
hold on
plot(t,IA(3,:),'r','linewidth',3) % estimated IAs
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('m3','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);	
axis([0 15 -1.5 1.5])
