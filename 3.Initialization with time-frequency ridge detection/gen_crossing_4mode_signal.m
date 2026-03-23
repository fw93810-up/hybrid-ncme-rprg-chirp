%% gen_crossing_4mode_signal_short.m
% 生成一个「短一点」的 4 分量交叉信号：
%   fs = 400 Hz, T = 1 s, N ≈ 401
%   其中两条强分量在 t=0.5 s 左右交叉，频率 80 Hz 附近

clc;
clear;
close all;

%% 1. 时间轴 & 采样率（缩短版）
fs = 400;                 % 采样频率从 1000 降到 400 Hz
T  = 1;                   % 总时长从 2 s 缩到 1 s
t  = 0 : 1/fs : T;        % 时间向量
N  = length(t);

%% 2. 设计 4 个分量的 IF & IA
% 注意：因为 T 变成 1s，我们适当调一下斜率，让交叉仍然发生在中间

% === 强分量 1：从 40 Hz 线性扫到 120 Hz（上升 chirp）===
%  IF1(t) = 40 + 80*t  (t∈[0,1], 40→120 Hz)
IF1 = 40 + 80*t;
IA1 = 1.5 * (1 + 0.2*sin(2*pi*0.3*t));   % 轻微调制振幅

% === 强分量 2：从 120 Hz 线性扫到 40 Hz（下降 chirp）===
%  IF2(t) = 120 - 80*t (t∈[0,1], 120→40 Hz)
%  两条在 t=0.5 s, f=80 Hz 交叉
IF2 = 120 - 80*t;
IA2 = 1.3 * (1 + 0.15*cos(2*pi*0.25*t));

% === 弱分量 3：低频，非交叉 ===
IF3 = 20 + 5*sin(2*pi*0.5*t);           % 大致在 15~25 Hz 之间
IA3 = 0.6 * (1 + 0.1*sin(2*pi*0.4*t));

% === 弱分量 4：高频，非交叉 ===
% Nyquist(fs=400) = 200 Hz，所以 180±10Hz 是安全的
IF4 = 180 + 10*sin(2*pi*0.4*t);         % 大致在 170~190 Hz
IA4 = 0.4 * (1 + 0.2*cos(2*pi*0.35*t));

%% 3. 由 IF 积分得到相位，再生成 4 个分量信号
phi1 = 2*pi*cumtrapz(t, IF1);
phi2 = 2*pi*cumtrapz(t, IF2);
phi3 = 2*pi*cumtrapz(t, IF3);
phi4 = 2*pi*cumtrapz(t, IF4);

Sig1 = IA1 .* cos(phi1);
Sig2 = IA2 .* cos(phi2);
Sig3 = IA3 .* cos(phi3);
Sig4 = IA4 .* cos(phi4);

Sig = Sig1 + Sig2 + Sig3 + Sig4;

%% 4. 时域：总信号 + 各分量（简单看一眼）

figure;
set(gcf,'Color','w','Position',[100 100 800 600]);

subplot(5,1,1);
plot(t, Sig, 'k', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amp');
title('Total signal (short, 4 components, 2 strong & crossing)');
grid on;

subplot(5,1,2);
plot(t, Sig1, 'r', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 1 (strong up-chirp)');
grid on;

subplot(5,1,3);
plot(t, Sig2, 'b', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 2 (strong down-chirp, crossing)');
grid on;

subplot(5,1,4);
plot(t, Sig3, 'm', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 3 (weak low-freq)');
grid on;

subplot(5,1,5);
plot(t, Sig4, 'g', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 4 (weak high-freq)');
grid on;

%% 5. 真 IF 曲线

figure;
set(gcf,'Color','w','Position',[200 100 800 400]);

plot(t, IF1, 'r', 'LineWidth', 1.5); hold on;
plot(t, IF2, 'b', 'LineWidth', 1.5);
plot(t, IF3, 'm--', 'LineWidth', 1.2);
plot(t, IF4, 'g--', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('True IFs of 4 components (1 & 2 strong & crossing, short signal)');
legend('IF1 (strong up-chirp)', ...
       'IF2 (strong down-chirp)', ...
       'IF3 (weak low)', ...
       'IF4 (weak high)', ...
       'Location','best');
ylim([0 220]);
grid on;

%% 6. STFT + 真 IF 叠加

window  = 128;        % 可以也相应减小一点
Nfrebin = 512;

[Spec, f_stft] = STFT(Sig.', fs, Nfrebin, window);

figure;
set(gcf,'Color','w','Position',[150 150 800 400]);
imagesc(t, f_stft, abs(Spec));
set(gca,'YDir','normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('STFT magnitude (short signal) with true IF curves');
axis([0 T 0 220]);
colorbar; hold on;

plot(t, IF1, 'w-',  'LineWidth', 1.5);
plot(t, IF2, 'w-',  'LineWidth', 1.5);
plot(t, IF3, 'w--', 'LineWidth', 1.2);
plot(t, IF4, 'w--', 'LineWidth', 1.2);

%% 7. 存成 mat，方便后面直接跑 ncme_multimode_rprg

save('crossing4_modes_short.mat', ...
     'fs','t', ...
     'Sig','Sig1','Sig2','Sig3','Sig4', ...
     'IF1','IF2','IF3','IF4', ...
     'IA1','IA2','IA3','IA4');

disp('Done. Saved as crossing4_modes_short.mat');
