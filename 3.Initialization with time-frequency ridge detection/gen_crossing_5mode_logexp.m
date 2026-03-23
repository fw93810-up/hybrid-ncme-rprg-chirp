%% gen_crossing5_modes_poly_cubic_noise.m
% 5 分量多项式 IF/IA 信号 + t^3 项 + 可调高斯白噪声
% - 2 条强分量：多项式 chirp，频率在中频段，相互交叉，IF 里有 t^2 / t^3
% - 3 条弱分量：多项式 IF，不交叉，能量依次减弱，同样带 t^2 / t^3

clc; clear; close all;

%% 1. 时间轴 & 采样率
fs = 400;              % 采样频率
T  = 1.0;              % 信号总时长 1 s
t  = 0 : 1/fs : T;     % 时间向量
N  = length(t);

%% 2. 设计 5 个分量的 IF（带 t^3） & IA（也带高阶）

% -------- 强 1：上升 chirp，多项式 IF（凸一点） --------
% IF1(t) ≈ 40 → 180 Hz
IF1 = 40  + 70*t + 40*t.^2 + 30*t.^3;

% -------- 强 2：下降 chirp，多项式 IF（有点弯），和 IF1 交叉 --------
% IF2(t) ≈ 160 → 10 Hz，和 IF1 在中间附近交叉
IF2 = 160 - 80*t - 60*t.^2 - 20*t.^3;

% -------- 弱 3：低频、多项式 IF，不交叉（一直很低） --------
% 大致 5 → 22 Hz
IF3 = 5   + 10*t +  5*t.^2 +  2*t.^3;

% -------- 弱 4：中高频、多项式 IF，不交叉 --------
% 大致 140 → 180 Hz
IF4 = 140 + 20*t + 10*t.^2 + 10*t.^3;

% -------- 弱 5：更高频，多项式 IF，不交叉 --------
% 大致 180 → 200 Hz（控制在 Nyquist=200 之内）
IF5 = 180 + 10*t +  5*t.^2 +  5*t.^3;

% -------- 对应的 IA（幅值），也用多项式 + t^2/t^3，且逐条变弱 --------
IA1 = 1.6 * (1 + 0.4*t + 0.2*t.^2);                % 最强
IA2 = 1.3 * (1 + 0.3*t + 0.1*t.^2 + 0.1*t.^3);     % 次强
IA3 = 0.9 * (1 + 0.3*t + 0.2*t.^2 + 0.1*t.^3);     % 弱 1
IA4 = 0.6 * (1 + 0.2*t + 0.1*t.^2 + 0.05*t.^3);    % 弱 2
IA5 = 0.4 * (1 + 0.15*t + 0.1*t.^2 + 0.05*t.^3);   % 最弱

%% 3. IF 积分得到相位，再生成分量信号
phi1 = 2*pi*cumtrapz(t, IF1);
phi2 = 2*pi*cumtrapz(t, IF2);
phi3 = 2*pi*cumtrapz(t, IF3);
phi4 = 2*pi*cumtrapz(t, IF4);
phi5 = 2*pi*cumtrapz(t, IF5);

Sig1 = IA1 .* cos(phi1);
Sig2 = IA2 .* cos(phi2);
Sig3 = IA3 .* cos(phi3);
Sig4 = IA4 .* cos(phi4);
Sig5 = IA5 .* cos(phi5);

Sig_clean = Sig1 + Sig2 + Sig3 + Sig4 + Sig5;

%% 4. 加噪声：高斯白噪声，按目标 SNR 设计
SNR_dB = 10;   % 你可以改成 20/10/5/0 等不同 SNR

P_sig   = mean(Sig_clean.^2);            % 信号平均功率
P_noise = P_sig / 10^(SNR_dB/10);        % 噪声目标功率
sigma_n = sqrt(P_noise);                 % 噪声标准差

noise   = sigma_n * randn(size(Sig_clean));
Sig     = Sig_clean + noise;             % 含噪信号（之后算法都用这个）

%% 5. 时域：总信号 + 各分量
figure;
set(gcf,'Color','w','Position',[100 100 900 650]);

subplot(6,1,1);
plot(t, Sig, 'k', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amp');
title(sprintf('Total noisy signal, SNR = %.1f dB', SNR_dB));
grid on;

subplot(6,1,2);
plot(t, Sig1, 'r', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 1 (strong up-chirp, poly+ t^3)');
grid on;

subplot(6,1,3);
plot(t, Sig2, 'b', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 2 (strong down-chirp, crossing, poly+ t^3)');
grid on;

subplot(6,1,4);
plot(t, Sig3, 'm', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 3 (weak low-freq poly+ t^3)');
grid on;

subplot(6,1,5);
plot(t, Sig4, 'g', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 4 (weak mid-high poly+ t^3)');
grid on;

subplot(6,1,6);
plot(t, Sig5, 'c', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amp'); title('Component 5 (weakest high poly+ t^3)');
grid on;

%% 6. 真 IF 曲线（方便后面对比 NCME / VNCMD）
figure;
set(gcf,'Color','w','Position',[200 120 800 400]);

plot(t, IF1, 'r', 'LineWidth', 1.6); hold on;
plot(t, IF2, 'b', 'LineWidth', 1.6);
plot(t, IF3, 'm--', 'LineWidth', 1.2);
plot(t, IF4, 'g--', 'LineWidth', 1.2);
plot(t, IF5, 'c--', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('True IFs of 5 polynomial (with t^3) components');
legend('IF1 (strong up-chirp)', ...
       'IF2 (strong down-chirp)', ...
       'IF3 (weak low)', ...
       'IF4 (weak mid-high)', ...
       'IF5 (weak high)', ...
       'Location','best');
ylim([0 220]);
grid on;

%% 7. STFT + 真 IF 叠加（看交叉 & 噪声下 TFR）
window  = 128;
Nfrebin = 512;

[Spec, f_stft] = STFT(Sig.', fs, Nfrebin, window);

figure;
set(gcf,'Color','w','Position',[180 150 900 400]);
imagesc(t, f_stft, abs(Spec));
set(gca,'YDir','normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('STFT magnitude (noisy poly t^3 signal, SNR = %.1f dB)', SNR_dB));
axis([0 T 0 220]);
colorbar; hold on;

plot(t, IF1, 'w-',  'LineWidth', 1.4);
plot(t, IF2, 'w-',  'LineWidth', 1.4);
plot(t, IF3, 'w--', 'LineWidth', 1.0);
plot(t, IF4, 'w--', 'LineWidth', 1.0);
plot(t, IF5, 'w--', 'LineWidth', 1.0);

%% 8. 保存 mat，方便后续用 hybrid_RPRG_VNCMD_NCME / NCME_multi 等分析
save('crossing5_modes_poly_cubic_noisy.mat', ...
     'fs','t', ...
     'Sig','Sig_clean', ...
     'Sig1','Sig2','Sig3','Sig4','Sig5', ...
     'IF1','IF2','IF3','IF4','IF5', ...
     'IA1','IA2','IA3','IA4','IA5', ...
     'SNR_dB','noise');

disp('Done. Saved as crossing5_modes_poly_cubic_noisy.mat');
