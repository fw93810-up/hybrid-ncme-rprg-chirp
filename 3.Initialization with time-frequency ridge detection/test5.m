clc; clear; close all;

%% 1. 载入信号 & 配置参数
load('crossing4_modes_short.mat');      % 假设变量名为 data
Sig = data(:)';        % 行向量
fs  = 250;             % 采样频率，和之前保持一致

cfg = struct;

% 先认为有两条强线用 RPRG 处理，最终想要 4 个模态
cfg.num_strong = 2;
cfg.num_total  = 4;

% STFT 参数
cfg.STFT.winLen  = 256;
cfg.STFT.Nfrebin = 1024;

% Ridge / Dechirp 参数
cfg.Ridge.delta      = 20;
cfg.Ridge.beta_ridge = 1e-4;
cfg.Ridge.bw         = fs/80;

% RPRG 参数
cfg.RPRG.thrf_scale = 30;    % thrf = length(f)/30

% VNCMD（仅用于强模态 IF 的微调）
cfg.VNCMD.use   = true;
cfg.VNCMD.alpha = 1e-5;
cfg.VNCMD.beta  = 1e-5;
cfg.VNCMD.var   = 0;
cfg.VNCMD.tol   = 1e-8;

% NCME_multi 参数
cfg.NCME.lambda = 1e3;
cfg.NCME.beta   = 1e-6;
cfg.NCME.tol    = 1e-8;

% 调用混合算法
out = hybrid_RPRG_VNCMD_NCME(Sig, fs, cfg);

% 取出结果
t        = out.t;
f        = out.f;
SpecMag  = abs(out.Spectrogram);
IF_rprg  = out.IF_rprg;      % 强模态 RPRG(+VNCMD) 的初始 IF（K_strong×N）
IF_ncme  = out.IF_ncme;      % 最终所有模态的 NCME_multi IF（num_total×N）
modes    = out.modes;        % 各模态时域波形（num_total×N）
IA       = out.IA;           % IA 估计（num_total×N）
residual = out.residual;     % 最终残差
num_strong = out.num_strong;
num_total  = out.num_total;

%% 2. 原始时域波形
figure;
set(gcf,'Position',[20 100 640 250],'Color','w');
plot(t, Sig, 'w', 'LineWidth', 1.5);   % 白色线，适合深色背景
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
title('Original signal');
set(gca,'FontSize',12,'LineWidth',1);
grid on;

%% 3. STFT + 脊线 & NCME IF 叠加
figure;
imagesc(t, f, SpecMag);
set(gca,'YDir','normal','FontSize',12,'LineWidth',1);
set(gcf,'Color','w');
axis([t(1) t(end) 0 fs/2]);
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
title('STFT with RPRG IFs (white dashed) and NCME IFs (colored)');
colorbar;
hold on;

clr = lines(num_total);   % 每条模态一条彩色线

% 1) RPRG (+VNCMD) 得到的强模态 IF：白虚线
if ~isempty(IF_rprg)
    plot(t, IF_rprg.', 'w--', 'LineWidth', 1.2);
end

% 2) NCME_multi 最终 IF：彩色实线
for k = 1:num_total
    plot(t, IF_ncme(k,:), 'Color', clr(k,:), 'LineWidth', 1.6);
end

%% 4. 每个模态单独在 T-F 图上叠加 IF（可视化每条线）
figure;
set(gcf,'Position',[100 100 900 600],'Color','w');

for k = 1:num_total
    if num_total == 4
        subplot(2,2,k);
    else
        subplot(num_total,1,k);
    end

    imagesc(t, f, SpecMag);
    set(gca,'YDir','normal','FontSize',10,'LineWidth',1);
    axis([t(1) t(end) 0 fs/2]);
    hold on;

    % 对应的 RPRG IF（如果这个索引在强模态范围内）
    if k <= num_strong
        plot(t, IF_rprg(k,:), 'w--', 'LineWidth', 1.2);
    end

    % 对应的 NCME IF
    plot(t, IF_ncme(k,:), 'Color', clr(k,:), 'LineWidth', 1.6);

    xlabel('Time / Sec');
    ylabel('Frequency / Hz');
    title(sprintf('Mode %d: IF (white = RPRG, color = NCME)', k));
end

%% 5. 各模态时域波形（类似之前的 data.mat 版 Mode1~Mode4）
figure;
set(gcf,'Position',[100 100 900 600],'Color','w');

for k = 1:num_total
    if num_total == 4
        subplot(2,2,k);
    else
        subplot(num_total,1,k);
    end

    plot(t, modes(k,:), 'Color', clr(k,:), 'LineWidth', 1.4); 
    hold on;
    % 背景上用淡灰色画一下原始信号
    plot(t, Sig, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8);

    xlabel('Time / Sec');
    ylabel('Amp');
    title(sprintf('Mode %d (estimated)', k));
    set(gca,'FontSize',12,'LineWidth',1);
    grid on;
end

%% 6. Sum of all modes vs Original & Error
sum_modes = sum(modes, 1);
err_sig   = Sig - sum_modes;

figure;
set(gcf,'Position',[50 50 800 500],'Color','w');

subplot(2,1,1);
plot(t, sum_modes, 'r', 'LineWidth', 1.4); 
hold on;
plot(t, Sig, 'c--', 'LineWidth', 1.0);   % 原始信号用青色虚线
xlabel('Time (s)');
ylabel('Amp');
legend('Sum of all modes','Original signal','Location','best');
title('(a) Sum of all estimated modes vs original');
set(gca,'FontSize',12,'LineWidth',1);
grid on;

subplot(2,1,2);
plot(t, err_sig, 'b', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Amp');
title('(b) Error = Original - Sum of modes');
set(gca,'FontSize',12,'LineWidth',1);
grid on;

%% 7. 残差信号本身（如果你想单独看）
figure;
set(gcf,'Position',[50 100 640 250],'Color','w');
plot(t, residual, 'm', 'LineWidth', 1.2);
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Residual','FontSize',12,'FontName','Times New Roman');
title('Final residual signal');
set(gca,'FontSize',12,'LineWidth',1);
grid on;
