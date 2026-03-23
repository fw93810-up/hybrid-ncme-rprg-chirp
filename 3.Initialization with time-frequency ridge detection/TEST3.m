%% TEST_crossing4_hybrid.m
% 用 hybrid_RPRG_VNCMD_NCME 处理 crossing4_modes_short.mat
% 依赖：
%   hybrid_RPRG_VNCMD_NCME.m
%   STFT.m, extridge_mult.m, Dechirp_filter.m, findridges.m,
%   RPRG.m, curvesmooth.m, VNCMD.m, NCME_multi.m （以及它们依赖的 gene_lasso 等）

clc; clear; close all;

%% 1. 载入交叉信号（你之前生成的 crossing4_modes_short.mat）
load('crossing5_modes_poly_noisy.mat');   % 里面有 fs, t, Sig, Sig1..Sig4, IF1..IF4, IA1..IA4

Sig_in = Sig(:)';    % 保证是行向量
fs     = fs;         % 采样频率就用文件里的

N = length(Sig_in);
t = (0:N-1)/fs;      % 跟 hybrid 函数内部一致的时间轴（也可以用文件里的 t，看上去差不多）

%% 2. 配置 cfg

cfg = struct;

% 先验：有 2 条强交叉分量（IF1、IF2），共 4 个模态
cfg.num_strong = 2;      % 交叉的两条强线
cfg.num_total  = 4;      % 总共希望分成 4 条

% STFT 参数（和生成脚本里差不多）
cfg.STFT.winLen  = 128;
cfg.STFT.Nfrebin = 512;

% Ridge / Dechirp 参数
cfg.Ridge.delta      = 20;        % 相邻时间允许的频率 bin 跳变
cfg.Ridge.beta_ridge = 1e-4;      % Dechirp_filter 平滑参数
cfg.Ridge.bw         = fs / 80;   % TF 滤波带宽（Hz）

% RPRG 参数：thrf = length(f) / thrf_scale
cfg.RPRG.thrf_scale = 30;

% VNCMD：只用来微调强分量 IF，不直接用它的重构
cfg.VNCMD.use   = true;
cfg.VNCMD.alpha = 1e-5;
cfg.VNCMD.beta  = 1e-5;
cfg.VNCMD.var   = 0;
cfg.VNCMD.tol   = 1e-8;

% NCME_multi 参数（你可以在 10~1000 之间调 lambda）
cfg.NCME.lambda = 1e3;
cfg.NCME.beta   = 1e-6;
cfg.NCME.tol    = 1e-8;

%% 3. 跑 hybrid 管线

out = hybrid_RPRG_VNCMD_NCME(Sig_in, fs, cfg);

% 取出结果
t_h   = out.t;                     % 1×N
f_h   = out.f;                     % Nf×1
Spec  = abs(out.Spectrogram);     % Nf×N

IF_rprg = out.IF_rprg;            % num_strong×N  (RPRG(+VNCMD) 强分量 IF)
IF_ncme = out.IF_ncme;            % num_total×N   (最终 NCME_multi IF)
modes   = out.modes;              % num_total×N   (各模态时域信号)
IA      = out.IA;                 % num_total×N   (各模态 IA)
resid   = out.residual;           % 1×N           (最终残差)

num_strong = out.num_strong;
num_total  = out.num_total;

% 真 IF（来自 crossing4_modes_short.mat）
IF_true = [IF1; IF2; IF3; IF4];   % 4×N

% 总重构与误差
Sig_rec = sum(modes, 1);
Sig_err = Sig_in - Sig_rec;

colors = lines(num_total);

%% 4. 原始 STFT + RPRG 强 IF + NCME IF（估计的）

figure;
set(gcf,'Color','w','Position',[100 100 900 450]);

imagesc(t_h, f_h, Spec);
set(gca,'YDir','normal');
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',12,'FontName','Times New Roman');
axis([t_h(1) t_h(end) 0 fs/2]);
set(gca,'FontSize',12,'LineWidth',1);
colormap(jet);
colorbar;
title('STFT with RPRG(+VNCMD) IF (white) and NCME-multi IF (colors)');
hold on;

% 强分量 RPRG(+VNCMD) 初始 IF：白虚线
if ~isempty(IF_rprg)
    plot(t_h, IF_rprg', 'w--', 'LineWidth', 1.2);
end

% NCME_multi 最终 IF：彩色实线
for k = 1:num_total
    plot(t_h, IF_ncme(k,:), 'Color', colors(k,:), 'LineWidth', 1.6);
end

%% 5. 真 IF vs 估计 IF（在同一个 TFR 上）

figure;
set(gcf,'Color','w','Position',[100 100 900 450]);

imagesc(t_h, f_h, Spec);
set(gca,'YDir','normal');
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',12,'FontName','Times New Roman');
axis([t_h(1) t_h(end) 0 fs/2]);
set(gca,'FontSize',12,'LineWidth',1);
colormap(jet);
colorbar;
title('STFT with true IFs (white) and NCME-multi IFs (colors)');
hold on;

% 真 IF：白色线（两条强分量实线，两个弱分量虚线）
plot(t, IF1, 'w-',  'LineWidth', 1.4);
plot(t, IF2, 'w-',  'LineWidth', 1.4);
plot(t, IF3, 'w--', 'LineWidth', 1.2);
plot(t, IF4, 'w--', 'LineWidth', 1.2);

% 估计 IF：彩色
for k = 1:num_total
    plot(t_h, IF_ncme(k,:), 'Color', colors(k,:), 'LineWidth', 1.6);
end

%% 6. 每个模态在 T–F 上单独看（强模态还附带 RPRG IF）

figure;
set(gcf,'Color','w','Position',[100 100 900 600]);

for k = 1:num_total
    if num_total == 4
        subplot(2,2,k);
    else
        subplot(num_total,1,k);
    end

    imagesc(t_h, f_h, Spec);
    set(gca,'YDir','normal');
    axis([t_h(1) t_h(end) 0 fs/2]);
    colormap(jet);
    hold on;

    % 对前 num_strong 条，加 RPRG(+VNCMD) 的初始 IF
    if k <= num_strong && ~isempty(IF_rprg)
        plot(t_h, IF_rprg(k,:), 'w--', 'LineWidth', 1.2);
    end

    % 该模态的 NCME IF
    plot(t_h, IF_ncme(k,:), 'Color', colors(k,:), 'LineWidth', 1.6);

    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    if k <= num_strong && ~isempty(IF_rprg)
        title(sprintf('Mode %d: RPRG IF (white) vs NCME IF (color)', k));
    else
        title(sprintf('Mode %d: NCME IF', k));
    end
    set(gca,'FontSize',10,'LineWidth',1);
end

%% 7. 总重构 vs 原始 + 残差

figure;
set(gcf,'Color','w','Position',[80 80 800 500]);

subplot(2,1,1);
plot(t, Sig_rec, 'r', 'LineWidth', 1.2); hold on;
plot(t, Sig_in, 'c--', 'LineWidth', 0.8);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Sum of all modes','Original signal','Location','best');
title('(a) Sum of all estimated modes vs original');
set(gca,'FontSize',12,'LineWidth',1);
grid on;

subplot(2,1,2);
plot(t, Sig_err, 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('(b) Error = Original - Sum of modes');
set(gca,'FontSize',12,'LineWidth',1);
grid on;

%% 8. 每个模态的时域波形（可和真分量对比）

figure;
set(gcf,'Color','w','Position',[120 120 900 600]);

true_modes = [Sig1; Sig2; Sig3; Sig4];   % 4×N

for k = 1:num_total
    if num_total == 4
        subplot(2,2,k);
    else
        subplot(num_total,1,k);
    end

    plot(t, modes(k,:), 'Color', colors(k,:), 'LineWidth', 1.3); hold on;

    % 把第 k 条真分量也画上去（注意：可能存在"模式编号对不上"的情况，
    % 这里只是直接按 1→1,2→2,3→3,4→4 画，具体对应关系你可以看 IF 形状来判断）
    if k <= size(true_modes,1)
        plot(t, true_modes(k,:), 'k--', 'LineWidth', 1.0);
        legend('Estimated mode','True component','Location','best');
    else
        legend('Estimated mode','Location','best');
    end

    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Mode %d (estimated vs true)', k));
    set(gca,'FontSize',12,'LineWidth',1);
    grid on;
end
