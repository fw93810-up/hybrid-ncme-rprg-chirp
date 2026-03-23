%% test_hybrid_RPRG_VNCMD_NCME.m
% 测试 hybrid_RPRG_VNCMD_NCME 在 3 模态（含交叉/非线性调频 + AM）信号上的表现

clear; clc; close all;

%% 0) 生成测试信号（你给的）
SampFreq = 100;
t = 0:1/SampFreq:15;            % 0~15s, 步长 0.01s
Sig1 = cos(2*pi*(0.23 + 15*t + 0.2*t.^2));
amp2 = 0.5*cos(2*pi*0.3*t) + 1;
Sig2 = amp2 .* cos(2*pi*(5*sin(pi/4*t) + 5*t + 1.2*t.^2));
Sig3 = cos(2*pi*(0.35 + 35*t - 0.8*t.^2));
Sig  = Sig1 + Sig2 + Sig3;

Sig = Sig(:)';                 % 确保行向量
N   = length(Sig);

%% 1) 理论真值 IF（由 phase(t) 的导数得到，单位 Hz）
% 注意：你的 cos(2*pi*(...)) 里括号内是"cycles"相位；对 t 求导就是 IF (Hz)
IF1_true = 15 + 0.4*t;                                  % d/dt (15t + 0.2t^2)
IF3_true = 35 - 1.6*t;                                  % d/dt (35t - 0.8t^2)
IF2_true = 5*(pi/4)*cos(pi/4*t) + 5 + 2.4*t;            % d/dt (5*sin(pi/4 t) + 5t + 1.2t^2)

%% 2) 配置 cfg（按你的 hybrid 函数头注释要求填全）
cfg = struct();

% 强/总模态数
cfg.num_strong = 3;    % 强模态数（建议 2 或 3）
cfg.num_total  = 3;    % 总模态数（这里就是 3）

% STFT 参数
cfg.STFT.winLen  = 256;
cfg.STFT.Nfrebin = 1024;

% Ridge 参数（extridge_mult / findridges / Dechirp_filter）
cfg.Ridge.delta      = 20;                  % 允许的 bin 跳变
cfg.Ridge.beta_ridge = 1e-4;                % ridge 平滑
cfg.Ridge.bw         = SampFreq/80;         % Dechirp_filter 带宽（Hz）

% RPRG 阈值：thrf = length(f) / thrf_scale  （越大 => thrf越小）
% 想要 thrf ~ 2 bins：length(f)=1024 => thrf_scale=512
cfg.RPRG.thrf_scale = 30;

% VNCMD（可选）
cfg.VNCMD = struct();
cfg.VNCMD.use   = true;      % 先关掉，确保流程可跑；想开就改 true
cfg.VNCMD.alpha = 1e-5;       % 下面这些参数需要你按你的 VNCMD 实现调
cfg.VNCMD.beta  = 1e-5;
cfg.VNCMD.var   = 0;
cfg.VNCMD.tol   = 1e-8;

% NCME 参数
cfg.NCME = struct();
cfg.NCME.lambda = 1e6;        % 稀疏正则（你之前脚本常用 1e4）
cfg.NCME.beta   = 1e-9;       % IF 增量平滑
cfg.NCME.tol    = 1e-6;

%% 3) 跑 hybrid
fprintf('\n===== Running hybrid_RPRG_VNCMD_NCME =====\n');
out = hybrid_RPRG_VNCMD_NCME(Sig, SampFreq, cfg);

%% 4) 画图：STFT + 估计 IF + 真值 IF（对比）
figure('Color','w');
imagesc(out.t, out.f, out.Spectrogram);
set(gca,'YDir','normal'); axis tight;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('STFT |S| with estimated IFs (NCME) and true IFs');
colorbar; hold on;

% 估计 IF（NCME 输出）
K = out.num_total;
clr = lines(K);
for k = 1:K
    plot(out.t, out.IF_ncme(k,:), 'LineWidth', 1.6, 'Color', clr(k,:));
end

% 真值 IF（虚线黑）
plot(t, IF1_true, 'k--', 'LineWidth', 1.2);
plot(t, IF2_true, 'k--', 'LineWidth', 1.2);
plot(t, IF3_true, 'k--', 'LineWidth', 1.2);

legend([arrayfun(@(k)sprintf('IF_{est} mode %d',k),1:K,'UniformOutput',false), ...
        {'IF_{true}'}], 'Location','best');

%% 5) 画图：各模态时域 + 残差
figure('Color','w');
K = out.num_total;

for k = 1:K
    subplot(K+1,1,k);
    plot(out.t, out.modes(k,:), 'LineWidth', 1.1);
    grid on;
    ylabel(sprintf('mode %d',k));
    if k==1
        title('Estimated modes (time domain)');
    end
end

subplot(K+1,1,K+1);
plot(out.t, out.residual, 'LineWidth', 1.1);
grid on;
xlabel('Time (s)'); ylabel('residual');
title('Residual = Sig - sum(modes)');

%% 6) 可选：打印一些简单指标
recon = sum(out.modes(1:K,:), 1);
rel_err = norm(Sig - recon) / norm(Sig);
fprintf('\nRelative reconstruction error ||Sig - sum(modes)|| / ||Sig|| = %.3e\n', rel_err);

%% 7) 如果你想看"强模态初值 IF（RPRG(+VNCMD)）vs NCME 最终 IF"
figure('Color','w');
imagesc(out.t, out.f, out.Spectrogram);
set(gca,'YDir','normal'); axis tight;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Strong IF init (RPRG) vs final IF (NCME)');
colorbar; hold on;

Ks = out.num_strong;

for k = 1:K
    plot(out.t, out.IF_ncme(k,:), 'LineWidth', 1.6, 'Color', clr(k,:)); % 最终
end
legend({'IF_{init,strong} (RPRG)', 'IF_{final} (NCME)'}, 'Location','best');


