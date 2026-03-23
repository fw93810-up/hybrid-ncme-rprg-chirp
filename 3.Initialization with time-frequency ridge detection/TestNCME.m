%%%%%%%%%%%%  Test script for real test signal in data.mat + NCME  %%%%%%%%%%%%
% 需要: STFT.m, extridge_mult.m, Dechirp_filter.m, findridges.m, NCME.m
%   （如果没有 Signal Processing Toolbox，请在 STFT / extridge_mult 里用 myhilbert 替代 hilbert）

clc;
clear;
close all;

%% 1. 载入测试信号
load('data.mat');          % 假设里面的变量名是 data
Sig = data(:)';            % 转成行向量，保证 extridge_mult / NCME 接受

% 采样频率：这里先用 250（只影响坐标刻度，不影响算法内部）
SampFreq = 250;

N = length(Sig);
t = (0:N-1) / SampFreq;    % 时间轴（秒）

%% 2. 画原始时域波形
figure;
set(gcf,'Position',[20 100 640 250]);
set(gcf,'Color','w');
plot(t, Sig, 'w', 'LineWidth', 1.5);
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
grid on;
title('Original signal (time domain)');

%% 3. STFT 时频图
window  = 256;             % 窗长
Nfrebin = 1024;            % 频率点数（FFT 长度）

[Spec, f] = STFT(Sig', SampFreq, Nfrebin, window);   % 注意 Sig' 变成列向量

figure;
imagesc(t, f, abs(Spec));
set(gca, 'YDir', 'normal');        % 频率轴向上
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
set(gcf,'Color','w');

% 只看非负频率
axis([t(1) t(end) 0 SampFreq/2]);
colorbar;
title('STFT magnitude');

%% 4. Ridge 提取 + 平滑（extridge_mult）
bw    = SampFreq / 80;     % TF 滤波带宽，可以酌情调
beta1 = 1e-4;              % 曲线平滑参数（extridge_mult 用）
num   = 2;                 % 预期分量个数：自己看 TFR 决定，比如 1/2/3...
delta = 20;                % 相邻时间允许的最大频率索引跳变

[fidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, ...
                                  beta1, bw, Nfrebin, window);
% fidexmult: [num × N]，每一行是一条 ridge 的频率索引
% tfdv:      [num × N]，对应 ridge 上的 TF 幅值

%% 5. 画提取到的 ridge 曲线（在频率轴上）
figure;
set(gcf,'Position',[20 100 640 400]);
set(gcf,'Color','w');

% f(fidexmult) 把频率索引映射成实际 Hz，形状仍为 [num × N]
plot(t, f(fidexmult), 'LineWidth', 2);
xlabel('Time / Sec','FontSize',14,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',14,'FontName','Times New Roman');
set(gca,'FontSize',14);
set(gca,'LineWidth',1.5);
axis([t(1) t(end) 0 SampFreq/2]);
grid on;
title('Detected ridge curves (initial IFs from STFT)');

%% 6. 构造 NCME 的初始 IF（iniIF），并按"重要程度"排序
% extridge_mult 已经给出每条 ridge 的频率索引 fidexmult，
% 对应频率 Hz 就是 iniIF_all = f(fidexmult)，尺寸 [num × N]。
iniIF_all = f(fidexmult);           % 每一行是一条 ridge 的 IF(t)，单位 Hz

% 用 tfdv 的能量评估"哪条 ridge 更主导"
% 这里用沿时间积分的能量 sum(|TFR|^2) 作为重要程度
ridge_energy = sum(tfdv.^2, 2);     % [num × 1]

% 按能量从大到小排序，保证"更主要的分量在前面"
[~, order] = sort(ridge_energy, 'descend');

iniIF_sorted = iniIF_all(order, :); % 排序后的 IF 初值，行顺序已从大到小

% 如果你想只让 NCME 先估计最强的那一条，可以只取第一行：
% iniIF = iniIF_sorted(1,:);   % 单模态用法（K=1）
% 下面我按多模态接口写（K = num），但注意你那版 NCME 代码目前只更新第 1 行

iniIF = iniIF_sorted;             % [K × N]，K=num，且主分量在前

%% 7. 调用 NCME 做 IF / IA / 模式估计
% 参数可以参考作者 demo 代码，之后再根据效果微调
lambda = 5e-4;        % LASSO 正则权重（稀疏性）
beta   = 1e-6;        % IF 增量平滑参数（越小越平滑）
tol    = 1e-8;        % 收敛阈值

[IFmset, IA, smset] = NCME(Sig, SampFreq, iniIF, lambda, beta, tol);

% 说明：
%   IFmset: [K × N × iter]，每次迭代的 IF 轨迹集合
%   smset:  [K × N × iter]，每次迭代的模式重构集合
%   IA:     [1 × N]，当前版本 NCME 里是"第 1 个模式"的瞬时幅度估计

K        = size(iniIF,1);
[~,~,it] = size(IFmset);

% 目前这版 NCME 代码实际上只在迭代里更新了 i=1 这一行，
% 所以我们主要看第 1 个模式（能量最大的那条 ridge）。
IF_est_final = squeeze(IFmset(1,:,end));   % [1 × N] 估计的 IF
s_est_final  = squeeze(smset(1,:,end));    % [1 × N] 估计的模式（重构信号）
IA_est       = IA;                         % [1 × N] 估计的 IA

%% 8. 画 NCME 结果：重构信号 / IF / IA
figure;
fsize = 12;

% (a) 重构信号 vs 原始信号
subplot(3,1,1);
plot(t, s_est_final, 'r', 'LineWidth', 1);   % 估计模式（主分量）
hold on;
plot(t, Sig, 'b', 'LineWidth', 1);          % 原始信号
xlabel({'Time (s)','(a)'});
ylabel('Amplitude');
legend('Estimated mode (NCME)','Original signal');
set(gca,'FontSize',fsize,'LineWidth',1);
grid on;

% (b) IF：估计 vs 初始 ridge IF（无真值，只能对比初值）
subplot(3,1,2);
plot(t, IF_est_final, 'r', 'LineWidth', 1);      % NCME 最终 IF
hold on;
plot(t, iniIF(1,:), 'w--', 'LineWidth', 1.2);      % 最强 ridge 的初始 IF
xlabel({'Time (s)','(b)'});
ylabel('Frequency (Hz)');
legend('Estimated IF (NCME)','Initial IF from ridge');
ylim([0 SampFreq/2]);
set(gca,'FontSize',fsize,'LineWidth',1);
grid on;

% (c) IA：估计的瞬时幅度
subplot(3,1,3);
plot(t, IA_est, 'r', 'LineWidth', 1);
xlabel({'Time (s)','(c)'});
ylabel('Amplitude');
title('Estimated instantaneous amplitude (dominant mode)');
set(gca,'FontSize',fsize,'LineWidth',1);
grid on;

%% 9. （可选）重构 vs 原始，再单独画一张
figure;
plot(t, Sig, 'b', 'LineWidth', 1); 
hold on;
plot(t, s_est_final, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original signal','Estimated dominant mode');
set(gca,'FontSize',fsize,'LineWidth',1);
grid on;
title('Reconstructed dominant mode vs original signal');

% 如果你有 SNR.m，可以试一下：
% SNR_signal = SNR(Sig, s_est_final)
% SNR_IF     = SNR(iniIF(1,:), IF_est_final)
