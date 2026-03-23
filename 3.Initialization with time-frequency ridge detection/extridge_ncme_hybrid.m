function [fidexmult, tfdv, modes, IFs, IAs, Sig_res] = ...
    extridge_ncme_hybrid(Sig, SampFreq, num_modes, delta, ...
                          lambda_ncme, bw, beta_sIF, ...
                          Nfrebin, window)
% 结合 ridge 提取 + Dechirp_filter 平滑 IF + 加权 NCME_amp_only 的多模态分解
%
% 输入:
%   Sig        : 原始信号（行/列向量）
%   SampFreq   : 采样频率
%   num_modes  : 希望提取的模态数
%   delta      : findridges 中允许的最大频率跳变（bin）
%   lambda_ncme: NCME_amp_only 的 LASSO 正则
%   bw         : Dechirp_filter 的 TF 滤波带宽 (Hz)，如 SampFreq/80
%   beta_sIF   : Dechirp_filter 中 IF 平滑参数
%   Nfrebin    : STFT 频率点数
%   window     : STFT 窗长
%
% 输出:
%   fidexmult  : [M × N] 每个模态对应的 ridge 频率索引（M 为实际模态数）
%   tfdv       : [M × N] ridge 上的 TF 幅度
%   modes      : [M × N] 各模态的重构信号
%   IFs        : [M × N] 各模态使用的 IF(t)（Hz）
%   IAs        : [M × N] 各模态 IA
%   Sig_res    : 最后一轮后的残差信号

    % ---------- 全局参数：权重相关 ----------
    band_Hz   = 5;     % 在 sIF 附近 +/- band_Hz 的频带内看能量
    sigma_w   = 0.0001;   % 高斯权重的标准差，控制"重/轻"的过渡
    w_min     = 0.1;   % 权重下界，避免某些点完全被忽略

    % ---------- 预处理 ----------
    if isreal(Sig)
        Sig = hilbert(Sig);
    end
    Sig = Sig(:).';
    N   = length(Sig);
    t   = (0:N-1)/SampFreq;

    fidexmult = zeros(num_modes, N);
    tfdv      = zeros(num_modes, N);
    modes     = zeros(num_modes, N);
    IFs       = zeros(num_modes, N);
    IAs       = zeros(num_modes, N);

    Sig_res = Sig;
    clr     = lines(num_modes);

    for k = 1:num_modes

        % ----- 1) 当前残差信号的 STFT -----
        [Spec, f] = STFT(Sig_res(:), SampFreq, Nfrebin, window);
        [nF, Ntime] = size(Spec);
        if Ntime ~= N
            error('STFT time length mismatch.');
        end

        % 画当前残差 TFR
        figure;
        imagesc(t, f, abs(Spec));
        set(gca,'YDir','normal');
        xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
        ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
        set(gca,'FontSize',12,'LineWidth',1);
        set(gcf,'Color','w');
        axis([t(1) t(end) 0 SampFreq/2]);
        colorbar;
        title(sprintf('Residual TFR before extracting mode %d', k));
        hold on;

        % ----- 2) 用 findridges 在当前 TFR 上找一条 ridge -----
        c = findridges(Spec, delta);   % 1×N

        if all(c == 0)
            fprintf('Mode %d: no valid ridge found, stop decomposition.\n', k);
            num_eff   = k - 1;
            fidexmult = fidexmult(1:num_eff, :);
            tfdv      = tfdv(1:num_eff, :);
            modes     = modes(1:num_eff, :);
            IFs       = IFs(1:num_eff, :);
            IAs       = IAs(1:num_eff, :);
            return;
        end

        ridge_freq = f(c);             % 粗略 IF，Hz

        % ----- 3) 用 Dechirp_filter 得到平滑 IF sIF -----
        [sIF, ~] = Dechirp_filter(Sig_res, SampFreq, bw, ridge_freq, beta_sIF);
        sIF = sIF(:).';
        IFs(k,:) = sIF;

        % 在 TFR 上画出平滑后的 IF
        plot(t, sIF, '--', 'Color', clr(k,:), 'LineWidth', 1.5);

        % ----- 4) 计算这一条脊线的权重 w_n（高斯函数控制） -----
        % 4.1 把 "Hz 带宽" 转成 "频率 bin 半径"
        if length(f) > 1
            df = mean(diff(f));
        else
            df = 1;
        end
        band_bins = max(1, round(band_Hz / df));

        ridge_amp = zeros(1, N);
        Spec_abs  = abs(Spec);

        for j = 1:N
            % 找到最接近 sIF(j) 的频率 bin
            [~, idx_center] = min(abs(f - sIF(j)));
            low = max(1, idx_center - band_bins);
            up  = min(nF, idx_center + band_bins);
            % 该时间点沿脊线带状区域的能量（这里用最大值，也可以用求和）
            ridge_amp(j) = max(Spec_abs(low:up, j));
        end

        % 归一化到 [0,1]
        maxA = max(ridge_amp);
        if maxA <= 0
            w = ones(N,1);   % 极端情况下就不加权
        else
            a_norm = ridge_amp / maxA;
            % 高斯权重：离主脊线能量越强，权重越接近 1
            w = exp(-((1 - a_norm).^2) / (2 * sigma_w^2));
            w = max(w, w_min);     % 保底，避免完全 0
            w = w(:);
        end

        % ----- 5) 由平滑 IF + 权重 w 调用 NCME_amp_only -----
        [mode_k, IA_k, ~] = NCME_amp_only(Sig_res, SampFreq, sIF, lambda_ncme, w);

        modes(k,:) = mode_k;
        IAs(k,:)   = IA_k;

        % 也顺便更新 fidexmult / tfdv，用 sIF 对应的索引
        findex = zeros(1,N);
        for j = 1:N
            [~, findex(j)] = min(abs(f - sIF(j)));
            tfdv(k,j) = Spec_abs(findex(j), j);
        end
        fidexmult(k,:) = findex;

        % ----- 6) 更新残差信号 -----
        Sig_res = Sig_res - mode_k;
    end

    % 如果能完整跑完 num_modes 次，就全部有效
end
