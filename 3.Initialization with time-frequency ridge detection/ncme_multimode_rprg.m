function out = ncme_multimode_rprg(Sig, SampFreq, num_modes, lambda_ncme, thrf)
% ncme_multimode_rprg
% 多模态 NCME + RPRG 联合分解管线（函数版本）
%
% out = ncme_multimode_rprg(Sig, SampFreq, num_modes, lambda_ncme, thrf)
%
% 输入：
%   Sig         : 1×N 或 N×1 实/复数信号
%   SampFreq    : 采样频率（Hz），只影响坐标刻度
%   num_modes   : 预期分解成的模态个数（例如 4）
%   lambda_ncme : NCME_multi 的 λ（LASSO 稀疏正则），例如 1e4
%   thrf        : RPRG 的频率索引差阈值（bin），控制"多近算交叉"，
%                 非交叉信号建议从 0 或 1 开始试
%
% 输出（结构体 out）：
%   out.Sig_orig     : 原始信号（1×N）
%   out.t            : 时间轴（1×N）
%   out.f            : 频率轴（length(f)×1）
%   out.ridgeIF      : 每轮 extridge_mult 找到的 ridge 初始 IF (num_modes×N)
%   out.IFs_init     : RPRG/最终微调前，NCME_multi 得到的 IF (num_modes×N)
%   out.modes_init   : 同上，对应的时域模态 (num_modes×N)
%   out.IAs_init     : 同上，对应 IA (num_modes×N)
%   out.eIF_reg      : RPRG 重组后的 IF 初值 (K_reg×N)
%   out.findex_reg   : RPRG 重组后的索引矩阵 (K_reg×N)
%   out.IFs_final    : 最终 NCME_multi 微调后的 IF (K_reg×N)
%   out.modes_final  : 最终时域模态 (K_reg×N)
%   out.IAs_final    : 最终 IA (K_reg×N)
%   out.err_sig      : 原始信号 - 模态和 (1×N)
%   out.interset     : RPRG 返回的交叉区间信息
%
% 依赖函数：
%   STFT.m, extridge_mult.m, Dechirp_filter.m, findridges.m,
%   NCME_multi.m, RPRG.m 及 NCME_multi 的内部依赖 gene_lasso.m, Differ.m 等
%
% 调用示例：
%   load data.mat             % 假设里面变量叫 data
%   Sig = data(:)';
%   out = ncme_multimode_rprg(Sig, 250, 4, 1e4, 1);
%

    %% 0. 参数默认值处理
    if nargin < 2 || isempty(SampFreq)
        SampFreq = 1;   % 没给就当 1 Hz
    end
    if nargin < 3 || isempty(num_modes)
        num_modes = 4;  % 默认分 4 个模态
    end
    if nargin < 4 || isempty(lambda_ncme)
        lambda_ncme = 1e4;
    end
    if nargin < 5 || isempty(thrf)
        thrf = 1;       % RPRG 阈值：非交叉信号先取 0 或 1
    end

    % 强制成行向量
    Sig = Sig(:).';  

    %% 1. 基本变量
    Sig_orig = Sig;               % 原始信号
    Sig_res  = Sig;               % 残差信号，用于逐轮 ridge 检测

    N = length(Sig);
    t = (0:N-1) / SampFreq;

    %% 2. 原始信号 STFT（做 TFR 背景）
    window  = 256;
    Nfrebin = 1024;

    [Spec0, f] = STFT(Sig_orig.', SampFreq, Nfrebin, window);

    %% 3. 主循环：每次 1 条 ridge + 一次 NCME_multi 对所有已发现模态联合优化
    delta      = 20;            % findridges 中允许频率跳变（bin 数）
    bw         = SampFreq / 80; % Dechirp_filter 滤波带宽
    beta_ridge = 1e-4;          % ridge 平滑参数

    beta_ncme  = 1e-6; 
    tol_ncme   = 1e-8;

    % 结果存储
    modes    = zeros(num_modes, N);   % 每模态重构信号
    IFs      = zeros(num_modes, N);   % 每模态 IF
    IAs      = zeros(num_modes, N);   % 每模态 IA
    ridgeIF  = zeros(num_modes, N);   % 每轮 extridge_mult 提供的 ridge 初值
    eIF_all  = [];                    % 当前所有分量的 IF，用于下一轮初值

    for k = 1:num_modes
        fprintf('==== 提取第 %d 条模态 ====\n', k);

        % --- (1) 在当前残差上用 extridge_mult 找最强 ridge（num=1） --- %
        [fidex_k, ~] = extridge_mult(Sig_res, SampFreq, 1, delta, ...
                                     beta_ridge, bw, Nfrebin, window);
        ridgeIF_k = f(fidex_k(1,:));      % 这一轮新分量的平滑 IF 初值（Hz）
        ridgeIF(k,:) = ridgeIF_k;

        % --- (2) 构造本轮 NCME_multi 的 IF 初值矩阵 eIF_init (k×N) --- %
        if k == 1
            eIF_init = ridgeIF_k;         % 1×N
        else
            eIF_init          = zeros(k, N);
            eIF_init(1:k-1,:) = eIF_all(1:k-1,:);
            eIF_init(k,:)     = ridgeIF_k;
        end

        % --- (3) 在原始信号上跑一次 k 模态 NCME_multi，全局重构 --- %
        [IFmset_k, IA_k, smset_k] = NCME_multi(Sig_orig, SampFreq, ...
                                               eIF_init, lambda_ncme, ...
                                               beta_ncme, tol_ncme);
        % IFmset_k: k×N×T
        % smset_k : k×N×T
        % IA_k    : k×N

        IF_ref    = IFmset_k(:,:,end);   % k×N
        modes_ref = smset_k(:,:,end);    % k×N
        IA_ref    = IA_k;                % k×N

        % 保存当前所有分量的结果（覆盖前一轮）
        eIF_all(1:k,:) = IF_ref(1:k,:);
        IFs(1:k,:)     = IF_ref(1:k,:);
        modes(1:k,:)   = modes_ref(1:k,:);
        IAs(1:k,:)     = IA_ref(1:k,:);

        % --- (4) 用当前所有 k 个分量的重构结果更新残差 --- %
        Sig_res = Sig_orig - sum(modes(1:k,:), 1);
    end

    %% 4. 把 IFs 映射成频率索引矩阵，交给 RPRG 重组路径

    [K_init, ~] = size(IFs);      % K_init = num_modes
    M_f = length(f);
    f = f(:);                     % 列向量

    findex_all = zeros(K_init, N);
    for k = 1:K_init
        for n = 1:N
            [~, idx] = min(abs(f - IFs(k,n)));
            idx = max(1, min(M_f, idx)); % 防止越界
            findex_all(k,n) = idx;
        end
    end

    % RPRG：路径重组
    [findex_reg, interset] = RPRG(findex_all, thrf);
    [K_reg, ~] = size(findex_reg);

    % 重组后索引 → Hz，得到 eIF_reg
    eIF_reg = zeros(K_reg, N);
    for k = 1:K_reg
        idx_k = findex_reg(k,:);
        idx_k = max(1, min(M_f, idx_k));
        eIF_reg(k,:) = f(idx_k);
    end

    %% 5. 用 eIF_reg 再跑一次 NCME_multi 微调

    lambda_fine = lambda_ncme;    % 你可以单独调一个 lambda_fine
    beta_fine   = beta_ncme;
    tol_fine    = tol_ncme;

    [IFmset_final, IA_final, smset_final] = NCME_multi(Sig_orig, SampFreq, ...
                                                       eIF_reg, lambda_fine, ...
                                                       beta_fine, tol_fine);

    IFs_final   = IFmset_final(:,:,end);   % K_reg×N
    modes_final = smset_final(:,:,end);    % K_reg×N

    % 更新 num_modes 为重组后的条数
    num_modes_final = size(IFs_final,1);

    %% 6. 计算误差
    sum_modes = sum(modes_final, 1);
    err_sig   = Sig_orig - sum_modes;

    %% 7. 一些标准图像（和你脚本类似）

    % 7.1 原始 STFT + ridge IF + 初始 NCME IF
    figure;
    imagesc(t, f, abs(Spec0));
    set(gca,'YDir','normal');
    xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
    ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
    set(gca,'FontSize',12,'LineWidth',1);
    set(gcf,'Color','w');
    axis([t(1) t(end) 0 SampFreq/2]);
    colorbar;
    title('STFT with ridge IFs (white dashed) and NCME-multi IFs (colored)');
    hold on;

    clr0 = lines(K_init);
    plot(t, ridgeIF', 'w--', 'LineWidth', 1.2);   % ridge 初值
    for k = 1:K_init
        plot(t, IFs(k,:), 'Color', clr0(k,:), 'LineWidth', 1.6);
    end

    % 7.2 STFT + RPRG 后 IF + 最终微调 IF
    figure;
    imagesc(t, f, abs(Spec0));
    set(gca,'YDir','normal');
    xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
    ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
    set(gca,'FontSize',12,'LineWidth',1);
    set(gcf,'Color','w');
    axis([t(1) t(end) 0 SampFreq/2]);
    colorbar;
    title('TFR + RPRG IF (yellow) + final NCME IF (colored)');
    hold on;

    clr = lines(num_modes_final);
    plot(t, eIF_reg', 'y-.', 'LineWidth', 1.0);  % RPRG 后 IF
    for k = 1:num_modes_final
        plot(t, IFs_final(k,:), 'Color', clr(k,:), 'LineWidth', 1.6);
    end

    % 7.3 Sum of modes vs Error
    figure;
    set(gcf,'Position',[50 50 700 500]);
    set(gcf,'Color','w');

    subplot(2,1,1);
    plot(t, sum_modes, 'r', 'LineWidth', 1); 
    hold on;
    plot(t, Sig_orig, 'c--', 'LineWidth', 0.8); 
    xlabel('Time (s)');
    ylabel('Amp');
    legend('Sum of all modes','Original signal','Location','best');
    title('(a) Sum of all estimated modes vs original');
    set(gca,'FontSize',12,'LineWidth',1);
    grid on;

    subplot(2,1,2);
    plot(t, err_sig, 'b', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Amp');
    title('(b) Error = Original - Sum of modes');
    set(gca,'FontSize',12,'LineWidth',1);
    grid on;

    % 7.4 每个最终模态单独画一格（时域）
    figure;
    set(gcf,'Position',[100 100 800 600]);
    set(gcf,'Color','w');

    for k = 1:num_modes_final
        if num_modes_final == 4
            subplot(2,2,k);
        else
            subplot(num_modes_final,1,k);
        end

        plot(t, modes_final(k,:), 'r', 'LineWidth', 1.0); 
        hold on;
        plot(t, Sig_orig, 'b:', 'LineWidth', 0.5);

        xlabel('Time (s)');
        ylabel('Amp');
        title(sprintf('Final Mode %d (NCME\\_multi + RPRG)', k));
        set(gca,'FontSize',12,'LineWidth',1);
        grid on;
    end

    %% 8. 打包输出结构体
    out.Sig_orig    = Sig_orig;
    out.t           = t;
    out.f           = f;
    out.ridgeIF     = ridgeIF;
    out.IFs_init    = IFs;
    out.modes_init  = modes;
    out.IAs_init    = IAs;

    out.eIF_reg     = eIF_reg;
    out.findex_reg  = findex_reg;
    out.interset    = interset;

    out.IFs_final   = IFs_final;
    out.modes_final = modes_final;
    out.IAs_final   = IA_final;
    out.err_sig     = err_sig;

end
