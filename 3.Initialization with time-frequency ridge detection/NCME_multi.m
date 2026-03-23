function [IFmset, IA, smset] = NCME_multi(s, fs, eIF, lambda, beta, tol, maxIter)
% 多模态 NCME：给定 K 条初始 IF，一次性全局估计 K 个非线性chirp分量
%
% 输入:
%   s      : 测量信号（一维，行/列向量都可以）
%   fs     : 采样频率
%   eIF    : 初始 IF，大小 K×N（每行一条 IF 曲线，单位 Hz）
%   lambda : LASSO 正则系数（典型 10~1000，根据信号能量调）
%   beta   : IF 增量平滑的上限（原文里的 beta）
%   tol    : 收敛阈值
%
% 输出:
%   IFmset : K×N×T，T 是迭代轮数（包括初始）
%   IA     : K×N，最终的瞬时幅值 sqrt(x_k^2 + y_k^2)
%   smset  : K×N×T，每轮迭代的模式重构信号

    %---------------- 预处理 ----------------%
    
    % 确保是行向量
    s = s(:).';
    [K, N] = size(eIF);   % K：分量数, N：样本数

    % --------- 迭代上限设置 ---------
    if nargin < 7 || isempty(maxIter)
        iternum = 60;           % 默认就是 50 步（弱分量那里不传参，用这个）
    else
        iternum = maxIter;      % 外面传多少就用多少
    end
    % -------------------------------


    if length(s) ~= N
        error('NCME_multi: length(s) must equal N = size(eIF,2).');
    end

    t = (0:N-1) / fs;

    % 二阶差分算子 oper 和 opedoub（和原 NCME 相同）
    e  = ones(N,1);
    e2 = -2*e;
    e2(1)   = -1;
    e2(end) = -1;
    oper    = spdiags([e e2 e], -1:1, N, N);   % N×N
    opedoub = oper' * oper;                   % N×N

    % 容器
    sinm = zeros(K, N);
    cosm = zeros(K, N);

    IFsetiter = zeros(K, N, iternum+1);
    IFsetiter(:,:,1) = eIF;

    ssetiter = zeros(K, N, iternum+1);   % 每轮的模式重构
    sDif     = tol + 1;                  % 收敛指标
    iter     = 1;

    % 初始 sin/cos
    for i = 1:K
        phi_i     = cumtrapz(t, eIF(i,:));          % 相位 φ_i(t) = ∫ f_i(τ)dτ
        sinm(i,:) = sin(2*pi*phi_i);
        cosm(i,:) = cos(2*pi*phi_i);
    end

    b   = s.';   % 列向量
    rho = 1;     % gene_lasso 中的 ADMM 参数

    %---------------- 迭代 ----------------%
    while (sDif > tol*1e-2 && iter <= iternum)

        % 逐步增大 beta（原代码的调节策略）
        betathr = 10^(iter/36 - 10);
        if betathr > beta
            betathr = beta;
        end

        %----- 1) 构造大字典矩阵 A = [A1 ... AK B1 ... BK] -----%
        % 每个 Ai/Bi 是 N×N，对角矩阵
        A_blocks = cell(1, 2*K);
        for i = 1:K
            Am_i         = spdiags(cosm(i,:).', 0, N, N);
            A_blocks{i}  = Am_i;
        end
        for i = 1:K
            Bm_i         = spdiags(sinm(i,:).', 0, N, N);
            A_blocks{K+i}= Bm_i;
        end
        A = [A_blocks{:}];     % N × (2*K*N)

        %----- 2) 构造正则矩阵 D = blkdiag(oper,...,oper) (2K块) -----%
        D_blocks = cell(1, 2*K);
        for i = 1:2*K
            D_blocks{i} = oper;
        end
        D = blkdiag(D_blocks{:});      % (2*K*N) × (2*K*N)

        %----- 3) 用 gene_lasso 解 LASSO -----%
        %   min 1/2 ||A u - b||^2 + lambda ||D u||_1
        [CVX_u_v, ~] = gene_lasso(A, D, b, lambda, rho, 1);

        % 把 u 拆回每个分量的 x_k, y_k
        xm = zeros(K, N);
        ym = zeros(K, N);

        for i = 1:K
            xm(i,:) = CVX_u_v((i-1)*N + 1 : i*N).';
        end
        for i = 1:K
            ym(i,:) = CVX_u_v(K*N + (i-1)*N + 1 : K*N + i*N).';
        end

        %----- 4) 更新每个分量的 IF / 模式 -----%
        sDif = 0;   % 本轮收敛指标

        for i = 1:K
            % (1) 计算 x,y 的时间导数（用 Differ 函数）
            ybar = Differ(ym(i,:), 1/fs);   % dy/dt
            xbar = Differ(xm(i,:), 1/fs);   % dx/dt

            % (2) 频率增量 deltaIF = (x*y' - y*x') / (x^2 + y^2) / (2π)
            denom   = xm(i,:).^2 + ym(i,:).^2;
            % 为避免除零，做个极小正则
            denom(denom < 1e-12) = 1e-12;

            deltaIF = (xm(i,:).*ybar - ym(i,:).*xbar) ./ denom / (2*pi);
            deltaIF = deltaIF(:);   % 列向量

            % (3) 用(2/betathr*opedoub + I)^{-1} 低通滤波平滑 deltaIF
            L      = (2/betathr)*opedoub + speye(N);
            deltaIF_smooth = L \ deltaIF;   % N×1

            % (4) 更新 IF：eIF_new = eIF_old - 0.5 * deltaIF_smooth
            eIF(i,:) = eIF(i,:) - 0.5 * deltaIF_smooth.';   % 更新第 i 条 IF

            % (5) 用更新后的 IF 重算 sin/cos
            phi_i     = cumtrapz(t, eIF(i,:));
            sinm(i,:) = sin(2*pi*phi_i);
            cosm(i,:) = cos(2*pi*phi_i);

            % (6) 更新该分量的重构信号（本轮结果）
            ssetiter(i,:,iter+1) = xm(i,:).*cosm(i,:) + ym(i,:).*sinm(i,:);

            % (7) 累加收敛指标
            if iter == 1
                % 第一轮，上一轮是全 0，防止除以 0
                prev_norm = max(norm(ssetiter(i,:,1)), 1e-12);
            else
                prev_norm = max(norm(ssetiter(i,:,iter)), 1e-12);
            end
            sDif = sDif + (norm(ssetiter(i,:,iter+1) - ssetiter(i,:,iter)) / prev_norm)^2;
        end

        % 保存本轮更新后的 IF
        IFsetiter(:,:,iter+1) = eIF;

        fprintf('NCME_multi: iter = %d, sDif = %.3e\n', iter, sDif);

        iter = iter + 1;
    end

    %---------------- 打包输出 ----------------%
    T = iter;   % 实际迭代到 iter-1，因此切到 1:T
    IFmset = IFsetiter(:,:,1:T);
    smset  = ssetiter(:,:,1:T);

    % 最终瞬时幅值 IA = sqrt(x^2 + y^2)，用最后一轮的 xm, ym
    IA = sqrt(xm.^2 + ym.^2);
end
