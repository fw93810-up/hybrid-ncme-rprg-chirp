function [mode, IA, IF_used, xm, ym] = NCME_amp_only(s, fs, eIF, lambda, w)
% 固定 IF 的简化版 NCME（只估计振幅），带样本权重 w
%
% 输入:
%   s     : 信号（行向量或列向量）
%   fs    : 采样频率
%   eIF   : IF(t) 轨迹（Hz），1×N
%   lambda: LASSO 正则参数
%   w     : (可选) 长度为 N 的权重向量；越大表示该时刻误差越重要
%
% 输出:
%   mode    : 重构的模态信号
%   IA      : 瞬时幅度 sqrt(x^2 + y^2)
%   IF_used : 实际使用的 IF（就是 eIF）
%   xm, ym  : cos / sin 通道的系数

    % --- 预处理 ---
    s   = s(:).';    % 行向量
    eIF = eIF(:).';  % 行向量
    N   = length(s);

    if length(eIF) ~= N
        error('NCME_amp_only: eIF length must match signal length.');
    end

    % 如果没给权重，就默认全 1（和原来完全一致）
    if nargin < 5 || isempty(w)
        w = ones(N, 1);
    else
        w = w(:);
        if length(w) ~= N
            error('NCME_amp_only: weight vector w must have length N.');
        end
    end

    t = (0:N-1)/fs;

    % --- 由 IF 构造相位并生成 cos/sin 基底 ---
    phi  = cumtrapz(t, eIF);          % φ(t) = ∫ f(τ) dτ
    cosm = cos(2*pi*phi);
    sinm = sin(2*pi*phi);

    % 对角矩阵形式
    Am = spdiags(cosm.', 0, N, N);
    Bm = spdiags(sinm.', 0, N, N);
    A  = [Am Bm];

    % --- 二阶差分算子，用作 D ---
    e  = ones(N,1);
    e2 = -2*e;
    e2(1)   = -1;
    e2(end) = -1;
    oper    = spdiags([e e2 e], -1:1, N, N);
    D       = blkdiag(oper, oper);

    % --- 构造加权后的 Ã, b̃ ---
    b   = s.';                      % 列向量
    W   = spdiags(w, 0, N, N);      % 对角权重矩阵
    A_w = W * A;
    b_w = W * b;

    % --- ADMM LASSO 求解 ---
    rho = 1;
    [CVX_u_v, ~] = gene_lasso(A_w, D, b_w, lambda, rho, 1);

    xm = CVX_u_v(1:N).';
    ym = CVX_u_v(N+1:end).';

    % --- 模态重构与瞬时幅度 ---
    mode     = xm .* cosm + ym .* sinm;
    IA       = sqrt(xm.^2 + ym.^2);
    IF_used  = eIF;
end
