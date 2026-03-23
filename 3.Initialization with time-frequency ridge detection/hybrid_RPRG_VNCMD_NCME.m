function out = hybrid_RPRG_VNCMD_NCME(Sig, fs, cfg)
% HYBRID_RPRG_VNCME_NCME
% 1) 用 extridge_mult + RPRG (+ 可选 VNCMD) 找 num_strong 条强模态 IF 初值
% 2) 在原始信号上对这 num_strong 条跑一次 NCME_multi，多模态收敛
% 3) 从第 num_strong+1 条开始：
%       - 在残差信号上找新的 ridge IF
%       - 把"已有 IF + 新 IF"一起交给 NCME_multi，在原始信号上做全局多模态优化
%
% 输入：
%   Sig : 1×N 实信号
%   fs  : 采样频率
%   cfg : 结构体，需至少包含：
%       cfg.num_strong, cfg.num_total
%       cfg.STFT.winLen, cfg.STFT.Nfrebin
%       cfg.Ridge.delta, cfg.Ridge.beta_ridge, cfg.Ridge.bw
%       cfg.RPRG.thrf_scale
%       cfg.VNCMD.use, cfg.VNCMD.alpha, cfg.VNCMD.beta, cfg.VNCMD.var, cfg.VNCMD.tol
%       cfg.NCME.lambda, cfg.NCME.beta, cfg.NCME.tol, cfg.NCME.maxIter
%
% 输出 out 结构体（见函数末尾）

Sig = Sig(:)';       % 确保行向量
N   = length(Sig);
t   = (0:N-1)/fs;

num_strong = cfg.num_strong;
num_total  = cfg.num_total;

% ---------- 1. STFT ----------
winLen  = cfg.STFT.winLen;
Nfrebin = cfg.STFT.Nfrebin;

[Spec0, f] = STFT(Sig(:), fs, Nfrebin, winLen);   % Spec0: Nfrebin × N
SpecMag = abs(Spec0);

% ---------- 2. 用 extridge_mult + RPRG 找强模态 IF 初值 ----------
% 先用 extridge_mult 提 num_strong 条 ridge（顺序提取）
bw         = cfg.Ridge.bw;
delta      = cfg.Ridge.delta;
beta_ridge = cfg.Ridge.beta_ridge;

[fidexmult, ~] = extridge_mult(Sig, fs, num_strong, delta, ...
                               beta_ridge, bw, Nfrebin, winLen);
% fidexmult: num_strong × N, 每行是一个分量的频率索引

% RPRG 重新归组，处理交叉
thrf = length(f) / cfg.RPRG.thrf_scale;
[findex_rprg, interset] = RPRG(fidexmult, thrf); %#ok<NASGU>

K_rprg = size(findex_rprg,1);
K_strong = min(num_strong, K_rprg);      % 实际可用强模态条数

IF_init_strong = f(findex_rprg(1:K_strong,:));   % K_strong × N, Hz

% ---------- 3. 可选：用 VNCMD 对强模态 IF 进行微调 ----------
if isfield(cfg,'VNCMD') && isfield(cfg.VNCMD,'use') && cfg.VNCMD.use
    alpha_v = cfg.VNCMD.alpha;
    beta_v  = cfg.VNCMD.beta;
    var_v   = cfg.VNCMD.var;
    tol_v   = cfg.VNCMD.tol;

    % 按原论文做法，先 smooth 一下 IF
    iniIF_v = curvesmooth(IF_init_strong, beta_v);
    fprintf('Running VNCMD on strong components...\n');
    [IFmset_v, IA_v, smset_v] = VNCMD(Sig, fs, iniIF_v, ...
                                      alpha_v, beta_v, var_v, tol_v);
    IF_strong0 = IFmset_v(:,:,end);   % K_strong × N
else
    IF_strong0 = IF_init_strong;
    IA_v = [];
    smset_v = [];
end

% ---------- 4. 第一次：只对强模态做多模态 NCME_multi ----------
lambda_ncme = cfg.NCME.lambda;
beta_ncme   = cfg.NCME.beta;
tol_ncme    = cfg.NCME.tol;

fprintf('NCME_multi: strong phase, K = %d\n', K_strong);
[IFmset_s, IA_s, smset_s] = NCME_multi(Sig, fs, IF_strong0, ...
                                       lambda_ncme, beta_ncme, tol_ncme, num_strong*60);

IF_strong_ncme = IFmset_s(:,:,end);   % K_strong × N
modes_strong   = smset_s(:,:,end);    % K_strong × N

% 初始化总结果数组（后面会逐步填满）
IF_ncme_all   = zeros(num_total, N);
modes_all     = zeros(num_total, N);
IA_ncme_all   = zeros(num_total, N);

IF_ncme_all(1:K_strong,:) = IF_strong_ncme;
modes_all(1:K_strong,:)   = modes_strong;
IA_ncme_all(1:K_strong,:) = IA_s(1:K_strong,:);

% 残差（仅供后面脊线检测时用）
Sig_res = Sig - sum(modes_all(1:K_strong,:), 1);

% 记录强模态的初始 IF（方便对比）
IF_rprg_strong = IF_init_strong;

% ---------- 5. 从第 K_strong+1 条模态开始：逐条加入弱模态 ----------
for k = (K_strong+1) : num_total
    fprintf('--- Adding weak mode %d (global NCME_multi) ---\n', k);

    % (1) 在当前残差信号上提 1 条 ridge（用 Dechirp 的平滑 IF）
    [Spec_res, f_res] = STFT(Sig_res(:), fs, Nfrebin, winLen);
    Spec_res = abs(Spec_res);
    c = findridges(Spec_res, delta);     % 1×N 索引
    % 用 Dechirp_filter 获得平滑 sIF
    [sIF, extr_Sig_dummy] = Dechirp_filter(Sig_res, fs, ...
                                           bw, f_res(c), beta_ridge); %#ok<NASGU>
    IF_new_init = sIF(1,:);             % 1×N

    % (2) 构造这轮 NCME_multi 的初始 IF 矩阵（前 k-1 条用上一次 IF_ncme_all）
    K_cur = k;
    eIF_init_k          = zeros(K_cur, N);
    eIF_init_k(1:K_cur-1,:) = IF_ncme_all(1:K_cur-1,:);
    eIF_init_k(K_cur,:)     = IF_new_init;

    % (3) 在原始信号上，对 K_cur 条模态做一次全局多模态 NCME_multi
    fprintf('NCME_multi: global phase, K = %d\n', K_cur);
    [IFmset_k, IA_k, smset_k] = NCME_multi(Sig, fs, eIF_init_k, ...
                                           lambda_ncme, beta_ncme, tol_ncme, cfg.NCME.maxIter);

    IF_ref    = IFmset_k(:,:,end);   % K_cur × N
    modes_ref = smset_k(:,:,end);    % K_cur × N

    IF_ncme_all(1:K_cur,:) = IF_ref(1:K_cur,:);
    modes_all(1:K_cur,:)   = modes_ref(1:K_cur,:);
    IA_ncme_all(1:K_cur,:) = IA_k(1:K_cur,:);

    % (4) 更新残差，供下一轮脊线检测使用
    Sig_res = Sig - sum(modes_all(1:K_cur,:), 1);
end

% ---------- 6. 打包输出 ----------
out = struct();
out.Sig          = Sig;
out.fs           = fs;
out.t            = t;
out.f            = f;
out.Spectrogram  = SpecMag;           % |STFT|
out.IF_rprg      = IF_rprg_strong;    % RPRG(+VNCMD) 得到的强模态初始 IF (K_strong×N)
out.IF_ncme      = IF_ncme_all;       % 最终所有 NCME_multi IF (num_total×N)
out.modes        = modes_all;         % 最终各模态时域波形 (num_total×N)
out.IA           = IA_ncme_all;       % IA 估计 (num_total×N)
out.residual     = Sig_res;           % 最终残差
out.num_strong   = K_strong;
out.num_total    = num_total;

end
