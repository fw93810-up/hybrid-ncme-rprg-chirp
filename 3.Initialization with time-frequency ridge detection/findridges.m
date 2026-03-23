function index = findridges(Spec, delta)
% Ridge detection algorithm
% 输入:
%   Spec  : 复数 TFR 矩阵 (freq × time)
%   delta : 相邻时间点之间允许的最大频率 bin 跳变
% 输出:
%   index : 1×N，每个时间点的频率索引

    % 用幅度做 ridge 检测
    Spec = abs(Spec);

    % 把 NaN / Inf 清理掉，避免 max 得到 NaN
    Spec(~isfinite(Spec)) = 0;

    [M, N] = size(Spec);
    index = zeros(1, N);   % 先初始化为 0（"没找到"）

    % 空矩阵或全 0：直接返回，交给外面处理
    if M == 0 || N == 0 || ~any(Spec(:) > 0)
        return;
    end

    % 不用 ==max(...)，直接拿 max 的索引
    [~, lin_idx] = max(Spec(:));            % 线性索引
    [fmax, tmax] = ind2sub([M, N], lin_idx);

    index(tmax) = fmax;

    % 向右边（未来）延伸 ridge
    f0 = fmax;
    for j = (min(tmax+1, N)):N
        low = max(1, f0-delta);
        up  = min(M, f0+delta);
        [~, rel_idx] = max(Spec(low:up, j));
        f0 = rel_idx + low - 1;
        index(j) = f0;
    end

    % 向左边（过去）延伸 ridge
    f1 = fmax;
    for j = (max(1, tmax-1)):-1:1
        low = max(1, f1-delta);
        up  = min(M, f1+delta);
        [~, rel_idx] = max(Spec(low:up, j));
        f1 = rel_idx + low - 1;
        index(j) = f1;
    end
end
