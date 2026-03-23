function [fidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, beta, bw, Nfrebin, window)
% Extract ridges for multi-component signals.
% In each iteration, the signal component associated with the extracted ridge is
% reconstructed by a time-frequency (TF) filter and then removed from the signal,
% so that ridge curves of other components with smaller energies can be extracted.
%
% 输入:
%   Sig       : measured signal, row vector
%   SampFreq  : sampling frequency
%   num       : number of signal components to extract
%   delta     : max frequency-bin jump between consecutive ridge points
%   beta      : smoothing parameter in Dechirp_filter
%   bw        : TF filter bandwidth (Hz)
%   Nfrebin   : number of frequency bins (FFT length)
%   window    : STFT window length
%
% 输出:
%   fidexmult : [num × N] ridge frequency *index* for each component
%   tfdv      : [num × N] TF magnitude along each ridge

    % 保证是解析信号
    if isreal(Sig)
        Sig = hilbert(Sig);
    end

    L = length(Sig);
    t = (0:L-1) / SampFreq;   % 时间轴

    fidexmult = zeros(num, L);
    tfdv      = zeros(num, L);

    for i = 1:num
        %----------- 1) 当前残差信号的 STFT (Sig 已是残差) ----------
        [Spec, f] = STFT(Sig(:), SampFreq, Nfrebin, window);

        % 在这里直接画"当前剩余信号"的时频图
        figure;
        imagesc(t, f, abs(Spec));
        set(gca, 'YDir', 'normal');
        xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
        ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
        set(gca,'FontSize',12,'LineWidth',1);
        set(gcf,'Color','w');
        axis([t(1) t(end) 0 SampFreq/2]);   % 只看非负频率
        colorbar;
        title(sprintf('Residual STFT before extracting component %d', i));

        %----------- 2) 在当前残差信号上找 ridge 并提取第 i 个分量 ----------
        c = findridges(Spec, delta);   % ridge detection on current residual

        % TF 滤波，重构第 i 个分量
        [sIF, extr_Sig] = Dechirp_filter(Sig, SampFreq, bw, f(c), beta);
        % sIF: 平滑后的 ridge IF
        % extr_Sig: 时域中的第 i 个提取分量

        % 把 sIF 对应回频率索引，并记录 TF 幅值
        findex = zeros(1, L);
        for j = 1:L
            [~, findex(j)] = min(abs(f - sIF(1,j)));
            tfdv(i,j)      = abs(Spec(findex(j), j));
        end
        fidexmult(i,:) = findex;

        % 如果你也想把这条 sIF 的 ridge 叠加到刚才画的 TFR 上，可以加：
        figure(gcf);   % 回到上一张图
        hold on;
        plot(t, sIF(1,:), 'w--', 'LineWidth', 1.2);  % 白虚线标出第 i 条 ridge

        %----------- 3) 从信号中减去第 i 个分量，更新残差信号 ----------
        Sig = Sig - extr_Sig;   % 这一句就是"去掉第 i 个分量"

        % 下一轮 for 循环时，Sig 就已经是"去掉前 i 个分量后的残差"了，
        % 顶部的 STFT 和图就是当前真正的"剩余部分"的时频图。
    end
end
