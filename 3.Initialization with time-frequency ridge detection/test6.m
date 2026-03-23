% demo_strongweak.m
% Reproduce the "Strong/weak separation + reconstruction" example from data.mat
% Output 3 standalone figures for README:
%   figures/fig_strongweak_stft.png
%   figures/fig_modes_if.png
%   figures/fig_recon.png

clear; close all; clc;

%% 1) Load data
S = load('data.mat');

% --- find signal ---
if isfield(S,'Sig'), Sig = S.Sig;
elseif isfield(S,'sig'), Sig = S.sig;
elseif isfield(S,'x'),   Sig = S.x;
elseif isfield(S,'y'),   Sig = S.y;
else
    fn = fieldnames(S);
    Sig = [];
    for i = 1:numel(fn)
        v = S.(fn{i});
        if isnumeric(v) && isvector(v) && numel(v) > 200
            Sig = v; fprintf('Using "%s" as Sig\n', fn{i}); break;
        end
    end
    assert(~isempty(Sig), 'Cannot find signal vector in data.mat');
end
Sig = Sig(:)';   % row vector

% --- sampling rate ---
fs = 250;  % set manually

% ===== cfg preset for data.mat (strong/weak example) =====
cfg = struct();


cfg.num_strong = 2;   
cfg.num_total  = 4;   

% STFT
cfg.STFT.winLen  = 256;
cfg.STFT.Nfrebin = 1024;

% Ridge / Dechirp 
cfg.Ridge.delta      = 20;          
cfg.Ridge.beta_ridge = 1e-4;        
cfg.Ridge.bw         = fs/80;       

% RPRG
% thrf = length(f)/thrf_scale，thrf_scale 越大 -> thrf越小
cfg.RPRG.thrf_scale = 256;          

% VNCMD
cfg.VNCMD = struct();
cfg.VNCMD.use   = 1;                
cfg.VNCMD.alpha = 1e-5;
cfg.VNCMD.beta  = 1e-5;
cfg.VNCMD.var   = 0.1;
cfg.VNCMD.tol   = 1e-8;

% NCME
cfg.NCME = struct();
cfg.NCME.lambda = 1000;
cfg.NCME.beta   = 1e-6;
cfg.NCME.tol    = 1e-8;
cfg.NCME.maxIter = 30;

%% 3) Run hybrid pipeline
out = hybrid_RPRG_VNCMD_NCME(Sig, fs, cfg);

% sanity checks (avoid "only two plots" because data empty)
assert(isfield(out,'modes') && ~isempty(out.modes), 'out.modes is empty.');
assert(isfield(out,'IF_ncme') && ~isempty(out.IF_ncme), 'out.IF_ncme is empty.');
K = out.num_total;
K = min(K, size(out.modes,1));
K = min(K, size(out.IF_ncme,1));

%% Make output folder
if ~exist('figures','dir'), mkdir('figures'); end

%% 4A) Standalone STFT + IFs (not flat)
fig1 = figure('Color','w','Position',[80 80 1100 700]);  % good aspect ratio
imagesc(out.t, out.f, out.Spectrogram); axis xy;
colormap(turbo); colorbar;
hold on;

if isfield(out,'IF_rprg') && ~isempty(out.IF_rprg)
    plot(out.t, out.IF_rprg, 'w--', 'LineWidth', 1.6);
end
plot(out.t, out.IF_ncme(1:K,:), 'LineWidth', 1.6);

title('STFT + ridge IF init (white dashed) and refined IFs (color)');
xlabel('Time (s)'); ylabel('Frequency (Hz)');

% optional: zoom frequency band if your figure looks too "empty"
% ylim([0 150]);

hold off;

exportgraphics(fig1, fullfile('figures','fig_strongweak_stft.png'), 'Resolution', 250);
fprintf('Saved: figures/fig_strongweak_stft.png\n');

%% 4B) Per-mode IF (on STFT) + estimated modes (time-domain)
Kshow = min(K, 4);

fig2 = figure('Color','w','Position',[80 80 1200 850]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for k = 1:Kshow
    nexttile;
    imagesc(out.t, out.f, out.Spectrogram); axis xy;
    colormap(turbo); colorbar;
    hold on;

    % show ridge init for that mode if available
    if isfield(out,'IF_rprg') && ~isempty(out.IF_rprg) && size(out.IF_rprg,1) >= k
        plot(out.t, out.IF_rprg(k,:), 'w--', 'LineWidth', 1.2);
    end
    plot(out.t, out.IF_ncme(k,:), 'LineWidth', 1.4);

    title(sprintf('Mode %d: ridge IF (white) vs refined IF (color)', k));
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    hold off;
end

exportgraphics(fig2, fullfile('figures','fig_modes_if.png'), 'Resolution', 250);
fprintf('Saved: figures/fig_modes_if.png\n');

% Optional separate time-domain modes plot (if you want)
% fig2b = figure('Color','w','Position',[80 80 1200 700]);
% tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
% for k=1:Kshow
%     nexttile; plot(out.t, out.modes(k,:)); grid on;
%     title(sprintf('Mode %d (estimated)',k)); xlabel('Time (s)'); ylabel('Amp');
% end
% exportgraphics(fig2b, fullfile('figures','fig_modes_time.png'), 'Resolution', 250);

%% 4C) Reconstruction + residual (standalone)
fig3 = figure('Color','w','Position',[80 80 1100 700]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

recon = sum(out.modes(1:K,:), 1);
err = Sig - recon;

nexttile;
plot(out.t, recon, '-', 'LineWidth', 1.2); hold on;
plot(out.t, Sig, '--', 'LineWidth', 1.0);
grid on;
legend('Sum of all modes','Original signal','Location','best');
title('(a) Sum of all estimated modes vs original');
xlabel('Time (s)'); ylabel('Amp');

nexttile;
plot(out.t, err, '-', 'LineWidth', 1.0);
grid on;
title('(b) Error = Original - Sum of modes');
xlabel('Time (s)'); ylabel('Amp');

exportgraphics(fig3, fullfile('figures','fig_recon.png'), 'Resolution', 250);
fprintf('Saved: figures/fig_recon.png\n');