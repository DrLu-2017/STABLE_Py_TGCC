% =========================================================================
% plot_scan_map.m
% 
% 目的：读取 STABLE 扫描输出的 scan_summary.mat，绘制高分辨率 2D 稳定性图谱，
%       包含完整色标、理论平势阱曲线（xi=1.0）、以及每个仿真点的离散形状标记：
%         o  (圆形)  : Stable (稳定)
%         x  (叉形)  : PTBL (周期瞬态束长抖动)
%         ★  (五角星): Coupled D-Q (偶极-四极耦合失稳)
%         ▲  (三角形): LCBI (纵向耦合束团失稳)
%         ◆  (菱形)  : Other / Unknown
% =========================================================================

clc; clear; close all;

%% 1. 查找或指定扫描结果路径
% 自动寻找 scan_data 目录下包含 scan_summary.mat 的最新扫描文件夹
d_list = dir('scan_data/Scan_*');
if isempty(d_list)
    error('未找到任何扫描结果文件夹，请确认 scan_data 路径。');
end
[~, s_idx] = sort([d_list.datenum], 'descend');
scan_file = '';
for k = 1:length(s_idx)
    cand = fullfile(d_list(s_idx(k)).folder, d_list(s_idx(k)).name, 'scan_summary.mat');
    if exist(cand, 'file')
        scan_file = cand;
        break;
    end
end
if isempty(scan_file)
    error('未找到任何包含 scan_summary.mat 的扫描文件夹。');
end

fprintf('正在加载扫描数据: %s ...\n', scan_file);
data = load(scan_file);

Df_scan_kHz    = data.Df_scan_kHz;
I0_scan_mA     = data.I0_scan_mA;
stability_code = data.stability_code;

%% 2. 创建独立绘图窗口
fig = figure(101); clf;
set(fig, 'Color', 'w', 'Position', [150 150 900 620]);
ax = axes('Parent', fig);
hold(ax, 'on');

% 绘制背景 2D 热图
imagesc(ax, Df_scan_kHz, I0_scan_mA, stability_code);
set(ax, 'YDir', 'normal');

% 定义与 mbtrack2 对应的专属配色方案
cmap = [0.20 0.63 0.17;   % 0: 绿色 (Stable)
        0.99 0.75 0.05;   % 1: 琥珀金 (PTBL)
        0.85 0.33 0.10;   % 2: 橙红色 (Coupled D-Q)
        0.55 0.10 0.65;   % 3: 深紫色 (LCBI)
        0.50 0.50 0.50];  % 4: 灰色 (Other)

colormap(ax, cmap);
clim(ax, [-0.5 4.5]);

% 配置右侧颜色条 (Colorbar)
cb = colorbar(ax);
cb.Ticks = 0:4;
cb.TickLabels = {'Stable (稳定)', 'PTBL (束长振荡)', 'Coupled D-Q (偶极四极耦合)', 'LCBI (耦合束团失稳)', 'Other / Unknown'};
cb.Label.String = 'Instability Classification';
cb.Label.FontSize = 11;
cb.Label.FontName = 'Times New Roman';

%% 3. 绘制理论平势阱曲线 (\xi = 1.0)
try
    C_val   = 354.097; h_val = 416; n_hc = 4;
    U0_val  = 470e3;   V_mc  = 1.7e6;
    R_hc    = 1.86e6;  Q_hc  = 31e3;
    if isfield(data, 'config')
        if isfield(data.config, 'C'), C_val = data.config.C; end
        if isfield(data.config, 'h'), h_val = round(data.config.h); end
        if isfield(data.config, 'n_hc'), n_hc = round(data.config.n_hc); end
        if isfield(data.config, 'U0'), U0_val = data.config.U0; end
        if isfield(data.config, 'V_mc'), V_mc = data.config.V_mc; end
        if isfield(data.config, 'R_hc'), R_hc = data.config.R_hc; end
        if isfield(data.config, 'Q_hc'), Q_hc = data.config.Q_hc; end
    end
    I0_line = linspace(min(I0_scan_mA), max(I0_scan_mA), 100) * 1e-3;
    Df_opt  = arrayfun(@(Ib) detune_HC_calc(Ib, n_hc, C_val, h_val, U0_val, V_mc, R_hc, Q_hc), I0_line) * 1e-3;
    plot(ax, Df_opt, I0_line * 1e3, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Flat-pot. (\xi=1.0)');
catch ME_opt
    fprintf('提示：未绘制平势阱参考线（%s）\n', ME_opt.message);
end

%% 4. 在每个仿真网格点上绘制不同形状的标记
[Df_grid, I0_grid] = meshgrid(Df_scan_kHz, I0_scan_mA);

% 用于生成清晰图例的虚拟图层
plot(ax, nan, nan, 'o', 'MarkerFaceColor', [0.20 0.70 0.20], 'MarkerEdgeColor', [0 0.30 0], 'MarkerSize', 7, 'DisplayName', 'Stable (o)');
plot(ax, nan, nan, 'x', 'Color', [0.85 0.50 0], 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'PTBL (x)');
plot(ax, nan, nan, 'p', 'MarkerFaceColor', [0.90 0.20 0.10], 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'DisplayName', 'Coupled D-Q (★)');
plot(ax, nan, nan, '^', 'MarkerFaceColor', [0.60 0.10 0.70], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', 'LCBI (▲)');
plot(ax, nan, nan, 'd', 'MarkerFaceColor', [0.50 0.50 0.50], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Other (◆)');

legend(ax, 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman', 'Box', 'on', 'Color', [1 1 1 0.9]);

% 绘制真实网格点标记
m0 = (stability_code == 0);
m1 = (stability_code == 1);
m2 = (stability_code == 2);
m3 = (stability_code == 3);
m4 = (stability_code == 4);

if any(m0(:)), plot(ax, Df_grid(m0), I0_grid(m0), 'o', 'MarkerFaceColor', [0.20 0.70 0.20], 'MarkerEdgeColor', [0 0.30 0], 'MarkerSize', 6, 'HandleVisibility', 'off'); end
if any(m1(:)), plot(ax, Df_grid(m1), I0_grid(m1), 'x', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'HandleVisibility', 'off'); end
if any(m2(:)), plot(ax, Df_grid(m2), I0_grid(m2), 'p', 'MarkerFaceColor', [0.90 0.20 0.10], 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'HandleVisibility', 'off'); end
if any(m3(:)), plot(ax, Df_grid(m3), I0_grid(m3), '^', 'MarkerFaceColor', [0.60 0.10 0.70], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility', 'off'); end
if any(m4(:)), plot(ax, Df_grid(m4), I0_grid(m4), 'd', 'MarkerFaceColor', [0.50 0.50 0.50], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off'); end

%% 5. 格式化坐标轴
title(ax, '2D Longitudinal Stability Map (STABLE Simulation)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel(ax, 'HHC Detuning \Delta f_{HC}  [kHz]', 'FontSize', 12);
ylabel(ax, 'Beam Current I_0  [mA]', 'FontSize', 12);

dDf = 2.5; if length(Df_scan_kHz) > 1, dDf = (Df_scan_kHz(2)-Df_scan_kHz(1))/2; end
dI0 = 12.5; if length(I0_scan_mA) > 1, dI0 = (I0_scan_mA(2)-I0_scan_mA(1))/2; end
xlim(ax, [min(Df_scan_kHz)-dDf, max(Df_scan_kHz)+dDf]);
ylim(ax, [min(I0_scan_mA)-dI0, max(I0_scan_mA)+dI0]);

grid(ax, 'on');
set(ax, 'FontSize', 11, 'FontName', 'Times New Roman', 'Box', 'on');

hold(ax, 'off');

% 保存图片
out_png = 'stability_map_with_markers.png';
saveas(fig, out_png);
fprintf('失谐频率图谱已保存至: %s\n', out_png);

%% =========================================================================
%% 6. 创建调谐角 \psi_{HC} 图谱（精确对标 mbtrack2）
%% =========================================================================
cspeed  = 299792458;
f0_val  = cspeed / C_val;
f_rf    = h_val * f0_val;
f_hc    = n_hc * f_rf;

% 转换为调谐角 \psi_{HC} [deg]
psi_scan_deg = atand(2 * Q_hc * (Df_scan_kHz * 1e3) / f_hc);

fig2 = figure(102); clf;
set(fig2, 'Color', 'w', 'Position', [180 180 900 620]);
ax2 = axes('Parent', fig2);
hold(ax2, 'on');

% 绘制非线性网格背景色块
d_psi = diff(psi_scan_deg);
psi_edges = [psi_scan_deg(1) - d_psi(1)/2, ...
             (psi_scan_deg(1:end-1) + psi_scan_deg(2:end))/2, ...
             psi_scan_deg(end) + d_psi(end)/2];
d_I0 = diff(I0_scan_mA);
I0_edges = [I0_scan_mA(1) - d_I0(1)/2, ...
            (I0_scan_mA(1:end-1) + I0_scan_mA(2:end))/2, ...
            I0_scan_mA(end) + d_I0(end)/2];

[X_mesh, Y_mesh] = meshgrid(psi_edges, I0_edges);
pcolor(ax2, X_mesh, Y_mesh, [stability_code, stability_code(:,end); stability_code(end,:), stability_code(end,end)]);
shading(ax2, 'flat');
colormap(ax2, cmap);
clim(ax2, [-0.5 4.5]);

cb2 = colorbar(ax2);
cb2.Ticks = 0:4;
cb2.TickLabels = {'Stable (稳定)', 'PTBL (束长振荡)', 'Coupled D-Q (偶极四极耦合)', 'LCBI (耦合束团失稳)', 'Other / Unknown'};
cb2.Label.String = 'Instability Classification';
cb2.Label.FontSize = 11;
cb2.Label.FontName = 'Times New Roman';

% 理论平势阱曲线 (\xi = 1.0) 及 PTBL 理论边界 (\xi = 1.05) vs \psi
k_fp = sqrt(1/n_hc^2 - 1/(n_hc^2-1) * (U0_val/V_mc)^2);
cos_psi_opt = (k_fp * V_mc) ./ (2 * I0_line * R_hc);
cos_psi_opt = max(-1, min(1, cos_psi_opt));
psi_opt_line = acosd(cos_psi_opt);
plot(ax2, psi_opt_line, I0_line * 1e3, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Flat-pot. (\xi=1.0)');

cos_psi_105 = (1.05 * k_fp * V_mc) ./ (2 * I0_line * R_hc);
cos_psi_105 = max(-1, min(1, cos_psi_105));
psi_105_line = acosd(cos_psi_105);
plot(ax2, psi_105_line, I0_line * 1e3, 'g:', 'LineWidth', 1.8, 'DisplayName', '\xi=1.05 (PTBL thresh.)');

% 绘制离散形状标记
[Psi_grid, ~] = meshgrid(psi_scan_deg, I0_scan_mA);
plot(ax2, nan, nan, 'o', 'MarkerFaceColor', [0.20 0.70 0.20], 'MarkerEdgeColor', [0 0.30 0], 'MarkerSize', 7, 'DisplayName', 'Stable (o)');
plot(ax2, nan, nan, 'x', 'Color', 'k', 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'PTBL (x)');
plot(ax2, nan, nan, 'p', 'MarkerFaceColor', [0.90 0.20 0.10], 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'DisplayName', 'Coupled D-Q (★)');
plot(ax2, nan, nan, '^', 'MarkerFaceColor', [0.60 0.10 0.70], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', 'LCBI (▲)');

if any(m0(:)), plot(ax2, Psi_grid(m0), I0_grid(m0), 'o', 'MarkerFaceColor', [0.20 0.70 0.20], 'MarkerEdgeColor', [0 0.30 0], 'MarkerSize', 6, 'HandleVisibility', 'off'); end
if any(m1(:)), plot(ax2, Psi_grid(m1), I0_grid(m1), 'x', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'HandleVisibility', 'off'); end
if any(m2(:)), plot(ax2, Psi_grid(m2), I0_grid(m2), 'p', 'MarkerFaceColor', [0.90 0.20 0.10], 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'HandleVisibility', 'off'); end
if any(m3(:)), plot(ax2, Psi_grid(m3), I0_grid(m3), '^', 'MarkerFaceColor', [0.60 0.10 0.70], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility', 'off'); end
if any(m4(:)), plot(ax2, Psi_grid(m4), I0_grid(m4), 'd', 'MarkerFaceColor', [0.50 0.50 0.50], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off'); end

legend(ax2, 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman', 'Box', 'on', 'Color', [1 1 1 0.9]);
title(ax2, '2D Longitudinal Stability Map vs HC Tuning Angle \psi_{HC} (SOLEIL II)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel(ax2, 'HC Tuning Angle \psi_{HC}  [^\circ]', 'FontSize', 12);
ylabel(ax2, 'Beam Current I_0  [mA]', 'FontSize', 12);
xlim(ax2, [65.0, 85.0]);
ylim(ax2, [min(I0_scan_mA)-12.5, max(I0_scan_mA)+12.5]);
grid(ax2, 'on'); set(ax2, 'FontSize', 11, 'FontName', 'Times New Roman', 'Box', 'on');
hold(ax2, 'off');

out_psi_png = 'stability_map_psi.png';
saveas(fig2, out_psi_png);
fprintf('调谐角图谱已保存至: %s\n', out_psi_png);

