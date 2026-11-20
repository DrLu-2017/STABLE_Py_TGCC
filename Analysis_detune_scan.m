% =========================================================================
% Analysis_detune_scan.m (全自动化动态扫描分析版)
% 目的：自动检索当前目录下的失谐量数据，分析拉伸效果并寻找最优工作点
% =========================================================================
clc; clear; close all;

fprintf('正在自动检索当前目录下的扫描数据...\n');
fprintf('--------------------------------------------------\n');

% 1. 自动寻找当前目录下所有符合 Soleil 2 命名规则的 .mat 文件
file_pattern = 'SOLEIL_I0500mA_RLfp930000_QLfp31000_detune*_kp1_ki1e-05_1.mat';
dir_info = dir(file_pattern);

if isempty(dir_info)
    % 备用匹配：防止最后的数字不是 1 或者有微调
    file_pattern = 'SOLEIL_I0*_RLfp*_QLfp*_detune*.mat';
    dir_info = dir(file_pattern);
end

if isempty(dir_info)
    error(['错误：在当前目录下未找到任何符合命名规则的扫描数据文件！\n' ...
           '请检查当前工作路径下是否存在形如 SOLEIL_I0500mA... 的文件。']);
end

% 2. 从文件名中动态解析出所有的失谐量 (detune) 数值
num_files = length(dir_info);
fre_shift_scan = zeros(1, num_files);
valid_file_flags = true(1, num_files); % 用于过滤未能匹配成功的文件

for f = 1:num_files
    fname = dir_info(f).name;
    % 运用正则表达式提取 detune 后面的数字（包含正负号和科学计数法，如 detune-50000 或 detune50000）
    tokens = regexp(fname, 'detune(-?\d+\.?\d*[eE]?[-+]?\d*)', 'tokens');
    if ~isempty(tokens)
        fre_shift_scan(f) = str2double(tokens{1}{1});
    else
        valid_file_flags(f) = false; % 如果提取失败，标记为无效
    end
end

% 过滤并重新排序，确保曲线绘制连续性
dir_info = dir_info(valid_file_flags);
fre_shift_scan = fre_shift_scan(valid_file_flags);

[fre_shift_scan, sort_idx] = sort(fre_shift_scan);
dir_info = dir_info(sort_idx);

num_scans = length(fre_shift_scan);
mean_bunch_length_ps = zeros(1, num_scans);  % 稳定态平均束长
centroid_jitter_ps = zeros(1, num_scans);    % 质心抖动

% 3. 循环读取文件并提取物理数据
for i = 1:num_scans
    df_khz = fre_shift_scan(i) / 1e3;
    filename = dir_info(i).name;
    
    % 载入数据
    data = load(filename);
    
    % 动态获取主程序在 HALF 结构体中保存的自然束长 sigma_t0
    if isfield(data, 'HALF') && isfield(data.HALF, 'sigma_t0')
        sigma_t0 = data.HALF.sigma_t0;
    else
        sigma_t0 = 40e-12; % 兜底默认值
    end
    
    % 取最后 10% 的追踪圈数（稳定态）进行统计
    num_records = size(data.record_Q_std, 1);
    stable_range = round(num_records*0.9):num_records;
    
    % 计算所有束团在稳定态下的平均束长 (单位: ps)
    mean_bunch_length_ps(i) = mean(mean(data.record_Q_std(stable_range, :))) * sigma_t0 * 1e12;
    
    % 计算稳定态下的束团质心位置的整体标准差
    centroid_jitter_ps(i) = mean(std(data.record_Q_mean(stable_range, :))) * sigma_t0 * 1e12;
    
    fprintf('找到文件: %-75s | 失谐量: %4.0f kHz | 束长: %6.2f ps\n', ...
        filename, df_khz, mean_bunch_length_ps(i));
end

% 移除可能解析失败的无效数据点
valid_idx = ~isnan(mean_bunch_length_ps) & (fre_shift_scan ~= 0 | mean_bunch_length_ps ~= 0);
scan_plot = fre_shift_scan(valid_idx) / 1e3;
bl_plot = mean_bunch_length_ps(valid_idx);
jitter_plot = centroid_jitter_ps(valid_idx);

% 4. 寻找最优工作点
stability_threshold = 1.5; % ps (束流纵向质心抖动的安全阈值，可根据实际情况微调)
stable_indices = find(jitter_plot < stability_threshold);

if ~isempty(stable_indices)
    [max_bl, best_sub_idx] = max(bl_plot(stable_indices));
    best_idx = stable_indices(best_sub_idx);
    best_detune = scan_plot(best_idx);
    best_bunch_length = bl_plot(best_idx);
else
    [~, best_idx] = min(jitter_plot);
    best_detune = scan_plot(best_idx);
    best_bunch_length = bl_plot(best_idx);
    warning('所有扫描工作点均伴随一定程度的纵向震荡！');
end

fprintf('--------------------------------------------------\n');
fprintf('【自动分析结论】当前文件夹下共有 %d 个有效扫描点。\n', length(scan_plot));
fprintf('最优工作点处于 HHC 失谐量: **%.0f kHz**，此时稳定束长为 **%.2f ps**。\n', ...
    best_detune, best_bunch_length);

% 5. 绘图展示
figure('Position', [200, 200, 850, 500]);
yyaxis left
plot(scan_plot, bl_plot, '-o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
ylabel('Bunch Length \sigma_{\tau} (ps)', 'FontSize', 12);
hold on;

plot(best_detune, best_bunch_length, 'p', 'MarkerSize', 14, ...
     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(best_detune + 5, best_bunch_length, sprintf('Optimal Point (%.0f kHz)', best_detune), ...
     'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');

yyaxis right
plot(scan_plot, jitter_plot, '--s', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
ylabel('Centroid Jitter (ps)', 'FontSize', 12);
yline(stability_threshold, ':r', 'Stability Limit', 'LineWidth', 1.2);

xlabel('HHC Detuning \Delta f_{hc} (kHz)', 'FontSize', 12);
title('Bunch Lengthening effect vs. HHC Detuning (Automated Survey)', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);