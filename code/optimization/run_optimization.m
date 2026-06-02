%% run_optimization.m
% 优化流水线 v3.1 入口脚本 — 多目标 (TF+位移+隔离带) + NIC负参数搜索
% 用法：在 MATLAB 中运行此脚本即可完成全部优化流程
%   >> run_optimization
%
% 目标函数: J = w1*TF_peak + w2*I1(位移) + w3*I3(0dB隔离) + w4*I4(-40dB) + pen_stab
%
% 流程：
%   Phase 1: HBM/弧长延拓快筛 (800 样本, Budget=1500, 含负参数)
%   Phase 2: 系统层 HBM FRF 精细验证 (Top-30, Budget=3000)
%   Phase 3: fminsearch 局部精修 (Top-8, sigmoid 边界约束)
%   Phase 4: 候选验证 + Wang BG 基线 + 算子解释 + 出版级图表
%   Phase 5: 自动保存结果 .mat + 导出 PDF 图片
%
% 搜索边界: sigma∈[-3,3], kap_e∈[0.02,3], kap_c∈[-3,3]
% 预计运行时间: ~30-60 分钟 (取决于 CPU)

% 添加路径
code_root = fileparts(fileparts(mfilename('fullpath')));
addpath(code_root);
init_path;

% 运行主优化脚本
run(fullfile(fileparts(mfilename('fullpath')), 'unified_optimization.m'));
