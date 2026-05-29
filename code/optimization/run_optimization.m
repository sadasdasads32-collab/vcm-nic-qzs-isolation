%% run_optimization.m
% 两级优化流水线入口脚本
% 用法：在 MATLAB 中运行此脚本即可完成全部优化流程
%   >> run_optimization
%
% 流程：
%   Phase 1: 算子层 LHS 快筛 (5000 样本)
%   Phase 2: 系统层 HBM FRF 验证 (Top-30)
%   Phase 3: fminsearch 精修 (Top-8)
%   Phase 4: Floquet 稳定性全扫描 + 出版级对比图
%   Phase 5: 自动保存结果 .mat + 导出 PDF 图片
%
% 预计运行时间: 1-2 小时 (取决于 CPU)

% 添加路径
code_root = fileparts(fileparts(mfilename('fullpath')));
addpath(code_root);
init_path;

% 运行主优化脚本
run(fullfile(fileparts(mfilename('fullpath')), 'unified_optimization.m'));
