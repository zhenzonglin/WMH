# 工作站版本发布验证

验证日期：2026-09-13。环境：Ubuntu-20.04 / Python 3.11.13。此记录对应新增的工作站接入流程；既有统计数值验证见`validation_report.md`。

- `ruff check src tests`通过。
- `pytest -q`：48项通过，包括原有34项和14项新增检查。
- 公共SAS文件验证未压缩、压缩、日期时间和特殊缺失的读取；不是CNSR-III患者文件。
- 新检查覆盖字段存在但全部缺失、非法编码、重复ID、前导零、大小写、来源歧义、数值ID特殊缺失、影像精确连接、QC未审核、两步配置以及停止阶段。
- 700例合成数据通过`run --through prepare`，完成临床审计、CSV入口、影像连接和病例准备。
- 同一合成数据通过`run --through report --hypothesis H1`，验证从审计到H1模型、敏感性分析及HTML报告的完整命令入口。
- 已测试临床审计无需影像目录、审计失败停止后续步骤、默认准备阶段不调用统计模型。

未在本机读取真实患者SAS或连接真实WMH目录。实际编码、文件间覆盖、患者数、影像路径布局和模型可估计性，在工作站接入后依据审计结果确定。上述通过状态仅代表软件检查。

## Conda兼容性验证（2026-09-13）

- 使用Conda 25.7.0按`environment.yml`从头创建独立Linux环境，Python 3.11.16；原有uv环境保留。
- 33个适用于Linux的Python依赖版本与`uv.lock`一致，`pip check`通过。
- 在激活的Conda环境中，`ruff check src tests`通过，`python -m pytest -q`为48项通过。
- 直接调用`wmh-hcy`命令成功；700例合成数据完成`run --through prepare`，不需要`uv run`。
- 环境文件、固定依赖清单和操作说明已纳入仓库；验证环境和数据输出不纳入Git。
