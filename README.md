# 人脑皮层空间转录组分析

陈辰的 Python 健康数据科学课程项目。公开展示页：https://ccipad123.github.io/spatial-transcriptomics/

基于 Maynard et al. (2021) 的 DLPFC Visium 数据，使用同一供体的相邻切片 151673、151674，展示空间表达、跨切片预测与误差分析。这里是已计算结果的静态展示；网页不重新训练模型。

## 页面内容

- 两张组织切片与参考层、12 个基因表达图层。
- 留出切片的预测层、错误位置、预测熵及捕获点信息。
- 25 组矢量分析图、模型指标 CSV、分析代码与 35 篇参考文献。
- 课程报告公开版 PDF；内容与提交版一致，封面不含学号。

## 维护与发布

GitHub Pages 使用 `main` 分支的 `/docs`。页面、样式、脚本和静态数据均位于 `docs/`，无需安装构建依赖。

可在仓库目录运行 `python3 -m http.server 8769 --bind 127.0.0.1 --directory docs` 预览。更新后检查电脑/手机、切片与基因切换、图版弹窗及报告下载，再提交推送。子路径为 `/spatial-transcriptomics/`，资源使用相对地址；不要写成本机路径。

原始数据和完整分析工程在本地课程项目保存。本地 `src/prepare_project_github.py` 可从已核验的静态页面和公开 PDF 更新本仓库。不要将项目记忆、学号、认证配置或原 Sites 仓库历史复制进来。

## 数据与解释范围

原始数据研究：[Maynard et al., Nature Neuroscience (2021)](https://doi.org/10.1038/s41593-020-00787-0)。
数据来源：[LieberInstitute/HumanPilot](https://github.com/LieberInstitute/HumanPilot)。

H&E 为原作者提供的 2000×2000 缩小图。组织图与原数据权利及许可遵循原始来源；本仓库未将其重新授权。两张切片同属一名供体，参考标签含表达证据，结果只用于教学复分析，不能代表跨供体泛化能力。
