# PyMOL Siteview Tool

一个纯 Python 的 PyMOL 自动作图工具，用于快速读取 `pdb`、`cif/mmcif`、`mol2`、`sdf/mol` 等结构文件，生成蛋白-小分子结合位点视角，保存 `.pse` 会话文件，并可同步输出 `.png` 图片。

## 功能

- 自动加载结构文件并识别蛋白、小分子配体。
- 自动生成 ligand-site 的 siteview 视角。
- 如果结构中有小分子：
  - 蛋白 cartoon transparency 默认为 `0.8`。
  - 小分子和结合位点残基以 sticks 显示，stick radius 默认为 `0.2`。
- 自动绘制相互作用：
  - 氢键 / polar contacts：黄色 dash line。
  - pi-pi 相互作用：蓝色 dash line。
  - pi-cation 相互作用：绿色 dash line。
- 使用 Arial-like 的 PyMOL label 字体设置标记相关氨基酸。
- 输出 `*_siteview.pse` 和 `*_siteview.png`。

## 安装

推荐使用 conda，因为 PyMOL 的 Python API 最稳定的安装方式通常来自 conda-forge。

```bash
git clone https://github.com/your-name/pymol-siteview.git
cd pymol-siteview
conda env create -f environment.yml
conda activate pymol-siteview
```

安装后检查命令行入口：

```bash
pymol-siteview --help
```

也可以使用安装脚本：

```bash
bash scripts/install_conda.sh
```

Windows PowerShell：

```powershell
.\scripts\install_conda.ps1
```

## 使用

安装后推荐使用：

```bash
pymol-siteview input.pdb -o siteview_out
```

也可以直接用 PyMOL 运行仓库中的兼容脚本：

```bash
pymol -cq pymol_siteview.py -- input.pdb -o siteview_out
```

批量处理：

```bash
pymol-siteview *.pdb *.cif *.mol2 -o siteview_out
```

如果自动识别的小分子不正确，可以手动指定 ligand：

```bash
pymol-siteview complex.pdb -o siteview_out --ligand-selection "resn LIG"
```

## 常用参数

```bash
pymol-siteview input.pdb \
  -o siteview_out \
  --ligand-selection "resn LIG" \
  --site-cutoff 4.5 \
  --cartoon-transparency 0.8 \
  --stick-radius 0.2
```

- `--ligand-selection`：手动指定配体选择，例如 `"resn LIG"`、`"chain A and resi 501"`。
- `--site-cutoff`：结合位点残基距离阈值，默认 `4.5` A。
- `--hbond-cutoff`：氢键 / polar contacts 阈值，默认 `3.6` A。
- `--pi-cutoff`：pi-pi 环中心距离阈值，默认 `5.5` A。
- `--pi-cation-cutoff`：pi-cation 距离阈值，默认 `6.0` A。
- `--no-png`：只保存 `.pse`，不渲染图片。
- `--no-ray`：PNG 使用 OpenGL 截图，不做 ray tracing，速度更快。

## 输出

默认输出到 `siteview_out/`：

```text
siteview_out/
  input_siteview.pse
  input_siteview.png
```

`.pse` 文件可直接用 PyMOL 打开并继续手动调整。

## 代码结构

```text
pymol-siteview/
  src/pymol_siteview/cli.py      # 主程序
  pymol_siteview.py              # 兼容 PyMOL 直接运行的入口
  environment.yml                # conda 环境
  pyproject.toml                 # Python 包配置
  scripts/                       # 安装脚本
  examples/                      # 使用示例
  tests/                         # 不依赖 PyMOL 的基础测试
```

## 开发检查

基础测试不依赖 PyMOL：

```bash
python -m pytest
python -m compileall src pymol_siteview.py
```

真实渲染需要 PyMOL：

```bash
pymol-siteview your_structure.pdb -o siteview_out
```

## 注意

`.pse` 是 PyMOL 会话文件，所以必须在带 PyMOL Python API 的环境中运行。pi-pi 和 pi-cation 是基于几何规则的自动近似识别；如果结构缺少键级、氢或电荷信息，建议使用 `--ligand-selection` 手动指定配体，并在 PyMOL 中复核关键相互作用。
# pymol-siteview
