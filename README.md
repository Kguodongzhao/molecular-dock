# Protein–Ligand Docking Workflow

一个可复现的蛋白质受体、配体预处理及 AutoDock Vina 批量对接流程，主要面向 Windows 和 Conda 环境。

本仓库不包含真实受体、配体、对接结果或软件二进制文件。

## 功能

- 使用 Open Babel、RDKit 和 Meeko 预处理配体；

- 使用 Meeko 预处理 PDB/mmCIF 受体；

- 配置对接盒、随机种子和 Vina 参数；

- 批量运行 AutoDock Vina；

- 生成最佳构象、复合物、元数据及 TSV/Markdown 汇总；

- 支持并行运行和断点续跑。

## 环境要求

已验证的核心版本：

- Python 3.12

- RDKit 2026.03

- Biopython 1.88

- Open Babel 3.1

- Meeko 0.7.1

- AutoDock Vina 1.2.7

创建环境：

```powershell

conda env create -f .\environment.yml

conda activate molecular-docking

```

AutoDock Vina 不随仓库分发。安装 Vina 后，可在当前 PowerShell 会话中指定其路径：

```powershell

$env:VINA_EXE = 'C:\path\to\vina.exe'

```

也可以把 `vina` 加入系统 `PATH`，或在本地配置中设置 `vina_relative_path`。

## 本地数据目录

所有分子结构和结果均由 `.gitignore` 排除，不会上传到 GitHub。推荐使用以下本地结构：

```text

data/

├── raw_ligands/

│   └── DEFAULT/

│       └── ligand.sdf

├── raw_receptor/

│   └── receptor.pdb

├── ligands_pdbqt/

│   └── DEFAULT/

│       └── ligand.pdbqt

└── receptor/

    ├── receptor_prep.pdb

    └── receptor_prep.pdbqt

```

`config` 中的分组名称必须与 `ligands_pdbqt` 下的文件夹名称一致。

## 创建本地配置

复制公开模板：

```powershell

Copy-Item `

    .\config\docking.example.json `

    .\config\docking.local.json

```

编辑 `config\docking.local.json`，至少确认：

- 受体 PDBQT 和显示用 PDB 路径；

- 配体 PDBQT 根目录；

- `groups` 分组；

- 对接盒中心 `center`；

- 对接盒尺寸 `size`；

- `exhaustiveness`、`num_modes` 和 CPU 参数。

模板中的对接盒数值只是占位值，不能直接用于正式计算。

## 配体预处理

如果原始配体按 `DEFAULT` 等分组存放，脚本会保留子目录结构：

```powershell

.\preprocessing\run.ps1 ligands `

    --input .\data\raw_ligands `

    --output-root .\data\ligands_pdbqt `

    --workers 4

```

查看全部参数：

```powershell

.\preprocessing\run.ps1 ligands --help

```

## 受体预处理

下面的对接盒坐标仅为命令格式示例，必须替换成实际结合位点坐标：

```powershell

.\preprocessing\run.ps1 receptor `

    --input .\data\raw_receptor\receptor.pdb `

    --output-root .\data\receptor `

    --name receptor `

    --box-center 10.0 20.0 30.0 `

    --box-size 20.0 20.0 20.0

```

也可以根据参考配体生成对接盒：

```powershell

.\preprocessing\run.ps1 receptor `

    --input .\data\raw_receptor\receptor.pdb `

    --output-root .\data\receptor `

    --name receptor `

    --box-enveloping .\data\raw_ligands\reference.sdf `

    --padding 5

```

查看全部参数：

```powershell

.\preprocessing\run.ps1 receptor --help

```

## 运行对接

```powershell

python .\scripts\03_run_docking.py `

    --config .\config\docking.local.json `

    --workers 4

```

断点续跑：

```powershell

python .\scripts\03_run_docking.py `

    --config .\config\docking.local.json `

    --workers 4 `

    --resume

```

默认结果写入 `results/current_run/`。主要输出包括：

- `results.tsv`

- `SUMMARY.md`

- 全部对接构象

- 最佳构象

- 最佳受体–配体复合物

- 每个任务的元数据和 Vina 日志

## 科学使用注意事项

自动预处理不能替代人工结构检查。正式对接前应检查：

- 缺失残基和侧链；

- 蛋白链、替代构象和模型选择；

- 金属、辅因子、结晶水和其他异质残基；

- His、Asp、Glu、Lys 等残基的质子化状态；

- 二硫键和结合位点氢键网络；

- 配体质子化、互变异构、手性和电荷状态；

- 对接盒是否覆盖真实结合位点。

对接分数适合用于同一工作流中的相对比较，不应直接解释为实验结合自由能。

## 隐私检查

提交前始终运行：

```powershell

git status --short

git check-ignore -v -- path\to\private_structure.pdb

```

确认分子结构、本地配置、结果和二进制文件均未进入 Git。
