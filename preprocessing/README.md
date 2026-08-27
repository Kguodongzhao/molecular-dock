# 受体与配体自动预处理

本目录把原始分子文件转换为 `scripts/03_run_docking.py` 可直接使用的 PDBQT 文件。

环境安装和完整工作流参见[项目 README](../README.md)。

## 运行入口

所有命令建议从项目根目录运行：

```powershell
.\preprocessing\run.ps1 ligands --help
.\preprocessing\run.ps1 receptor --help
```

入口默认使用名为 `molecular-docking` 的 Conda 环境。

## 配体预处理

支持以下输入格式：

- SDF、SD
- MOL、MOL2
- PDB
- SMI、SMILES
- InChI

目录输入会递归扫描，并在输出中保留原有子目录结构。推荐将配体按对接分组存放：

```text
data/raw_ligands/
└── DEFAULT/
    ├── ligand_001.sdf
    └── ligand_002.mol2
```

运行：

```powershell
.\preprocessing\run.ps1 ligands `
    --input .\data\raw_ligands `
    --output-root .\data\ligands_pdbqt `
    --workers 4
```

默认处理过程：

1. Open Babel 按指定 pH 处理质子化状态；
2. RDKit 选择主要有机片段并生成或优化三维构象；
3. Meeko 分配电荷、可旋转键并输出 PDBQT；
4. 写入元数据和预处理汇总。

常用参数：

- `--ph 7.4`
- `--seed 42`
- `--keep-salts`
- `--keep-existing-3d`
- `--rigid-macrocycles`
- `--rigid-if-over 32`
- `--overwrite`
- `--keep-intermediate`

## 受体预处理

支持 PDB、ENT、CIF 和 mmCIF。默认选择第一个模型和 altloc A，并删除水、未指定保留的 HETATM 以及输入氢原子。

下面的对接盒坐标只是命令格式示例，必须替换为实际结合位点参数：

```powershell
.\preprocessing\run.ps1 receptor `
    --input .\data\raw_receptor\receptor.pdb `
    --output-root .\data\receptor `
    --name receptor `
    --box-center 10.0 20.0 30.0 `
    --box-size 20.0 20.0 20.0
```

也可以使用参考配体定义对接盒：

```powershell
.\preprocessing\run.ps1 receptor `
    --input .\data\raw_receptor\receptor.pdb `
    --output-root .\data\receptor `
    --name receptor `
    --box-enveloping .\data\raw_ligands\reference.sdf `
    --padding 5
```

常用参数：

- `--chains A,B`
- `--model 1`
- `--altloc A`
- `--keep-residue HEM`
- `--keep-hetero`
- `--keep-water`
- `--set-template A:42=HID`
- `--allow-bad-res`
- `--keep-intermediate`

主要输出：

- `receptor_prep.pdbqt`：Vina 使用的刚性受体；
- `receptor_prep.pdb`：构建和显示复合物；
- `receptor_metadata.json`：处理参数和校验信息；
- `receptor_docking_config_snippet.json`：配置路径片段；
- `receptor_vina_box.txt`：对接盒参数。

## 接入批量对接

复制配置模板：

```powershell
Copy-Item `
    .\config\docking.example.json `
    .\config\docking.local.json
```

修改本地配置后运行：

```powershell
python .\scripts\03_run_docking.py `
    --config .\config\docking.local.json `
    --input-root .\data\ligands_pdbqt `
    --output-root .\results\new_run `
    --workers 4
```

配置中的 `groups` 必须与配体 PDBQT 根目录下的分组文件夹一致。

## 科学注意事项

自动预处理不能替代人工结构检查。正式计算前应重点检查缺失残基或侧链、链选择、金属和辅因子、二硫键、质子化状态、配体手性及结合位点氢键网络。

不要仅为绕过结构错误而直接启用 `--allow-bad-res`。