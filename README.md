# 20220624-virtual-clustering
本项目用于对小分子库进行降重
已经完成分子指纹计算+并行化

文件目录:
```
- data_dir
  - AA
    - AABAAC.tar
    - AABABC.tar
    ...
  - AB
  ...
- data.tools
  - pipeline.py
  - utils.py
- scripts
  - rdkit2pdbqt.py
- parallel_test.py
```

用法：python .\parallel_test.py
参数整合在这个文件中了，之后会有一个完整说明的程序，在run.py中。
mid_format不要改（已废弃），可以改fp_type, n_cpu, data_dir, out_dir
```python
out_dir = './debug/output'
data_pipeline = pipeline.DataPipeline(
    data_dir='./debug/example',
    n_cpu=16,
    mid_format='pdbqt',
    fp_type='maccs',
    out_dir=out_dir
)
```
可以改delete=False参数来保留临时文件夹（用于解压）
```python
with utils.tmpdir_manager('.', delete=True) as tmp_dir:
```