# 20220831-virtual-clustering
- [x] 完成tar包解压 进程池并行
- [x] 完成指纹计算base64位存储 mpi并行，分成了 576 block，每个200w分子(200M的存储空间)
- [ ] 读取base64指纹到numpy数组
- [ ] 分层kmeans聚类

当前问题：
读取base64到numpy数组出现OOM问题
sklearn的多节点计算支持，计划采用dask

