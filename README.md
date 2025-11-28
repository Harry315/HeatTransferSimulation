# Heat Transfer Simulation

## 使用流程
1. 用户给出**stp**文件
2. 运行脚本generateMeshArray4homo，输入stp文件，生成对称网格，输出inp文件
3. 运行脚本calculateHFL4symmetryMesh，输入inp文件，设置物理条件，进行传热仿真计算，输出odb文件
4. 运行脚本odb2thermoConductivity，输入odb文件，输出选定表面点集的平均热流密度

## 脚本详细说明

### generateMeshArray4homo

### calculateHFL4symmetryMesh

### odb2thermoConductivity