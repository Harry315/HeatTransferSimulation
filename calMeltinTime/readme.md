# 利用瞬态热传递进行石蜡包裹合金结构的熔化温度仿真

## 使用流程：
stp2Mesh4thermoMix.py 输入：（合金骨架）stp文件 将石蜡与合金骨架结合，赋予不同的材料属性，并生成八分之一对称网格 输出：inp文件
calculateMeltingTime4symmetryMesh.py 输入：inp文件（上一步生成） 施加分析步、边界条件、提交计算 输出：odb文件
odb2TimeTemp.py 输入：odb文件（上一步生成） 选择合金骨架/填充石蜡/任意点集，导出时间-平均温度关系 输出：np数组
postProcess.py 输入：np数组 选择插值方法和阈值，绘制曲线、确定熔化时间。

## 具体脚本

### stp2Mesh4thermoMix.py
### calculateMeltingTime4symmetryMesh.py
### odb2TimeTemp.py
### postProcess.py