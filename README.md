# RANSAC

本项目分别用python和C++实现了ransac算法，并且通过ransac算法对数据集中的点云进行平面拟合。

## RANSAC简介

**RANSAC**(RAndom SAmple Consensus,随机采样一致)算法是从一组含有“外点”(outliers)的数据中正确估计数学模型参数的迭代算法。“外点”一般指的的数据中的噪声，比如说匹配中的误匹配和估计曲线中的离群点。所以，RANSAC也是一种“外点”检测算法。RANSAC算法是一种不确定算法，它只能在一种概率下产生结果，并且这个概率会随着迭代次数的增加而加大（之后会解释为什么这个算法是这样的）。RANSAC算最早是由Fischler和Bolles在SRI上提出用来解决LDP(Location Determination Proble)问题的。

**ransac的基本假设**
整个数据集中同时包含好的点和不好的点，我们将它们分别称为局内点和局外点；
数据的分布可以通过某个数学模型来描述，而局内点就是可以适应该模型的点，局外点是不可以适应该模型的点；
随意给定一组点（可能是局内点也有可能是局外点），我们假设这组点都是局内点，它们满足某个数学模型，我们利用这个数学模型去估算其他的点，如果有足够多的点通过计算能被归类为假设的局内点，那么这个模型就足够合理。

## 算法基本思想和流程

RANSAC是通过反复选择数据集去估计出模型，一直迭代到估计出认为比较好的模型。
具体的实现步骤可以分为以下几步：

1. 选择出可以估计出模型的最小数据集；(对于直线拟合来说就是两个点，对于计算Homography矩阵就是4个点)
2. 使用这个数据集来计算出数据模型；
3. 将所有数据带入这个模型，计算出“内点”的数目；(累加在一定误差范围内的适合当前迭代推出模型的数据)
4. 比较当前模型和之前推出的最好的模型的“内点“的数量，记录最大“内点”数的模型参数和“内点”数；
5. 重复1-4步，直到迭代结束或者当前模型已经足够好了(“内点数目大于一定数量”)。