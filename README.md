# FarfieldCalculation
# 基于Matlab_根据全波仿真计算远场出射图
## 这个模型给出两个计算远场出射的方法
1. 自己使用FFT，考虑透镜的傅立叶变换，给出一个远场出射代码
  3. 夫琅禾费衍射本身带球面相位，而纯FFT不带球面相位，和直接FFT不一样，所以注释掉了夫琅禾费衍射
2. 计算bluesteinDFT代码，已经调通并能够运行
  1. 需要调用一个子程序

## 计算时需要给出，电场分量ExEyEz和坐标xy，
1. 使用fdtd给出计算结果的模型文件在./model/fdtd里
  1. 需要先运行全波仿真案例pc_fullsize_BIC.fsp
  2. 运行完结果后通过脚本ExportEfield.lsf输出探测器的
2. comsol文件太大没法放进来
### COMSOL运行完一个仿真后，
  1. 在Results-Datasets，右键More 3D Datasets-Parameterized Surface
  2. 在新产生的Datasets-Parameterized Surface 1中
    1. 选择数据集Dataset为求解完的结果
    2. Parameters里设定参数
      1. s1设置x的输出范围，s2设置y的输出范围
        1. 比如设置x范围是-p~p
    3. Expressions里设置
      1. x设s1，y设s2
      2. z设需要的平面的坐标
    4. Resolution里设置
      1. 分辨率，最好设个奇数比如201
  3. 选择一个画图模块 Results-Electric Field (ewfd)
    1. Dataset选择刚才设定的面，如Parameterized Surface 1
    2. 右键Electric Field (ewfd)选择surface，再Plot就可以绘图
    3. 绘制surface时，需要注意绘图给出ewfd.Ex只会给出实部
      1. 如果输出虚部，需要额外imag(ewfd.Ex)输出
    4. 右键选择surface，点击Add Plot Data to Export
  4. 在Results-Export中找到最新的Plot
    1. 在Output中点击Browse选择输出的文件夹和文件名
    2. 最后点击Settings下面的Export，就可以输出文件
  5. 再联合Matlab模型根据全波仿真计算远场出射图
