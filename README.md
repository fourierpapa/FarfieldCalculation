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
