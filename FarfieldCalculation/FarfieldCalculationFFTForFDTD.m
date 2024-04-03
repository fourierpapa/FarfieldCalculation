% function fdtdFarfieldCalculation
% 2024年3月25日 由fyb创建
% 功能：根据从FDTD输出的电场，给出辐射远场图
% 特点：使用了夫琅禾费衍射
clear all;
close all;
clc;

% load functions
addpath(genpath('./utils'))
addpath(genpath('../src'))

%% 定义参数
% physical parameters
% 载入待计算的场
load('E_Field.mat')
xx = E.x;
yy = E.y;

params.pxsize = xx(2)-xx(1);            % pixel size (m)
params.wavlen = 0.6e-6;                 % wavelength (m)
params.dist   = 5e-3;                      % imagicng distance (m)
k0 = 2*pi/params.wavlen;

kernelsize = 4;
nullpixels = ceil(kernelsize / params.pxsize);                             % number of padding pixels
nullpixels = 8000;

f_FT = 5e-3;                      % lens focal length (m)
NA = 0.8;

um = 10^-6;

numOfWavlen = 79;

% 定义操作

Z  = @(x) zeropad(x,nullpixels);                                           % zero-padded sample
C  = @(x) imgcrop(x,nullpixels);                                           % 裁剪操作

%% 载入计算结果
% E中储存着从FDTD的monitor中导出的所有结果，为(0,3)张量
% 第一个维度放平面的电场
% 第二个维度索引xyz分量
% 第三个维度索引波长

%% 取出各个分量
Ex=reshape(E.E(:,1,numOfWavlen),length(E.x),length(E.y));
Ey=reshape(E.E(:,2,numOfWavlen),length(E.x),length(E.y));
Ez=reshape(E.E(:,3,numOfWavlen),length(E.x),length(E.y));

[XX,YY] = meshgrid(xx,yy);

% 需要旋转后才和结构方向对得上
Ex=imrotate(Ex,90);
Ey=imrotate(Ey,90);
Ez=imrotate(Ez,90);

xum=xx./um;
yum=yy./um;

%% 计算farfield pattern

ENear = Ey;

figure
subplot(1,2,1)
imagesc(xum,yum,abs(ENear))
xlabel('um');ylabel('um');
title('intensityNearField')
axis equal
colorbar
subplot(1,2,2)
imagesc(xum,yum,angle(ENear))
xlabel('um');ylabel('um');
title('phaseNearField')
axis equal
colorbar


% 拓展图片
ENearZeropad    = Z(ENear);
[N,M] = size(ENearZeropad);    % size of the wavefront

% 定义感兴趣的需要计算的角度

kx = pi/params.pxsize*(-1:2/(M-1):1);
ky = pi/params.pxsize*(-1:2/(N-1):1);
[kX,kY] = meshgrid(kx,ky);

thetax = asind(kx/k0);
thetay = asind(ky/k0);
[thetaX,thetaY] = meshgrid(thetax,thetay);
thetaOrig = thetax;

tic
% 计算远场
EFar = fftshift(fft2(fftshift(ENearZeropad)));

% EFar = ENear;
% for ii = 1:numel(ENear)
%     % 夫琅禾费衍射计算远场，给出latex公式
%     % U_f(u,v)&=\frac{\text{A exp}[\left.j\frac k{2f}\left(1-\frac df\right)(u^2+v^2)\right]}{j\lambda f}\times\int_{-\infty}^{\infty}t_A(\xi,\eta)\exp\left[-j\frac{2\pi}{\lambda f}(\xi u+\eta v)\right]d\xi d\eta
%     
%     FTIntegral = ENear(:)'*exp(-1i*(kX(ii).*XX(:)+kY(ii).*YY(:))/(params.wavlen*f_FT));
%     EFar(ii) = (1/(1i*params.wavlen*f_FT))^2.*FTIntegral;
% end
disp(['远场计算完成'])
toc

% 寻找角度截断index
% 计算每个元素与a的绝对差值
NALimitAngle=asind(NA);
diff = abs(thetaOrig - NALimitAngle);
% 找到最小差值对应的索引
[~, idx] = min(diff);

% 重新按照角度截断裁剪
C  = @(x) imgcrop(x,numel(thetaOrig)-idx);                                  % 裁剪操作

ENear = C(ENearZeropad);             % zero-padded sample


EFar(tand(thetaX).^2+tand(thetaY).^2>tand(NALimitAngle).^2) = 0;
EFar  = C(EFar);                     % zero-padded sample
thetax = C(thetax);
thetay = C(thetay);

%% 画图
figure
subplot(1,2,1)
plotField2D(thetax,thetay,(abs(EFar)))
xlabel('deg'),ylabel('deg')
anglemap = angle(EFar);
subplot(1,2,2)
plotField2D(thetax,thetay,anglemap)
xlabel('deg'),ylabel('deg')

% =========================================================================

%% Auxiliary functions

% =========================================================================

function u = imgcrop(x,cropsize)

% =========================================================================

% Crop the central part of the image.

% -------------------------------------------------------------------------

% Input:    - x        : Original image.

%           - cropsize : Cropping pixel number along each dimension.

% Output:   - u        : Cropped image.

% =========================================================================
if numel(x)==length(x)
    u = x(cropsize+1:end-cropsize);
else
    u = x(cropsize+1:end-cropsize,cropsize+1:end-cropsize);
end
end

function u = zeropad(x,padsize)

% =========================================================================

% Zero-pad the image.

% -------------------------------------------------------------------------

% Input:    - x        : Original image.

%           - padsize  : Padding pixel number along each dimension.

% Output:   - u        : Zero-padded image.

% =========================================================================
u = padarray(x,[padsize,padsize],0);
end

function plotField2D(xx,yy,field)

imagesc(xx,yy,field)

% 获取当前坐标轴
ax = gca;
% 获取当前图像的大小
img_extent = get(ax, 'XLim'); % 获取x轴范围
xmin = img_extent(1);
xmax = img_extent(2);
img_extent = get(ax, 'YLim'); % 获取y轴范围
ymin = img_extent(1);  
ymax = img_extent(2);
axis([xmin xmax ymin ymax]); % 设置坐标轴范围
colorbar
axis equal

end