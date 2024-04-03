% function fdtdFarfieldCalculation
% 2024年3月25日 由fyb创建
% 功能：根据从FDTD输出的电场，给出辐射远场图
% 特点：使用了FFT计算远场
% 
clear all;
close all;
clc;
dbstop if error

tic

% load functions
addpath(genpath('./utils'))
addpath(genpath('../src'))

%% 定义参数
% physical parameters
% 载入待计算的场
dataAddr = "./model/comsol/data/";
dataName = "1x8Bloch2Unit/";
% dataName = "Bloch";

ExRealorig = load(strcat(dataAddr,dataName,'ExReal.txt'));
EyRealorig = load(strcat(dataAddr,dataName,'EyReal.txt'));
EzRealorig = load(strcat(dataAddr,dataName,'EzReal.txt'));
ExImagorig = load(strcat(dataAddr,dataName,'ExImag.txt'));
EyImagorig = load(strcat(dataAddr,dataName,'EyImag.txt'));
EzImagorig = load(strcat(dataAddr,dataName,'EzImag.txt'));

% Exorig = load('./comsol/data/FullFieldEx.txt');
% Eyorig = load('./comsol/data/FullFieldEy.txt');
% Ezorig = load('./comsol/data/FullFieldEz.txt');

xraw = ExRealorig(:,1)*1e-9;
yraw = ExRealorig(:,2)*1e-9;
xx = unique(sort(xraw));
yy = unique(sort(yraw));

params.pxsize = xx(2)-xx(1);            % pixel size (m)
params.wavlen = 0.6e-6;                 % wavelength (m)
params.dist   = 5e-3;                      % imagicng distance (m)
k0 = 2*pi/params.wavlen;

kernelsize = 4;
params.nullpixels = ceil(kernelsize / params.pxsize);                             % number of padding pixels
params.nullpixels = 8000;

f_FT = 5e-3;                      % lens focal length (m)
NA = 0.8;

um = 10^-6;

numOfWavlen = 79;

%% 载入计算结果
% E中储存着从FDTD的monitor中导出的所有结果，为(0,3)张量
% 第一个维度放平面的电场
% 第二个维度索引xyz分量
% 第三个维度索引波长

%% 取出各个分量
Exraw = ExRealorig(:,4) + 1i*ExImagorig(:,4);
Eyraw = EyRealorig(:,4) + 1i*EyImagorig(:,4);
Ezraw = EzRealorig(:,4) + 1i*EzImagorig(:,4);

[XX,YY] = meshgrid(xx,yy);

Ex = griddata(xraw,yraw,Exraw,XX,YY);
Ey = griddata(xraw,yraw,Eyraw,XX,YY);
Ez = griddata(xraw,yraw,Ezraw,XX,YY);
Ex(isnan(Ex)) = 0;
Ey(isnan(Ey)) = 0;
Ez(isnan(Ez)) = 0;

xum=xx./um;
yum=yy./um;

ENear = Ex - 1i*Ey;

figure
subplot(3,2,1);imagesc(xum,yum,abs(Ex));xlabel('um');ylabel('um');title('intensityNearFieldEx');axis equal;colorbar
subplot(3,2,2);imagesc(xum,yum,angle(Ex));xlabel('um');ylabel('um');title('phaseNearFieldEx');axis equal;colorbar

subplot(3,2,3);imagesc(xum,yum,abs(Ey));xlabel('um');ylabel('um');title('intensityNearFieldEy');axis equal;colorbar
subplot(3,2,4);imagesc(xum,yum,angle(Ey));xlabel('um');ylabel('um');title('phaseNearFieldEy');axis equal;colorbar

subplot(3,2,5);imagesc(xum,yum,abs(Ez));xlabel('um');ylabel('um');title('intensityNearFieldEz');axis equal;colorbar
subplot(3,2,6);imagesc(xum,yum,angle(Ez));xlabel('um');ylabel('um');title('phaseNearFieldEz');axis equal;colorbar

%% 计算farfield pattern
% 定义操作

Z  = @(x) zeropad(x,params.nullpixels);                                           % zero-padded sample
% Q  = @(x) fftshift(fft2(fftshift(Z(x))));
Q  = @(x) gather(fftshift(fft2(fftshift(gpuArray(Z(x))))));

% 拓展图片
ENearZeropad    = Z(ENear);
[N,M] = size(ENearZeropad);    % size of the wavefront

% 定义感兴趣的需要计算的角度

kx = pi/params.pxsize*(-1:2/(M-1):1);
ky = pi/params.pxsize*(-1:2/(N-1):1);
[kX,kY] = meshgrid(kx,ky);

thetaxFarFullSize = asind(kx/k0);
thetayFarFullSize = asind(ky/k0);
[thetaX,thetaY] = meshgrid(thetaxFarFullSize,thetayFarFullSize);
thetaFarFullSize = thetaxFarFullSize;

disp(['开始远场计算'])
toc

% 计算远场
EFarFullSize  = Q(ENear);
ExFarFullSize = Q(Ex);
EyFarFullSize = Q(Ey);
EzFarFullSize = Q(Ez);

disp(['远场计算完成'])
toc

%% 寻找角度截断index
kx = 0.1;
p = 360e-9;
NA = kx*(2*pi/p)/k0;
NA = 0.99;
NA

% 计算每个元素与a的绝对差值
NALimitAngle = asind(NA)
diff = abs(thetaFarFullSize - NALimitAngle);
% 找到最小差值对应的索引
[~, idx] = min(diff);

% 重新按照角度截断裁剪
C  = @(x) imgcrop(x,numel(thetaFarFullSize)-idx);                           % 裁剪操作

EFarFullSize(tand(thetaX).^2+tand(thetaY).^2>tand(NALimitAngle).^2) = 0;
EFar  = C(EFarFullSize);                     % zero-padded sample
ExFar  = C(ExFarFullSize);                     % zero-padded sample
EyFar  = C(EyFarFullSize);                     % zero-padded sample
EzFar  = C(EzFarFullSize);                     % zero-padded sample
thetax = C(thetaxFarFullSize);
thetay = C(thetayFarFullSize);

%% 画图
figure
subplot(3,2,1);plotFarField2D(thetax,thetay,(abs(ExFar)));xlabel('deg'),ylabel('deg');
subplot(3,2,2);plotFarField2D(thetax,thetay,angle(ExFar));xlabel('deg'),ylabel('deg')

subplot(3,2,3);plotFarField2D(thetax,thetay,(abs(EyFar)));xlabel('deg'),ylabel('deg');
subplot(3,2,4);plotFarField2D(thetax,thetay,angle(EyFar));xlabel('deg'),ylabel('deg')

subplot(3,2,5);plotFarField2D(thetax,thetay,(abs(EzFar)));xlabel('deg'),ylabel('deg');
subplot(3,2,6);plotFarField2D(thetax,thetay,angle(EzFar));xlabel('deg'),ylabel('deg')

figure
subplot(1,2,1);plotFarField2D(thetax,thetay,(abs(EFar)));xlabel('deg'),ylabel('deg');
subplot(1,2,2);plotFarField2D(thetax,thetay,angle(EFar));xlabel('deg'),ylabel('deg')

%% 自定义画图

% 计算斯托克斯参量
ExS = trapz(yy,trapz(xx, Ex, 1))/trapz(yy,trapz(xx, ones(size(Ex)), 1));
EyS = trapz(yy,trapz(xx, Ey, 1))/trapz(yy,trapz(xx, ones(size(Ex)), 1));

ExS = ExFar;
EyS = EyFar;

SDelta = -(angle(ExS) - angle(EyS));
S0 = conj(ExS).*ExS + conj(EyS).*EyS;
S1 = (conj(ExS).*ExS - conj(EyS).*EyS)./S0;
S2 = 2*abs(ExS).*abs(EyS).*cos(SDelta)./S0;
S3 = 2*abs(ExS).*abs(EyS).*sin(SDelta)./S0;
S3GammaCalculate = S3(thetax==0,thetay==0)

EFarPlot = S0;
figure
title(EFarPlot)
subplot(1,2,1);plotFarField2D(thetax,thetay,(abs(EFarPlot)));xlabel('deg'), ylabel('deg');title('intensityFarField');axis equal;colorbar
subplot(1,2,2);plotFarField2D(thetax,thetay,angle(EFarPlot));xlabel('deg'), ylabel('deg');title('phaseFarField')    ;axis equal;colorbar

% 图如果太小,比如每个pix都要算，就需要把图插值放大
scalingRatio = 1;

cmap = sinebow(256); % 这个调色板好看
[img,cbar] = visualizeComplex(EFarPlot,cmap);
img = imresize(img,scalingRatio*size(img(:,:,1)),'bilinear');
figure,imshow(img)                                   % 画复振幅的幅值和相位

wavefront_jones = cat(3,ExS,EyS);

wavefront_jones = imresize(wavefront_jones,scalingRatio*size(wavefront_jones(:,:,1)),'bilinear');

visualizePolarization(wavefront_jones, cmap, 31);                          % 画偏振图

disp(['程序运行完成'])
toc

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

function plotFarField2D(xx,yy,field)
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


