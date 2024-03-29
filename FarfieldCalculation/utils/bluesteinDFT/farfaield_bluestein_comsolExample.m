clc
close all
clear all
format long
dbstop if error

tic
%% Parameters
lamda = 1550e-9;   % 波长
k=2*pi/lamda;   % 波矢

p = 10e-6;     % 如果是全场仿真就放完整尺寸，是周期单元仿真就放单个元胞的尺寸
reptime = [1,1];   % 周期性结构的重复次数，如果是全场仿真就放[1,1]，是周期单元仿真就在xy方向拓展

NA = 0.42*4.73/3.66; % 物镜的数值孔径
anglemax = asind(NA);  % 远场的角度
mout = 2048;    % 分辨率

tmp=readmatrix('Ex.txt');
xraw = tmp(:,1)*1e-9;
yraw = tmp(:,2)*1e-9;
Exraw = tmp(:,4);

tmp=readmatrix('Ey.txt');
Eyraw = tmp(:,4);

tmp=readmatrix('Ez.txt');
xrawz = tmp(:,1)*1e-9;
yrawz = tmp(:,2)*1e-9;
Ezraw = tmp(:,4);

x = linspace(-p/2,p/2,1000)';
[xx,yy] = meshgrid(x,x);
Exx = griddata(xraw,yraw,Exraw,xx,yy);
Eyy = griddata(xraw,yraw,Eyraw,xx,yy);
Ezz = griddata(xrawz,yrawz,Ezraw,xx,yy);
Exx(isnan(Exx)==1) = 0;
Eyy(isnan(Eyy)==1) = 0;
Ezz(isnan(Ezz)==1) = 0;

Exxfull = [Exx]; %单个场分布时，Exx就是全部的Exxfull
Eyyfull = [Eyy];
Ezzfull = [Ezz];

figure();
subplot(231);imagesc(real(Exxfull));axis equal;axis tight;title('Intensity of Ex in the nearfield');colorbar()
subplot(232);imagesc(real(Eyyfull));axis equal;axis tight;title('Intensity of Ey in the nearfield');colorbar()
subplot(233);imagesc(real(Ezzfull));axis equal;axis tight;title('Intensity of Ez in the nearfield');colorbar()
subplot(234);imagesc(angle(Exxfull));axis equal;axis tight;title('Phase of Ex in the nearfield');colorbar()
subplot(235);imagesc(angle(Eyyfull));axis equal;axis tight;title('Phase of Ey in the nearfield');colorbar()
subplot(236);imagesc(angle(Ezzfull));axis equal;axis tight;title('Phase of Ez in the nearfield');colorbar()
frame=getframe(gcf);
imwrite(frame.cdata,'Eraw.png');
%% Bluestein dft  %只是一种提高fft分辨率的方法
pixel = p(1)/length(Exx);
fs = 1/pixel;
freint = cosd(90-anglemax)/lamda;   % Spatial frequency of interests
freqregion = [-freint,freint;-freint,freint];

% 这个有空找曲歌扬讨论一下
free = linspace(-freint,freint,mout);
alpha = 90-acosd(free*lamda);

% fft变换的方法
fy1=freqregion(1,1)+fs/2;
fy2=freqregion(1,2)+fs/2;
fx1=freqregion(2,1)+fs/2;
fx2=freqregion(2,2)+fs/2;

farx_blue = Bluestein_dft(Exxfull,fy1,fy2,fs,mout);
farx_blue = Bluestein_dft(farx_blue,fx1,fx2,fs,mout);
fary_blue = Bluestein_dft(Eyyfull,fy1,fy2,fs,mout);
fary_blue = Bluestein_dft(fary_blue,fx1,fx2,fs,mout);
farz_blue = Bluestein_dft(Ezzfull,fy1,fy2,fs,mout);
farz_blue = Bluestein_dft(farz_blue,fx1,fx2,fs,mout);
Efarsum = abs(farx_blue).^2+abs(fary_blue).^2;

figure();
subplot(221);imagesc(abs(farx_blue).^2);axis equal;axis tight;
title('Intensity of Ex in the farfield');colorbar();
subplot(222);imagesc(abs(fary_blue).^2);axis equal;axis tight;
title('Intensity of Ey in the farfield');colorbar()
subplot(223);imagesc(abs(farz_blue).^2);axis equal;axis tight;
title('Intensity of Ez in the farfield');colorbar()
subplot(224);imagesc(Efarsum);axis equal;axis tight;
title('Intensity of E in the farfield');colorbar()
frame=getframe(gcf);
imwrite(frame.cdata,'Efar.png');


figure();surf(alpha,alpha,Efarsum);shading interp;axis tight;title('Total intensity in the farfield');colorbar();
xlabel('Angle (Degree)');ylabel('Angle (Degree)')

figure();
imagesc(Efarsum);axis equal;axis tight;axis off;
set(gca,'Ydir','normal');
set(gcf,'Position',[600 300 512 512]);%消除白边
set(gca,'Position',[0 0 1 1]);%消除白边
title('Intensity in the farfield with PL');colormap hot;
frame=getframe(gcf);
imwrite(frame.cdata,['Efarsum.png']);
close


toc
%% 求添加偏振片后的场图 用ExEy的数据去做坐标变换即可，jones矩阵
%线偏振片
%jones_pol= [cosd(theta)^2 1/2*sind(2*theta);
% 1/2*sind(2*theta) sind(theta)^2]
% farx_blue = ones(2048,2048);
% fary_blue = ones(2048,2048);


for i = 1:6
    theta = 30*i;
    farx_pol=cosd(theta)^2*farx_blue+1/2*sind(2*theta)*fary_blue;
    fary_pol=1/2*sind(2*theta)*farx_blue+sind(theta)^2*fary_blue;

    Efarsum_pol(:,:,i) = abs(farx_pol).^2+abs(fary_pol).^2;
    figure();
    imagesc(Efarsum_pol(:,:,i));axis equal;axis tight;axis off;
    set(gca,'Ydir','normal');
    set(gcf,'Position',[600 300 512 512]);%消除白边
    set(gca,'Position',[0 0 1 1]);%消除白边
    title('Intensity in the farfield with PL');colormap hot;

    frame=getframe(gcf);
    imwrite(frame.cdata,[num2str(theta),'pol.png']);
    close
end

toc