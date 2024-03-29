function [fig] = visualizePolarization(jvec, cmap, n_grid)
% 这个函数使用琼斯矢量画偏振分布
% 潜在的问题是左右旋都是圆的
% 但也可以通过看相位解决
ex = jvec(:,:,1);
ey = jvec(:,:,2);
amp = sqrt(abs(ex).^2 + abs(ey).^2);

amin = 0;
amax = max(amp(:));

% fig = figure;
% imshow(amp,[amin,amax],'border','tight');colormap(cmap);

nx = size(jvec,1); ny = size(jvec,2);
xc = linspace(0,nx,n_grid+2);
yc = linspace(0,ny,n_grid+2);
px = nx/n_grid; py = ny/n_grid;
ratio = 0.5;
scale = min([px,py]) * ratio / 2 / max([abs(ex(:));abs(ey(:))]);
for i = 1:n_grid
    for j = 1:n_grid
        ix = round(xc(i+1));
        iy = round(yc(j+1));
        th = linspace(0,2*pi,100);
        xe = ix + scale*abs(ex(iy,ix)).*cos(th+angle(ex(iy,ix)));
        ye = iy + scale*abs(ey(iy,ix)).*cos(th+angle(ey(iy,ix)));
        hold on,plot(xe,ye,'color',[1.0,0.3,0.3],'linewidth',1.5)
%         drawnow;
    end
end

end