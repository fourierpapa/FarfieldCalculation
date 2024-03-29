function [img,cbarimg] = visualizeComplex(cimg, cmap)

amp = abs(cimg);
pha = angle(cimg);

amin = 0;
amax = max(amp(:));

ncmap = length(cmap);
img = zeros(size(cimg,1),size(cimg,2),3);

for i = 1:size(cimg,1)
    for j = 1:size(cimg,2)
        img(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(pha(i,j)+pi)),:) .* amp(i,j)./amax;
    end
end

n = 256;
cbarimg = nan(n,n,3);
x = linspace(-1,1,n);
y = x;
[X,Y] = meshgrid(x,y);
[theta,rho] = cart2pol(X,Y);
for i = 1:n
    for j = 1:n
        if rho(i,j) <= 1
            % 2024年1月3日
            % 这里theta前给了负，以对上自己左旋和右旋的理解
            % 但没有什么根据
            cbarimg(i,j,:) = cmap(1+round((ncmap-1)/2/pi*(-theta(i,j)+pi)),:) .* rho(i,j);
        end
    end
end

end