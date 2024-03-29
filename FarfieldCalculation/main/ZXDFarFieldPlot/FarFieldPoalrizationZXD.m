clc
close all
clear all

kl = 0;
kh = 0.07;
step = 0.002;

Data = load("Polarization_mode2_kx_ky_2_small range.txt");

kcorx = Data(:,1);
kcory = Data(:,2);
Len = length(kcorx);

kx = Data(:,3);
ky = Data(:,4);
kz = Data(:,5);

kpall = [kx ky]./sqrt(kx.^2+ky.^2);
kpall = [kpall zeros(Len,1)];

re_dx = Data(:,6);
im_dx = Data(:,7);
re_dy = Data(:,8);
im_dy = Data(:,9);
re_dz = Data(:,10);
im_dz = Data(:,11);

Re = [re_dx re_dy re_dz];
Im = [im_dx im_dy im_dz];

k = [kx ky kz];
kunit = k./sqrt(kx.^2+ky.^2+kz.^2);
x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

%sp basis vector
S = zeros(Len,3);
P = zeros(Len,3);
%sp projected to xy plane
Sproj = zeros(Len,3);
Pproj = zeros(Len,3);

for i = 1:Len
    s = cross(z,kunit(i,:))./(norm(cross(z,kunit(i,:))));
    S(i,:) = s;
    p = cross(kunit(i,:),s);
    P(i,:) = p;
    sproj = cross(z,kpall(i,:));
    Sproj(i,:) = sproj;
    pproj = cross(z,sproj);
    Pproj(i,:) = pproj;
end

Reds = zeros(Len,1);
Imds = zeros(Len,1);
Redp = zeros(Len,1);
Imdp = zeros(Len,1);
Re1proj = zeros(Len,3);
Im1proj = zeros(Len,3);
Re2proj = zeros(Len,3);
Im2proj = zeros(Len,3);
reDxproj = zeros(Len,1);
imDxproj = zeros(Len,1);
reDyproj = zeros(Len,1);
imDyproj = zeros(Len,1);

for j = 1:Len
    Reds(j,:) = dot(S(j,:),Re(j,:));
    Re1proj(j,:) = Reds(j)*Sproj(j,:);
    Imds(j,:) = dot(S(j,:),Im(j,:));
    Im1proj(j,:) = Imds(j)*Sproj(j,:);
    Redp(j,:) = dot(P(j,:),Re(j,:));
    Re2proj(j,:) = Redp(j)*Pproj(j,:);
    Imdp(j,:) = dot(P(j,:),Im(j,:));
    Im2proj(j,:) = Imdp(j)*Pproj(j,:);

    reDxproj(j,1) = Re1proj(j,1) + Re2proj(j,1);
    reDyproj(j,1) = Re1proj(j,2) + Re2proj(j,2);
    imDxproj(j,1) = Im1proj(j,1) + Im2proj(j,1);
    imDyproj(j,1) = Im1proj(j,2) + Im2proj(j,2);

    if isnan(reDxproj(j,1))
        reDxproj(j,1) = re_dx(j);
        reDyproj(j,1) = re_dy(j);
        imDxproj(j,1) = im_dx(j);
        imDyproj(j,1) = im_dy(j);
    end
end

Dx = reDxproj + 1i*imDxproj;
Dy = reDyproj + 1i*imDyproj;

%x-y basis momentum space radiative electric field intensity
Radiation_EF = reshape(abs(Dx).^2+abs(Dy).^2,sqrt(Len),[]);

%Stokes parameters
theta = atan2(imDxproj,reDxproj) - atan2(imDyproj,reDyproj);
S0 = abs(Dx).^2 + abs(Dy).^2;
S1 = (abs(Dx).^2 - abs(Dy).^2)./S0;
S2 = 2*abs(Dx).*abs(Dy).*cos(theta)./S0;
S3 = 2*abs(Dx).*abs(Dy).*sin(theta)./S0;

S1r = reshape(S1,sqrt(Len),[]);
S2r = reshape(S2,sqrt(Len),[]);
S3r = reshape(S3,sqrt(Len),[]);
S1r = S1r';
S2r = S2r';
S3r = S3r';
% 
% phi = reshape(atan(S1 + 1i.*S2),sqrt(Len),[]);
% 

kxx = kl:step:kh;
kyy = kl:step:kh;

figure
pcolor(kxx,kyy,atan2(S1r,S2r))
axis square
colormap hsv
colorbar
shading flat
shading interp 
% % xlim([kl kh])
% % ylim([kl kh])
% xlabel('\itk_xa/2\pi')
% ylabel('\itk_ya/2\pi')
% % set(gca,'XTick',[-0.02:0.01:0.02])
% % set(gca,'xticklabel',[])
% % set(gca,'YTick',[-0.02:0.01:0.02])
% % set(gca,'yticklabel',[])
% set(gca,'FontName','Arial','FontSize',14)
% % set(gca,'unit','centimeters','position',[3 2 8 8])
% title('S1/S2')

xlim([0 0.07])
ylim([0 0.07])
xlabel('\itk_xa/2\pi')
ylabel('\itk_ya/2\pi')
set(gca,'XTick',[-0.07:0.07:0.07])
% set(gca,'xticklabel',[])
set(gca,'YTick',[-0.07:0.07:0.07])
% set(gca,'yticklabel',[])
set(gca,'FontName','Arial','FontSize',14)
set(gca,'unit','centimeters','position',[3 2 8 8])


figure
pcolor(kxx,kyy,S3r)
axis square
colormap hsv
colorbar
shading flat
shading interp 
% xlim([kl kh])
% ylim([kl kh])
xlabel('\itk_xa/2\pi')
ylabel('\itk_ya/2\pi')
% set(gca,'XTick',[-0.02:0.01:0.02])
% set(gca,'xticklabel',[])
% set(gca,'YTick',[-0.02:0.01:0.02])
% set(gca,'yticklabel',[])
set(gca,'FontName','Arial','FontSize',14)
% set(gca,'unit','centimeters','position',[3 2 8 8])
title('S3')

figure
for k=1:1:Len 
    if i == 1000000
        scatter(kcorx(k),kcory(k),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'LineWidth',7)
    elseif i==100000
        scatter(kcorx(k),kcory(k),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'LineWidth',7)
    else
        plot_ellipse(kcorx(k),kcory(k),abs(Dx(k)),abs(Dy(k))^2,atan2(imDxproj(k),reDxproj(k)),atan2(imDyproj(k),reDyproj(k)),S3(i))
    end
    hold on
end
axis square
% xlim([kl kh])
% ylim([kl kh])
% xlabel('\itk_xa/2\pi')
% ylabel('\itk_ya/2\pi')
% % set(gca,'XTick',[-0.02:0.01:0.02])
% % set(gca,'xticklabel',[])
% % set(gca,'YTick',[-0.02:0.01:0.02])
% % set(gca,'yticklabel',[])
% set(gca,'FontName','Arial','FontSize',14)
% % set(gca,'unit','centimeters','position',[3 2 8 8])
% title('Ellipse')

xlim([-0.07 0.07])
ylim([-0.07 0.07])
xlabel('\itk_xa/2\pi')
ylabel('\itk_ya/2\pi')
set(gca,'XTick',[-0.07:0.07:0.07])
% set(gca,'xticklabel',[])
set(gca,'YTick',[-0.07:0.07:0.07])
% set(gca,'yticklabel',[])
set(gca,'FontName','Arial','FontSize',14)
set(gca,'unit','centimeters','position',[3 2 8 8])