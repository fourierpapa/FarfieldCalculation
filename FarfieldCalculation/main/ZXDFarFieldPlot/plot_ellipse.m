function plot_ellipse(x0,y0,Ex,Ey,theta_x,theta_y,S3)
%%  ˵��
%������һ�������ڣ�x0��y0)������Բ���䳤����ֱ�Ϊa,b,��Բ��Z��תtheta��
%This is a function to draw a ellipse with the center at (x0,y0) 
%%
% x0=1;
% y0=1;
% a=2;
% b=1;
% theta=pi/3;
%%
num_t=1e3;
t=linspace(0,2*pi,num_t);
% c=sqrt(Ex^2+Ey^2);
c=100000000;
x=x0+Ex/0.4/10^3/c*cos(t+theta_x);
y=y0+Ey/0.4/10^3/c*cos(t+theta_y);
global h1; %�����ô�����ɾ��

if S3>0
    h1=plot(x,y,'red','LineWidth',1.5);
else 
    S3<0
    h1=plot(x,y,'blue','LineWidth',1.5);
end
end