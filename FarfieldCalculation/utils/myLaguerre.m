function y=myLaguerre(p,l,x)
%%This file generates the Laguerre functions stired in the variable "y" The
%%inputs are: the indices "p","l" and the vector "x".
y=zeros(p+1,1);
if p==0
    y=1;
else
for m=0:p
    y(p+1-m)=((-1).^m.*(factorial(p+l)))./(factorial(p-m).*factorial(l+m).*factorial(m));
end
end
y=polyval(y,x);