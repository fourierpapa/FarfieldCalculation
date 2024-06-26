function b=Bluestein_dft(x,f1,f2,fs,mout)
[m,n]=size(x);
%------- Shifting
f11=f1+(mout*fs+f2-f1)/(2*mout);
f22=f2+(mout*fs+f2-f1)/(2*mout);

a = exp(1j*2*pi*f11/fs);
w = exp(-1j*2*pi*(f22-f11)/(mout*fs));
%------- Premultiply data
h=(-m+1):max(mout-1,m-1);
mp=m+mout-1;
h=w.^((h.^2)./2);
b=a.^(-(0:(m-1))).*h(m:(2*m-1));
%------- Fast convolution via FFT
b=fft(x.*repmat(b.',[1 n]),2^nextpow2(mp));
ft=fft(1./h(1:mp),2^nextpow2(mp));
b=ifft(b.*repmat(ft.',[1 n]));
%------- Final multiply.
b=(b(m:mp,1:n)).'.*repmat((h(m:mp)),[n 1]);
%------- Shifting
l=linspace(0,mout-1,mout);
l=l./mout.*(f22-f11)+f11;
Mshift=-m/2; 
Mshift=repmat(exp(-1i.*2*pi.*l.*(Mshift+1/2)/fs),[n 1]);
b=b.*Mshift;
end
