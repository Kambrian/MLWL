T=4;   %data window <---> leakage in frequency space (determines spetrum convolution kernal width)
            % also sets spectra physical resolution, if this is the data
            % period: w0=2*pi/T;
N=64;   %sampling frequency, determines extension/period/(overlap seperation) in spectra: wp=2*pi/dt=N*w0
            % changing this has no further effect on data quality when wp is big enough: wp>2*wmax
            % for fixed dt, N sets T thus w0: spectra resolution and leakage
dt=T/N;  

dk=2*pi/(2*N*dt);

x=dt*(0:N);
y=exp(-x.^2/2)*sqrt(2*pi);
% yy=y;N=N/2;
yy=[y,y(end-1:-1:2)];    %expand to symmetric form, so that its FT is real: y(N-n)=y(n)
yyk=ifft(yy/dk,'symmetric');           %restore from DFT amplitude to FT amplitude
z=imag(yyk)./real(yyk);
% figure;plot(z)
yyk(1)
%%
figure;
plot(x,y);
hold on;
xk=2*pi/((2*N)*dt)*(0:N);
yk=yyk(1:N+1);
plot(xk,abs(yk),'g');
plot(xk,real(yk),'r--');
xlim([0,T]);
% set(gca,'yscale','log')

yp=exp(-xk.^2/2);
figure;
plot(xk,abs(real(yk))./yp,'-',xk,abs(yk)./yp,'r--')
xlim([0,T])
% set(gca,'yscale','log')
%%
% if sampling from (-T/2,T/2), this corresponds to a phase shift:
% exp(i*w*T/2)=exp(i*k*w0*T/2)=exp(i*k*pi)=+1/-1

T=4;   %data window <---> leakage in frequency space (determines spetrum convolution kernal width)
            % also sets spectra physical resolution, if this is the data
            % period: w0=2*pi/T;
N=32;   %sampling frequency, determines extension/period/(overlap seperation) in spectra: wp=2*pi/dt=N*w0
            % changing this has no further effect on data quality when wp is big enough: wp>2*wmax
dt=T/N;  


x=dt*(-N:N-1);
y=exp(-x.^2/2);
yy=y;
% yy=[y,y(end-1:-1:2)];    %expand to symmetric form, so that its FT is real: y(N-n)=y(n)
yyk=fft(yy)*dt;           %restore from DFT amplitude to FT amplitude
yyk=yyk/sqrt(2*pi);
z=imag(yyk)./real(yyk);
figure;plot(z)

figure;
plot(x,y);
hold on;
xk=2*pi/((2*N)*dt)*(0:N);
yk=yyk(1:N+1);
plot(xk,abs(yk),'g');
plot(xk,abs(real(yk)),'r--');
xlim([0,T]);
% set(gca,'yscale','log')

yp=exp(-xk.^2/2);
figure;
plot(xk,real(yk)./yp,'-',xk,abs(yk)./yp,'r--')
xlim([0,T])
% set(gca,'yscale','log')
%% without symmetric expansion (lose even parity); 
% left half of the signal is truncated, leading to a severely smeared spectrum (convolved with exp(i*w*T/2)*sinc(wT/2))
% this must be avoided!!

T=5;   %data window <---> leakage in frequency space (determines spetrum convolution kernal width)
            % also sets spectra physical resolution, if this is the data
            % period: w0=2*pi/T;
N=320;   %sampling frequency, determines extension/period/(overlap seperation) in spectra: wp=2*pi/dt=N*w0
            % changing this has no further effect on data quality when wp is big enough: wp>2*wmax
dt=T/N;  


x=dt*(0:N);
y=exp(-x.^2/2);
yy=y;N=N/2;
% yy=[y,y(end-1:-1:2)];    %expand to symmetric form, so that its FT is real: y(N-n)=y(n)
yyk=fft(yy)*dt;           %restore from DFT amplitude to FT amplitude
yyk=yyk/sqrt(2*pi);
z=imag(yyk)./real(yyk);
figure;plot(z)

figure;
plot(x,y);
hold on;
xk=2*pi/((2*N)*dt)*(0:N);
yk=yyk(1:N+1);
plot(xk,abs(yk),'g');
plot(xk,real(yk),'r--');
xlim([0,T]);
% set(gca,'yscale','log')

yp=exp(-xk.^2/2);
figure;
plot(xk,abs(real(yk))./yp,'-',xk,abs(yk)./yp,'r--')
xlim([0,T])
% set(gca,'yscale','log')