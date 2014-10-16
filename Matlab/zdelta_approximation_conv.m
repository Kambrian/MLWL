rp=2;
sigmz=300;
f=@(w) exp(-rp*w);
g=@(w) exp(-w.^2*sigmz^2/2);
figure;
w=0.001:0.001:0.1;
plot(w,g(w),'r');
hold on;
plot(w,f(w).*g(w),'g')

req=10;
p=@(z) req^2./(rp^2+z.^2);
figure;
z=0.1:0.1:10;
plot(z,p(z))