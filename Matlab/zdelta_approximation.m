function a=zdelta_approximation(req,rp)

sigmz=100;
g=@(z) exp(-z.^2./2/sigmz^2)/sqrt(2*pi)/sigmz;

a=quad(@(z) gs(z),300,1e5)/(pi*req/(rp/req)*quad(g,300,1e5));
r=0:10:1000;
figure;
plot(r,pi*req/(rp/req)*g(r),'r');
hold on;
plot(r,gs(r),'g');
% plot(r,req^2./(rp^2+r.^2),'b');

function q=gs(z)
q=z;
for i=1:numel(z)
    q(i)=quadgk(@(zp) g(z(i)-zp)*req^2./(rp^2+zp.^2),-1e5,1e5);
end
end
end