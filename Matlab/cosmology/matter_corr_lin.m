function xi=matter_corr_lin(OmegaM,h,fb,sig8,r,z)
% matter correlation function from linear PT
% xi
% r: Mpc/h, comoving
%

tol=1e-2;  %relative tolerance
abstol=1e-7;  %absolute tol
inter=1000;

ns=1;  %primordial index

wx=@(x) 3*(sin(x)-x.*cos(x))./x.^3;
delt8=quadgk(@(k) wx(k*8).^2.*k.^(2+ns).*TF_BBKS(OmegaM,h,fb,k).^2,0,inf,'RelTol',1e-2);  %this integral is contributed primarily by k from 0.1 to 1
Ap=sig8.^2/delt8*growth_factor(OmegaM,z).^2; %normalization

xi=r;
for i=1:numel(r)
    if r(i)>1000
        int1=0;
    else
        %     disp(r(i));
        % fun=@(k) power_spectrum(OmegaM,h,fb,sig8,k,z).*k.*sin(k.*r(i))/r(i)/2/pi.^2;
        fun=@(k) Ap*TF_BBKS(OmegaM,h,fb,k).^2.*k.^(ns+1).*sin(k*r(i))/r(i);
        kmax=100;
        maxinter=max(inter,floor(kmax*r(i)));
        int1=quadgk(fun,0,kmax,'RelTol',tol,'AbsTol',abstol,'MaxIntervalCount',maxinter);
        %        int1=quadgk(fun,0,100,'RelTol',tol,'MaxIntervalCount',inter);
        
        for j=1:10
            int2=quadgk(fun,kmax,2*kmax,'RelTol',tol,'AbsTol',abstol,'MaxIntervalCount',maxinter);
            %         int2=quadgk(fun,kmax,2*kmax,'RelTol',tol,'MaxIntervalCount',maxinter);
            int1=int1+int2;
            kmax=kmax*2;
            if abs(int2)<tol*abs(int1)||abs(int2)<abstol
                %             disp('succeed');
                break;
            end
        end
        if j==10
            warning('MATLAB:integration','\n=====\nmaximum iteration exceeded towards inf, accuracy %g\n ======\n',abs(int2/int1));
        end
    end
    xi(i)=int1;
    
    % xi(i)=quadgk(fun,0,inf,'RelTol',tol,'AbsTol',abstol,'MaxIntervalCount',maxinter);
end