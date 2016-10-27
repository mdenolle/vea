% function from deep irf to DC:
function  [Gdc]=get_dc_irf(M,bigC,G,FT,HS,az)
M = rot_nez2rtz(M,az);

fmax=FT.Fmax;
N=FT.nwin;
df=FT.df;
Gdc=zeros(3,N);
Gh=zeros(size(G));U=zeros(3,N);
for k=1:9
    Gh(k,:)=fft(G(k,:)*1E-14);  % this conversion comes from converting all velocities from km/s to m/s
end
w =zeros(N,1);
for ifreq=2:floor(fmax/df)
    C = bigC(ifreq).a;
    [~,ib]=unique(C(1).zz);
    uyi = interp1(C(1).zz(ib),C(1).uy(ib),HS,'spline');
    uxi = interp1(C(1).zz(ib),C(1).ux(ib),HS,'spline');
    uzi = interp1(C(1).zz(ib),C(1).uz(ib),HS,'spline');
    duyi = -interp1(C(1).zz(ib),C(1).Duy(ib),HS,'spline');
    duxi = -interp1(C(1).zz(ib),C(1).Dux(ib),HS,'spline');
    duzi = -interp1(C(1).zz(ib),C(1).Duz(ib),HS,'spline');
    ii=size(C(1).uy);
    if min(ii)==1
        uy0=C(1).uy(end);ux0=C(1).ux(end);uz0=C(1).uz(end);
    else
        uy0=C(1).uy(end,1);ux0=C(1).ux(end,1);uz0=C(1).uz(end,1);
    end
    plot(C(1).zz,C(1).ux)
%     plot(C(1).omega/(2*pi),uxi,'o');hold on;%pause
U(2,ifreq) = ( 1i*C(1).kl(1)*M(2,1)*uyi/uy0 + duyi/uy0*M(2,3)   )*Gh(5,ifreq);

U(1,ifreq) =  ((1i*C(1).kr(1)*M(3,1)*uzi/uz0 + duzi/uz0*M(3,3) )*Gh(3,ifreq) ...
        + ( 1i*C(1).kr(1)*M(1,1)*uxi/ux0 + duxi/ux0*M(1,3) )*Gh(1,ifreq) );
    
U(3,ifreq) = ( ( 1i*C(1).kr(1)*M(1,1)*uxi/ux0 + duxi/ux0*M(1,3) )*Gh(7,ifreq) ...
        + ( 1i*C(1).kr(1)*M(3,1)*uzi/uz0 + duzi/uz0*M(3,3) )*Gh(9,ifreq) ) ;


end
    
for k=1:3
    toto(1:N)=U(k,1:N);
    Gdc(k,:)=ifft(toto(:),'symmetric')*1E-3; % again converting because k = 1/km
end
