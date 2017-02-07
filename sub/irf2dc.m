% function from deep irf to DC:
function  [R,T,Z]=irf2dc(M,FT,bigC,ZZ,ZR,RZ,RR,TT,HS,az,Sh,dr,typ)
% function  [R,T,Z,corrR,corrT,corrZ]=irf2dc(M,FT,bigC,ZZ,ZR,RZ,RR,TT,HS,az,Sh,dr,typ)
% remember the conventions, R = radial, T transverse, Z vertical downward.
% if typ='f', provide results in frequency domain
% if typ = 't', provide results in time domain
% rotate the moment tensor at each point source.
for ir=1:length(ZZ)
    cosa=cos(az{ir}*pi/180);cosb=cosa;
    sina=sin(az{ir}*pi/180);sinb=sina;
    R1 = [cosa sina 0;-sina cosa 0; 0 0 1];
    R2 = [cosb sinb 0;-sinb cosb 0; 0 0 1];
%     M2{ir}(1,1) = cosa*cosb*M(1,1)+cosb*sina*M(1,2)+sinb*cosa*M(2,1)+sina*sinb*M(2,2);
%     M2{ir}(1,2)= -cosb*sina*M(1,1)+cosb*cosa*M(1,2)-sina*sinb*M(2,1)+sinb*cosa*M(2,2);
%     M2{ir}(2,1) = -cosa*sinb*M(1,1)-sina*sinb*M(1,2)+cosa*cosb*M(2,1)+cosb*sina*M(2,2);
%     M2{ir}(2,2) = sina*sinb*M(1,1)-sinb*cosa*M(1,2)-cosb*sina*M(2,1)+cosa*cosb*M(2,2);
%     M2{ir}(2,3) = -sinb*M(1,3)+cosb*M(2,3);
%     M2{ir}(1,3) = cosb*M(1,3)+sinb*M(2,3);
%     M2{ir}(3,1) =cosa*M(1,3)+sina*M(2,3) ;
%     M2{ir}(3,2) = -sina*M(1,3)+cosa*M(2,3);
%     M2{ir}(3,3) = M(3,3);
    M2{ir} = R1*M*R2';
end
for ifreq=2:length(bigC)
    if isempty(bigC(ifreq).a);continue;end
    C = bigC(ifreq).a;
    kr(ifreq)=C(1).kr(1);
    kl(ifreq)=C(1).kl(1);
    [~,ib]=unique(C(1).zz);
    omega(ifreq)=C(1).omega;
%     Nmode = min(min(size(C(1).Dux)),min(min(size(C(1).ux))));
    
%     if Nmode==1;
        ux0(ifreq)=C(1).ux(end);%uz0(ifreq)=C(1).uz(end);uy0(ifreq)=C(1).uy(end);
        uxi(ifreq) = interp1(C(1).zz(ib),C(1).ux(ib),HS,'linear');
        uzi(ifreq) = interp1(C(1).zz(ib),C(1).uz(ib),HS,'linear');
        duxi(ifreq) = -interp1(C(1).zz(ib),C(1).Dux(ib),HS,'linear');
        duzi(ifreq) = -interp1(C(1).zz(ib),C(1).Duz(ib),HS,'linear');
        uyi(ifreq) = interp1(C(1).zz(ib),C(1).uy(ib),HS,'linear');
        duyi(ifreq) = -interp1(C(1).zz(ib),C(1).Duy(ib),HS,'linear');
end
uy0=1;
uz0=1;
for ir=1:length(ZZ)  
    if isempty(ZZ{ir});continue;end
    
        Th{ir}(2:ifreq) =0;
        Rh{ir}(2:ifreq) =0;
        Zh{ir}(2:ifreq) =0;
        Th{ir}(2:ifreq) =  ( 1i*kl(2:ifreq).*M2{ir}(2,1).*uyi(2:ifreq)./uy0 ...
                + duyi(2:ifreq)./uy0*M2{ir}(2,3) ).*TT{ir}(2:ifreq).*exp(-1i.*dr{ir}.*kl(2:ifreq));
        Rh{ir}(2:ifreq) = ((1i*kr(2:ifreq)*M2{ir}(3,1).*uzi(2:ifreq)./uz0 + ...
                duzi(2:ifreq)./uz0*M2{ir}(3,3) ).*RZ{ir}(2:ifreq) ...
                + ( 1i*kr(2:ifreq)*M2{ir}(1,1).*uxi(2:ifreq)./ux0(2:ifreq) + ...
                duxi(2:ifreq)./ux0(2:ifreq)*M2{ir}(1,3) ).*RR{ir}(2:ifreq) ).*exp(-1i.*dr{ir}.*kr(2:ifreq));
        Zh{ir}(2:ifreq) = (( 1i*kr(2:ifreq).*M2{ir}(1,1).*uxi(2:ifreq)./ux0(2:ifreq) + ...
                duxi(2:ifreq)./ux0(2:ifreq)*M2{ir}(1,3) ).*ZR{ir}(2:ifreq)  ...
                + ( 1i*kr(2:ifreq)*M2{ir}(3,1).*uzi(2:ifreq)/uz0 + ...
                duzi(2:ifreq)/uz0*M2{ir}(3,3) ).*ZZ{ir}(2:ifreq) ).*exp(-1i.*dr{ir}.*kr(2:ifreq)) ;

end
if strcmp(typ,'t')==1
for ir=1:length(ZZ)
Th{ir}(1:ifreq-1) = Th{ir}(1:ifreq-1).*Sh(1:ifreq-1)*1E-3;
Rh{ir}(1:ifreq-1) = Rh{ir}(1:ifreq-1).*Sh(1:ifreq-1)*1E-3;
Zh{ir}(1:ifreq-1) = Zh{ir}(1:ifreq-1).*Sh(1:ifreq-1)*1E-3;
Th{ir}(ifreq:FT.nwin)=0;Th{ir}(1)=0;
T{ir}  = (ifft(Th{ir},'symmetric'));

Rh{ir}(ifreq:FT.nwin)=0;Rh{ir}(1)=0;
R{ir}  = (ifft(Rh{ir},'symmetric'));

Zh{ir}(ifreq:FT.nwin)=0;Zh{ir}(1)=0;
Z{ir}  = (ifft(Zh{ir},'symmetric'));
end

else
    T=Th;R=Rh;Z=Zh;
end
