% function from deep irf to DC:
function  [R,T,Z]=irf2dc_fund(M,FT,bigC,ZZ,ZR,RZ,RR,TT,HS,az,Sh,dr)
% remember the conventions, R = radial, T transverse, Z vertical downward.

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

imode=1;
for ifreq=2:length(bigC)
    if isempty(bigC(ifreq).a);continue;end
        C = bigC(ifreq).a;
        [~,ib]=unique(C(1).zz);
        ux0=C(1).ux(end);uz0=C(1).uz(end);
        uy0(1)=C(1).uy(end);
        uxi(imode) = interp1(C(1).zz(ib),C(1).ux(ib),HS,'linear');
        uyi(imode) = interp1(C(1).zz(ib),C(1).uy(ib),HS,'linear');
        uzi(imode) = interp1(C(1).zz(ib),C(1).uz(ib),HS,'linear');
        duxi(imode) = -interp1(C(1).zz(ib),C(1).Dux(ib),HS,'linear');
        duyi(imode) = -interp1(C(1).zz(ib),C(1).Duy(ib),HS,'linear');
        duzi(imode) = -interp1(C(1).zz(ib),C(1).Duz(ib),HS,'linear');
    
    for ir=1:length(ZZ)  
        if isempty(ZZ{ir});continue;end
        %% note
        % with matlab FFT convention, we had F(w) = exp(iwt-ikr)
        %% Love waves
        DT(imode)=dr{ir}/C(1).cl(imode);
        Th{ir}(ifreq) =  ( 1i*C(1).kl(imode)*M2{ir}(2,1)*uyi(imode)/uy0(imode) ...
            + duyi(imode)/uy0(imode)*M2{ir}(2,3)   )*TT{ir}(imode,ifreq).*exp(-1i*C(1).omega*DT(imode));
        %% Rayleigh waves
        DT(imode)=dr{ir}/C(1).cr(imode);
        Rh{ir}(ifreq) =   ((1i*C(1).kr(imode)*M2{ir}(3,1)*uzi(imode)/uz0(imode) + ...
            duzi(imode)/uz0(imode)*M2{ir}(3,3) )*RZ{ir}(imode,ifreq) ...
                + (1i*C(1).kr(imode)*M2{ir}(1,1)*uxi(imode)/ux0(imode) + ...
                duxi(imode)/ux0(imode)*M2{ir}(1,3) )*RR{ir}(imode,ifreq) ).*exp(-1i*C(1).omega*DT(imode) );

        Zh{ir}(ifreq) =  ( ( 1i*C(1).kr(imode)*M2{ir}(1,1)*uxi(imode)/ux0(imode) + ...
            duxi(imode)/ux0(imode)*M2{ir}(1,3) )*ZR{ir}(imode,ifreq)  ...
                + ( 1i*C(1).kr(imode)*M2{ir}(3,1)*uzi(imode)/uz0(imode) + ...
                duzi(imode)/uz0(imode)*M2{ir}(3,3) )*ZZ{ir}(imode,ifreq) ).*exp(-1i*C(1).omega*DT(imode) ) ;
                

        Th{ir}(ifreq) = Th{ir}(ifreq)*Sh(ifreq)*1E-3;
        Rh{ir}(ifreq) = Rh{ir}(ifreq)*Sh(ifreq)*1E-3;
        Zh{ir}(ifreq) = Zh{ir}(ifreq)*Sh(ifreq)*1E-3;
end
end
for ir=1:length(ZZ)
Th{ir}(ifreq:FT.nwin)=0;Th{ir}(1)=0;
T{ir}  = (ifft(Th{ir},'symmetric'));

Rh{ir}(ifreq:FT.nwin)=0;Rh{ir}(1)=0;
R{ir}  = (ifft(Rh{ir},'symmetric'));

Zh{ir}(ifreq:FT.nwin)=0;Zh{ir}(1)=0;
Z{ir}  = (ifft(Zh{ir},'symmetric'));
end