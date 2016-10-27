% function from deep irf to DC:
function  [U]=get_dt_shift(U,bigC,FT,rgs,dt_shift)

N=FT.nwin;omega=FT.omega;
for k=1:3
    Uh(k,:)=fft(U(k,:));
end
Uh(:,1)=0;
for ifreq=2:length(bigC)
    C=bigC(ifreq).a;
    Uh(1,ifreq) = Uh(1,ifreq).*exp(1i*omega(ifreq)*(rgs/C(1).cr(1)+dt_shift(1)));
    Uh(2,ifreq) = Uh(2,ifreq).*exp(1i*omega(ifreq)*(rgs/C(1).cl(1)+dt_shift(2)));
    Uh(3,ifreq) = Uh(3,ifreq).*exp(1i*omega(ifreq)*(rgs/C(1).cr(1)+dt_shift(3)));
end
for k=1:3
    toto(1:N)=Uh(k,1:N);
    U(k,:)=ifft(toto(:),'symmetric');
end
