% function from deep irf to DC:
function  [U]=get_dc_conv(U,FT,fc)

N=FT.nwin;
omega=FT.omega;
tau = 1/4*1/(2*fc);% pulse width prop to corner frequency
Shat=zeros(N,1);
Shat(2:N/2) = 4.*exp(-1i*omega(2:N/2)*2*tau).*(sin(omega(2:N/2)*tau/2).^2)...
    .*sin(omega(2:N/2)*tau)./((omega(2:N/2)*tau).^3);
a=size(U);
for k=1:min(a)
    toto(1:N)=fft(U(k,:));
    U(k,:)=ifft(toto(:).*Shat(:),'symmetric');
end
