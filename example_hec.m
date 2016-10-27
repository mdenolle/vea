%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        VIRTUAL EARTHQUAKES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% this script constucts the fundamental mode of surface wavesmakes virtual 
% earthquakes from the ambient noise green
% tensor. This is only for fundamental mode surface waves, please make sure
% that the green tensor obtained have one, dispersed wave packet. This
% method approxlimates that the ZT-TZ-RT-TR components are close to zero,
% i.e. that the Love and Rayleigh are not coupled.
% This code is based on :
%Denolle et al (2013), "Ground Motion Prediction of Realistic Earthquake Sources
%Using the Ambient Seismic Field" , J. Geophys. Res., Vol. 118, No. 5, pp. 
% 2102-2118, doi: 10.1029/2012JB009603
% and on:
% Denolle et al (2012), "Solving the Surface-Wave Eigenproblem with Chebyshev 
% Spectral Collocation", 2012, Bull. Seismol. Soc. Am., Vol. 102, No. 3, 
%pp. 1214-1223, doi: 10.1785/0120110183
% to solve the surface wave eigenproblem for an approxlimated 1D medium
% around the virtual source.

% note on manipulating Green's function. 
% 1) In this example, I average both sides of the cross correlation, but it
% remains your choice.
% 2) The Green's function is proportional to the time derivative of the 
% cross correlation. In this example, the cross-correlations are raw and 
% I apply a time derivative in frequency domain.

% the data is organized as follow: in ./tst/SOURCE/RECEIVER

% Written by Marine Denolle (10/25/16) be aware of bugs.
% Email me at mdenolle@fas.harvard.edu


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!rm -f *~              
clear all;close all;clc % clear all variables
addpath '../sub'         % add the path of all subroutines

%% Strings and directories
diran = './tst';           % directory of the data (greens' functions)
ff = '2013_001_00_00_00_2013_366_00_00_00';
info = importdata('./list_sta.dat');% lat-long file of receiver stations (net sta lat long el)
net=info.textdata(:,1);
chan=info.textdata(:,3);
sta=info.textdata(:,2);
lat=info.data(:,1);
lon=info.data(:,2);


%% source
snet={'CI'};              % station name of virtual source
ssta={'HEC'};              % station name of virtual source
SourceLat=34.813;               % real source latitude HEC
SourceLon=-116.419;             % real source longitude HEC

k = find(strcmp(ssta,sta)==1);
StaSourceLat=info.data(k,1);    % virtual source latitude
StaSourceLon=info.data(k,2);    % virtual source longitude
HS=5.5;                          % real source depth (km)
% from SCSN catalog (NED already)
Mo=[ -2.749e+23 -3.734e+23 -0.959e+23 ; -3.734e+23 3.052e+23 -0.802e+23;-0.959e+23  -0.802e+23 -0.304e+23 ]*1E-7; % HEC 

beta=3000 ; % average shear velocity m/s
Ds=3E6 ; % average stress drop 3MPa
fc = 0.491 * beta*(Ds/norm(Mo,2))^(1/3); % corner frequency

%% earthquake data:
dirdata = ['./data/Event_2008_12_06_04_18_42.85_Mw_5.06']; % hec

%% Time and Frequency
dt = 0.05;      % sampling rate in seconds
T1 = 4;         % shortest period to evaluate
T2 = 10;        % longer period to evaluate
Fnyq = 1/(2*dt);% Nyquist frequency (Hz)
[a,b]=butter(2,[1/T2 1/T1]*2*dt);
twin=500;

    
%% LOAD surface wave excitation
% in this section, you are supposed to have already solved the SW
% eigenproblem and just load the bigC variable:
% load HEC basin-type 1D velocity model
load('./mat/air_arbitrary_HEC.mat');

% source time function (moment rate function integrated to moment
% spectrum) offset to half duration.
S = 1./(1 + (FT.omega./2/pi/fc).^2).*exp(-1i*FT.omega/2/pi/fc)./FT.omega./1i;

%% Build Virtual Earthquake
for ista=1:length(net)
      list{ista} = [char(net(ista)) '.' char(sta(ista))];
        % distance between virtual source and receiver
       [rg(ista),az]=distance(StaSourceLat,StaSourceLon,lat(ista),lon(ista),almanac('earth','ellipsoid'));
       % distance between real source and receiver
       [b1,az1]=distance(SourceLat,SourceLon,lat(ista),lon(ista),almanac('earth','ellipsoid'));
       AZ{ista} = az1;
       dr{ista}=b1-rg(ista);
     %% raw Impulse Response Function: 
        direc = [diran '/' char(snet) '.' char(ssta)];
        
        disp(list{ista})
        %% only keep twin length
%         
        f1 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'E_' char(chan(ista)) 'E.tst'];
        f2 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'E_' char(chan(ista)) 'N.tst'];
        f3 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'E_' char(chan(ista)) 'Z.tst'];
        f4 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'N_' char(chan(ista)) 'E.tst'];
        f5 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'N_' char(chan(ista)) 'N.tst'];
        f6 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'N_' char(chan(ista)) 'Z.tst'];
        f7 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'Z_' char(chan(ista)) 'E.tst'];
        f8 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'Z_' char(chan(ista)) 'N.tst'];
        f9 = [direc '/' char(list{ista}) '/' ff '_' char(chan(k)) 'Z_' char(chan(ista)) 'Z.tst'];
        
        if exist(f1,'file')~=2||exist(f2,'file')~=2||exist(f3,'file')~=2||...
           exist(f4,'file')~=2||exist(f5,'file')~=2||exist(f6,'file')~=2||...
           exist(f7,'file')~=2||exist(f8,'file')~=2||exist(f9,'file')~=2;
      
    
        continue;end;
    
    
    %% read data
    ff1 = [dirdata '/' char(net(ista)) '.' char(sta(ista)) '.' char(chan(ista)) 'N.sac'];
    ff2 = [dirdata '/' char(net(ista)) '.' char(sta(ista)) '.' char(chan(ista)) 'E.sac'];
    ff3 = [dirdata '/' char(net(ista)) '.' char(sta(ista)) '.' char(chan(ista)) 'Z.sac'];
    if exist(ff1,'file')~=2||exist(ff2,'file')~=2||exist(ff3,'file')~=2;continue;end
    
        ee=readsac(f1);
        en=readsac(f2);
        ez=readsac(f3);
        ne=readsac(f4);
        nn=readsac(f5);
        nz=readsac(f6);
        ze=readsac(f7);
        zn=readsac(f8);
        zz=readsac(f9);
        
        n=readsac(ff1);
        e=readsac(ff2);
        z=readsac(ff3);
        
        L=floor(ee.npts/2);
        dt=nn.delta;N= (1:floor(twin/dt));
        NN = nn.trace;NN = NN(L+1:end-1) + flipdim(NN(1:L),1);
        NE = ne.trace;NE = NE(L+1:end-1) + flipdim(NE(1:L),1);
        NZ = nz.trace;NZ = NZ(L+1:end-1) + flipdim(NZ(1:L),1);
        
        EN = en.trace;EN = EN(L+1:end-1) + flipdim(EN(1:L),1);
        EE = ee.trace;EE = EE(L+1:end-1) + flipdim(EE(1:L),1);
        EZ = ez.trace;EZ = EZ(L+1:end-1) + flipdim(EZ(1:L),1);
        
        ZN = zn.trace;ZN = ZN(L+1:end-1) + flipdim(ZN(1:L),1);
        ZE = ze.trace;ZE = ZE(L+1:end-1) + flipdim(ZE(1:L),1);
        ZZ = zz.trace;ZZ = ZZ(L+1:end-1) + flipdim(ZZ(1:L),1);
        
        % azimuth between source and receiver
        phi = az1*pi/180;
        cosa=cos(phi);
        sina=sin(phi);
        cosb=cos(phi);
        sinb=sin(phi);
        %  Rotate green tensor from North-East-Down to
        %  Radial-Transverse-Down
        II=1:length(NN);
        G{ista}(1,II) = cosa*cosb*NN+cosa*sinb*NE+sina*cosb*EN+sina*sinb*EE;
        G{ista}(2,II) = -cosa*sinb*NN+cosa*cosb*NE-sina*sinb*EN+sina*cosb*EE;
        G{ista}(4,II) = -sina*cosb*NN-sina*sinb*NE+cosa*cosb*EN+cosa*sinb*EE;
        G{ista}(5,II) = sina*sinb*NN-sina*cosb*NE-cosa*sinb*EN+cosa*cosb*EE;
        G{ista}(3,II) = cosa*(NZ)+sina*(EZ);
        G{ista}(6,II) = -sina*(NZ)+cosa*(EZ);
        G{ista}(7,II) = cosb*(ZN)+sinb*(ZE);
        G{ista}(8,II) = -sinb*(ZN)+cosb*(ZE);
        G{ista}(9,II) = ZZ;
        
        % Fourier transform ZZ, ZR, RZ, RR, TT and take the time derivative
        % of cross correlation
        ZZh{ista}=fft(G{ista}(9,N));ZZh{ista}=1i*FT.omega.*ZZh{ista}(1:length(FT.omega));
        ZRh{ista}=fft(G{ista}(7,N));ZRh{ista}=1i*FT.omega.*ZRh{ista}(1:length(FT.omega));
        RRh{ista}=fft(G{ista}(1,N));RRh{ista}=1i*FT.omega.*RRh{ista}(1:length(FT.omega));
        RZh{ista}=fft(G{ista}(3,N));RZh{ista}=1i*FT.omega.*RZh{ista}(1:length(FT.omega));
        TTh{ista}=fft(G{ista}(5,N));TTh{ista}=1i*FT.omega.*TTh{ista}(1:length(FT.omega));
        
        % simply take the (low accuracy) time derivative of cross correlationd and filter
        % them just for comparison with the data waveforms.
        for ii=1:9
            G{ista}(ii,N(1:end-1))=diff(filtfilt(a,b,G{ista}(ii,N))); 
            G{ista}(ii,N(end))=0;
        end
        % Rotate earthquake data from North-East-Down to
        % Radial-Transvere-Down and downsample to Xcorr sampling rate to
        % 20sps
        Zdata{ista} = -filtfilt(a,b,resample(z.trace,20,floor(1/z.delta))); % data was originally up
        Rdata{ista} = cosb*filtfilt(a,b,resample(n.trace,20,floor(1/z.delta)))+sinb*filtfilt(a,b,resample(e.trace,20,floor(1/z.delta)));
        Tdata{ista} = -sinb*filtfilt(a,b,resample(n.trace,20,floor(1/z.delta)))+cosb*filtfilt(a,b,resample(e.trace,20,floor(1/z.delta)));
        Zdata{ista} = Zdata{ista}(floor(60/dt):end); % data starts 60s after earthquake
        Rdata{ista} = Rdata{ista}(floor(60/dt):end); % data starts 60s after earthquake
        Tdata{ista} = Tdata{ista}(floor(60/dt):end); % data starts 60s after earthquake
        max1 = max([max(abs( G{ista}(:,II)))])/5;
        max2 = max([max(abs(Zdata{ista})) max(abs(Rdata{ista})) max(abs(Tdata{ista})) ])/5;
        
          
end
    
% convert from surface impulse force to buried double couple
[R,T,Z]=irf2dc(Mo,FT,bigC,ZZh,ZRh,RZh,RRh,TTh,HS,AZ,S,dr);
for ir=1:length(R)
    if max(abs(R{ir})) ==0 ; continue;end
    % find the times of high seismic energy (signal)
    A = cumsum(abs(Zdata{ir}(N)).^2);
    IIZ=find(A>0.01*max(A)&A<0.95*max(A));
    A = cumsum(abs(Tdata{ir}(N)).^2);
    IIT=find(A>0.01*max(A)&A<0.95*max(A));
    A = cumsum(abs(Rdata{ir}(N)).^2);
    IIR=find(A>0.01*max(A)&A<0.95*max(A));
    
    % Filter predictions
    R{ir}=filtfilt(a,b,R{ir});
    T{ir}=filtfilt(a,b,T{ir});
    Z{ir}=filtfilt(a,b,Z{ir});
    
    % amplitude levels
    maxR(ir)=max(abs(Rdata{ir}));
    maxT(ir)=max(abs(Tdata{ir}));
    maxZ(ir)=max(abs(Zdata{ir}));
    facR(ir) = max(abs(Rdata{ir}))/max(abs(R{ir}));
    facT(ir) = max(abs(Tdata{ir}))/max(abs(T{ir}));
    facZ(ir) = max(abs(Zdata{ir}))/max(abs(Z{ir}));
    
    % correlation coefficients
    CCR(ir) =  sum(Rdata{ir}(IIR).*R{ir}(IIR)') / norm(Rdata{ir}(IIR))/norm(R{ir}(IIR));
    CCZ(ir) =  sum(Zdata{ir}(IIZ).*Z{ir}(IIZ)') / norm(Zdata{ir}(IIZ))/norm(Z{ir}(IIZ));
    CCT(ir) =  sum(Tdata{ir}(IIT).*T{ir}(IIT)') / norm(Tdata{ir}(IIT))/norm(T{ir}(IIT));
    
    CCRR(ir) =  sum(Rdata{ir}(IIR).* G{ir}(1,IIR)') / norm(Rdata{ir}(IIR))/norm(G{ir}(1,IIR));
    CCZZ(ir) =  sum(Zdata{ir}(IIZ).* G{ir}(9,IIZ)') / norm(Zdata{ir}(IIZ))/norm(G{ir}(9,IIZ));
    CCTT(ir) =  sum(Tdata{ir}(IIT).* G{ir}(5,IIT)') / norm(Tdata{ir}(IIT))/norm(G{ir}(5,IIT));
    
    
end
ik=find(facR~=0&isnan(facZ)==0);
figure(2)
subplot(311)
hist(log10(facR(ik)));set(gca,'Fontsize',14);title('factor (ratio of obs/pred) R');grid on;%xlim([-3 0])
subplot(312)
hist(log10(facT(ik)));set(gca,'Fontsize',14);title('fac (ratio of obs/pred) T');grid on;%xlim([-3 0])
subplot(313)
hist(log10(facZ(ik)));set(gca,'Fontsize',14);title('fac (ratio of obs/pred) Z');grid on;%xlim([-3 0])


figure(22)
subplot(311)
plot(cat(1,AZ{ik}),log10(facR(ik)),'o');set(gca,'Fontsize',14);title('fac (ratio of obs/pred) R');grid on;xlim([0 360])
subplot(312)
plot(cat(1,AZ{ik}),log10(facT(ik)),'o');set(gca,'Fontsize',14);title('fac (ratio of obs/pred) T');grid on;xlim([0 360])
subplot(313)
plot(cat(1,AZ{ik}),log10(facZ(ik)),'o');set(gca,'Fontsize',14);title('fac (ratio of obs/pred) Z');grid on;xlim([0 360])


figure(3)
subplot(311)
hist((CCR(ik)));set(gca,'Fontsize',14);title('CC R');grid on;xlim([-1 1])
subplot(312)
hist((CCT(ik)));set(gca,'Fontsize',14);title('CC T');grid on;xlim([-1 1])
subplot(313)
hist((CCZ(ik)));set(gca,'Fontsize',14);title('CC Z');grid on;xlim([-1 1])

figure

subplot(311)
plot(cat(1,AZ{ik}),CCR(ik),'o');set(gca,'Fontsize',14);title('CC R');grid on;xlim([0 360])
subplot(312)
plot(cat(1,AZ{ik}),CCT(ik),'o');set(gca,'Fontsize',14);title('CC T');grid on;xlim([0 360])
subplot(313)
plot(cat(1,AZ{ik}),CCZ(ik),'o');set(gca,'Fontsize',14);title('CC Z');grid on;xlim([0 360])

pause
close all

mkdir('plots')
for ii=1:length(ik)
    ir=ik(ii);  
    
    figure(1)
    subplot(311)
    max1=max(abs(Rdata{ir}(N)));
    plot(N*dt,Rdata{ir}(N),'b',N*dt,G{ir}(1,N)*max1/max(abs(G{ir}(1,N))),'r');
    title([char(sta(ir)) '  R and RR']);xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    
    subplot(312)
    plot(N*dt,Rdata{ir}(N),'b',N*dt,G{ir}(7,N)*max1/max(abs(G{ir}(7,N))),'r');
    title('R and ZR');xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    
    subplot(313)
    %% here we normalize the prediction with the mean of the ratios between predicted and observed amplitude.
    % by looking at the mean of the values, it effectively calibrate all
    % the stations to one factor.
    plot(N*dt,Rdata{ir}(N),'b',N*dt,R{ir}(N)*10.^mean(log10(facR(ik))),'r');
    title('R and modeled R');xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    print('-dpsc',['plots/HEC_Radial_' char(sta(ir)) '.ps'])
    
    
    figure(2)
    subplot(311)
    max1=max(abs(Zdata{ir}(N)));
    plot(N*dt,Zdata{ir}(N),'b',N*dt,G{ir}(3,N)*max1/max(abs(G{ir}(3,N))),'r');
    title([char(sta(ir)) '  Z and RZ']);xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    
    subplot(312)
    plot(N*dt,Zdata{ir}(N),'b',N*dt,G{ir}(9,N)*max1/max(abs(G{ir}(9,N))),'r');
    title('Z and ZZ');xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    
    subplot(313)
    plot(N*dt,Zdata{ir}(N),'b',N*dt,Z{ir}(N)*10.^mean(log10(facZ(ik))),'r');
    title('Z and modeled Z');xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    print('-dpsc',['plots/HEC_Vertical_' char(sta(ir)) '.ps'])
    
    
    figure(3)
    subplot(211)
    max1=max(abs(Tdata{ir}(N)));
    plot(N*dt,Tdata{ir}(N),'b',N*dt,G{ir}(5,N)*max1/max(abs(G{ir}(5,N))),'r');
    title([char(sta(ir)) '  T and TT']);xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    
    subplot(212)
    plot(N*dt,Tdata{ir}(N),'b',N*dt,T{ir}(N)*10.^mean(log10(facT(ir))),'r');
    title('T and modeled T');xlim([0 250]);grid on;set(gca,'Fontsize',14);xlabel('Time (s)')
    print('-dpsc',['plots/HEC_Transverse_' char(sta(ir)) '.ps'])
    
    disp(sta(ir))
    
    
    disp(['correlation coefficient ' num2str(CCR(ir))])
    pause
end

