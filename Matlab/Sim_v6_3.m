%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D Langmuir turbulence simulation
%% according to Guio and Forme 2006 paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
%clc;
addpath('D:\UNI\research\Langmuir simulation\LangmuirSim_v5\Langmuir_v6.3\Langmuir_v5\test17')
addpath('D:\UNI\research\Langmuir simulation\Faddeeva-MATLAB')
fclose('all');
fid3=fopen('parameters.bin');
params=fread(fid3,'float64');
fid4=fopen('output1.bin');
output1=fread(fid4,'float64');
pi=params(1);
me=params(2);
electroncharge=params(3);
mi=params(4);
Kb=params(5);
eV=params(6);
epsilon0=params(7);
Z=params(8);
Te=params(9);
Ti=params(10);
nuic=params(11);
nuec=params(12);
n0=params(13);
nbeam=params(14);
vbeam_ev=params(15);
vbeam=params(16);
tetabeam=params(17);
endTime=params(18);
Tstep=params(19);
TT=params(20);
res=params(21);
TT_res=params(22);
L=params(23);
N=params(24);
Xstep=params(25);
QW=params(26);
SEED=params(27);
eta=params(28);
ve=params(29);
Cs=params(30);
omegae=params(31);
lambdaD=params(32);
j=sqrt(-1);
omega_off=omegae+2*pi*300000;

p=zeros(1,N);
k=zeros(1,N);
Xsection_ion=zeros(1,N);
Xsection_pl=zeros(1,N);
E_thermal_k_squared=zeros(1,N);
n_thermal_k_squared=zeros(1,N);
omegaL=zeros(1,N);
gamas=zeros(1,N);
nui=zeros(1,N);
nue=zeros(1,N);
Source_factor_E=zeros(1,N);
Source_factor_n=zeros(1,N);


for i=1:N
    p(i)=output1((i-1)*12+1);
    k(i)=output1((i-1)*12+2);
    Xsection_ion(i)=output1((i-1)*12+3);
    Xsection_pl(i)=output1((i-1)*12+4);
    E_thermal_k_squared(i)=output1((i-1)*12+5);
    n_thermal_k_squared(i)=output1((i-1)*12+6);
    omegaL(i)=output1((i-1)*12+7);
    gamas(i)=output1((i-1)*12+8);
    nui(i)=output1((i-1)*12+9);
    nue(i)=output1((i-1)*12+10);
    Source_factor_E(i)=output1((i-1)*12+11);
    Source_factor_n(i)=output1((i-1)*12+12);
end

fclose(fid3);
fclose(fid4);
clear output1;
clear params;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      plot over time for all k to check the steady state

%kk_ss=[756 876 1100];
kk_ss=(1:length(k));
[aa bb]=size(kk_ss);
EE_ss=zeros(bb,TT_res);
nn_ss=zeros(bb,TT_res);
buffer1=zeros(1,TT_res);
buffer2=zeros(1,TT_res);
taxis=(1:TT_res)*Tstep*res*1000;
for ii=1:QW
    fid=fopen(['nn' num2str(SEED+ii-1) '.bin']);
    nn=fread(fid,'float64');
    fid2=fopen(['EE' num2str(SEED+ii-1) '.bin']);
    EE=fread(fid2,'float64');

    for jj=1:bb
        for ttt=1:TT_res
            buffer1(1,ttt)=abs(EE((ttt-1)*N*2+(kk_ss(jj)-1)*2+1)+j*EE((ttt-1)*N*2+(kk_ss(jj)-1)*2+2)).^2;
            buffer2(1,ttt)=abs(nn((ttt-1)*N*2+(kk_ss(jj)-1)*2+1)+j*nn((ttt-1)*N*2+(kk_ss(jj)-1)*2+2)).^2;
        end
        EE_ss(jj,:)=EE_ss(jj,:)+buffer1(1,:);    
        nn_ss(jj,:)=nn_ss(jj,:)+buffer2(1,:);
        %figure(2*jj-1);
        %hold on;
        %plot(taxis,buffer1(1,:),'b');
        %figure(2*jj);
        %hold on;
        %plot(taxis,buffer2(1,:),'b');
    end
    clear nn;
    clear EE;
    fclose(fid);
    fclose(fid2);
floor(ii/QW*100)

end
EE_ss=EE_ss/QW;
nn_ss=nn_ss/QW;

 

% for jj=1:bb
% figure(2*jj-1);
% hold on;
% plot(taxis,sqrt(EE_ss(jj,:)),'b');
% xlabel('Time (ms)');
% ylabel('Field Intensity (V/m)');
% title(['E (k=' num2str(k(kk_ss(1,jj))) ')']);
% figure(2*jj);
% hold on;
% plot(taxis,sqrt(nn_ss(jj,:)),'b');
% xlabel('Time (ms)');
% ylabel('Density Fluctuations (linear scale)');
% title(['n (k=' num2str(k(kk_ss(1,jj))) ')']);
% end

figure;
imagesc(taxis,k,10*log10(abs(EE_ss)))
xlabel('Time (ms)');
ylabel('Wave number');
title('|E|^2 (dB)')
colorbar;
figure;
imagesc(taxis,k,10*log10(abs(nn_ss)))
xlabel('Time (ms)');
ylabel('Wave number');
title('|n|^2 (dB)')
colorbar;
        
figure(100);
X1=linspace(0,1/3*sqrt(me/mi)/lambdaD,100);
Y1=(X1*lambdaD).^2;
% figure(101);
semilogy(X1,Y1,'k','LineWidth',3);
hold on;
% semilogy(-X1,Y1,'k','LineWidth',3);
X2=X1(end)+linspace(0,k(end),1000);
Y2=X2*lambdaD*sqrt(me/mi); Y2=Y2-Y2(1)+Y1(end);
semilogy(X2,Y2,'k','LineWidth',3);
% semilogy(-X2,Y2,'k','LineWidth',3);
Y2=(X2*lambdaD).^4*mi/me; Y2=Y2-Y2(1)+Y1(end);
semilogy(X2,Y2,'k','LineWidth',3);
% semilogy(-X2,Y2,'k','LineWidth',3);
Y2=Y1(end)*ones(1,length(X1));
semilogy(X1,Y2,'k','LineWidth',3);
% semilogy(-X1,Y2,'k','LineWidth',3);
Y2=linspace(0,1,length(X1))*Y1(end);
X2=ones(1,length(X1))*X1(end);
semilogy(X2,Y2,'k','LineWidth',3);
% semilogy(-X2,Y2,'k','LineWidth',3);
hold off;

resolution=1.0; %in ms
resolution=length(taxis)/(taxis(end)/resolution);
[EE_max kk_max]=max(EE_ss);
W=.5*epsilon0*(EE_max(resolution:resolution:end))/(n0*Kb*Te);
figure(100);
hold on;
semilogy(abs(k(kk_max(resolution:resolution:end))),(W));
semilogy(abs(k(kk_max(resolution:resolution:end))),(W),'b*');
hold off;
Xlabel('Wave number (m^-^1)');
ylabel('Electrostatic to thermal energy ration');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resolution=0.5; %in ms
resolution=length(taxis)/(taxis(end)/resolution);
line=zeros(1,length(k));
W2=zeros(1,length(k));

for iio=0:18
    line=EE_ss(:,iio*resolution+10);
    W2=.5*epsilon0*line/(n0*Kb*Te);
    
figure;
X1=linspace(0,1/3*sqrt(me/mi)/lambdaD,100);
Y1=(X1*lambdaD).^2;
% figure(101);
semilogy(X1,Y1,'k','LineWidth',3);
hold on;
% semilogy(-X1,Y1,'k','LineWidth',3);
X2=X1(end)+linspace(0,k(end),1000);
Y2=X2*lambdaD*sqrt(me/mi); Y2=Y2-Y2(1)+Y1(end);
semilogy(X2,Y2,'k','LineWidth',3);
% semilogy(-X2,Y2,'k','LineWidth',3);
Y2=(X2*lambdaD).^4*mi/me; Y2=Y2-Y2(1)+Y1(end);
semilogy(X2,Y2,'k','LineWidth',3);
% semilogy(-X2,Y2,'k','LineWidth',3);
Y2=Y1(end)*ones(1,length(X1));
semilogy(X1,Y2,'k','LineWidth',3);
% semilogy(-X1,Y2,'k','LineWidth',3);
Y2=linspace(0,1,length(X1))*Y1(end);
X2=ones(1,length(X1))*X1(end);
semilogy(X2,Y2,'k','LineWidth',3);
% semilogy(-X2,Y2,'k','LineWidth',3);
semilogy(abs(k),(W2),'r*');
Xlabel('Wave number (m^-^1)');
ylabel('Electrostatic to thermal energy ration');
title(['Time: ', num2str(taxis(iio*resolution+10)),' (ms)']);
hold off;
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   average over realizations to calculate the spectrum for the last 0.5ms, for all k

    Taverage=(1e-3);
    OO=4; %how many periods?
    AA=zeros(1,1,floor(Taverage/(Tstep*res)));
    BB=zeros(1,1,floor(Taverage/(Tstep*res)));
    ACF_E=zeros(OO,N,2*length(AA)-1);
    ACF_n=zeros(OO,N,2*length(BB)-1);
    Spec_E=zeros(N,2*length(AA)-1);
    Spec_n=zeros(N,2*length(BB)-1);
    taxis=(1:TT_res)*Tstep*res*1000;
    
for ii=1:QW
    fid=fopen(['nn' num2str(SEED+ii-1) '.bin']);
    nn=fread(fid,'float64');
    fid2=fopen(['EE' num2str(SEED+ii-1) '.bin']);
    EE=fread(fid2,'float64');
    TT_res=length(taxis);
%     TT_res=TT_res-length(AA)*30;
for oo=1:OO
 % for ii=1:QW

for kk=1:N
    %floor(kk/N*100)
    for ttt=TT_res-length(AA):TT_res-1
            AA(ttt-(TT_res-length(AA))+1)=EE((ttt-1)*N*2+(kk-1)*2+1)+j*EE((ttt-1)*N*2+(kk-1)*2+2);
            BB(ttt-(TT_res-length(AA))+1)=nn((ttt-1)*N*2+(kk-1)*2+1)+j*nn((ttt-1)*N*2+(kk-1)*2+2);
    end
    ACF_E(oo,kk,:)=ACF_E(oo,kk,:)+xcorr(AA);
    %Spec_E(kk,:)=fftshift(fft(ACF_E(kk,:)))/(length(AA)*2-1);
    ACF_n(oo,kk,:)=ACF_n(oo,kk,:)+xcorr(BB);
    %Spec_n(kk,:)=fftshift(fft(ACF_n(kk,:)))/(length(AA)*2-1);
end
TT_res=TT_res-length(AA);

end
clear EE;
clear nn;
floor(ii/QW*100)
end
ACF_E=ACF_E/QW/(length(AA)-1)*2;
ACF_n=ACF_n/QW/(length(AA)-1)*2;



thermal_pl=ACF_E(1,:,length(AA));
thermal_il=ACF_n(1,:,length(AA));
figure;
plot(k,10*log10(thermal_pl));
hold on;
plot(k,10*log10(E_thermal_k_squared),'r');
hold off;
xlabel('Wave number (m^-^1)');
ylabel('Intensity (dB)');
title('|E|^2 thermal level');
figure;
plot(k,10*log10(thermal_il));
hold on;
plot(k,10*log10(n_thermal_k_squared),'r');
hold off;
xlabel('Wave number (m^-^1)');
ylabel('Intensity (dB)');
title('|n|^2 thermal level');


faxis=linspace(-1/(res*Tstep)/2,1/(res*Tstep)/2,length(AA)*2-1);
TT_res=length(taxis);
% TT_res=TT_res-length(AA)*30;
for oo=1:OO
  
for kk=1:N
    floor(kk/N*100)
    Spec_E(kk,:)=fftshift(fft(ACF_E(oo,kk,:)))/(length(AA)*2-1)/(faxis(2)-faxis(1));
    Spec_n(kk,:)=fftshift(fft(ACF_n(oo,kk,:)))/(length(AA)*2-1)/(faxis(2)-faxis(1));
end


figure;
imagesc(-faxis,k,log10(abs(Spec_E)));
colorbar;
hold on;
plot((omegaL-omegae)/2/pi,k,'k');
xlabel('Frequency (Hz)');
ylabel('Wave number (m^-^1)');
title(['Time: ', num2str(taxis(TT_res-length(AA))), '--', num2str(taxis(TT_res)), ' (ms)']);
caxis([-8 -4]);

figure;
imagesc(-faxis,k,log10(abs(Spec_n)));
colorbar;
hold on;
plot(abs(Cs*k/2/pi),k,'k');
xlabel('Frequency (Hz)');
ylabel('Wave number (m^-^1)');
title(['Time: ', num2str(taxis(TT_res-length(AA))), '--', num2str(taxis(TT_res)), ' (ms)']);
caxis([8 12]);
TT_res=TT_res-length(AA);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   DO ion- and pl-spectra calculations from model

fmi=mi/(1.6726e-27);
kk=1100;  % choose the wave number
radarf=abs(k(kk))*3*10^8/4/pi;
dt=4e-5/(radarf/(500*10^6))/700;
Ni=max(60,60*(radarf/(500*10^6)))*3000;
[Sion wion]=sheffield(radarf/10^6,Te,Ti,n0,fmi,Z,dt,Ni,0);



figure;
% loglog(w/2/pi/1000,S,'b');
plot(wion/2/pi,Sion*2*pi*n0,'b');
xlabel('Frequency (Hz)');
ylabel('Ion-line PSD');

% 
% BW=3000000;  %in Hz
% ff1=omegaL(kk)/2/pi-BW/2;   %Hz
% df=2;   %Hz
% ff2=omegaL(kk)/2/pi+BW/2;   %Hz
% [Spl wpl]=sheffield_pl(radarf/10^6,Te,Ti,n0,fmi,Z,0,ff1,df,ff2);
% 
% 
% figure;
% plot(wpl/2/pi/10^6,Spl,'b');
% xlabel('Frequency (MHz)');
% ylabel('Plasma-line PSD');
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_star=2/3*Cs*omegae/ve^2;
E_thermal=sqrt(E_thermal_k_squared);
E_thereshold=sqrt(4*n0*Kb*Ti*nuec/omegae/0.58/epsilon0);

figure
plot(taxis,sqrt(EE_ss(441,:)),'g')
hold on
plot(taxis,sqrt(EE_ss(700,:)))
plot(taxis,sqrt(EE_ss(668,:)),'k')
plot(taxis,E_thermal(700)*ones(1,length(taxis)),'r')
plot(taxis,E_thereshold*ones(1,length(taxis)),'r')
ylabel('Electric field (V/m)');
xlabel('Time (ms)');


