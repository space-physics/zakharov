%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D Langmuir turbulence simulation
%% according to Guio and Forme 2006 paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function Sim_v6_3_Linux(doWriteResults,newpath,useSinglePrec)
clc, close all

% if ~matlabpool('size')
%     matlabpool('local',4)
% end

if nargin < 1
    doWriteResults = false;
end

if nargin < 2 || isempty(newpath)
newpath='/tmp/test';    useSinglePrec = false;
end

if nargin < 3
    useSinglePrec = false;
end

noCaxis = true;


for Nnbeam=1:1
for Nvbeam=1:1
    close all;
    instantiation = [newpath,': n',num2str(Nnbeam,'%03d'),'_v',num2str(Nvbeam,'%03d')];
    
    if doWriteResults
     mkdir([newpath,filesep, 'results'], ['n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d')]);
    end


fnParam = [newpath,filesep,'parameters_n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d') '.bin'];
P = Langmuir_read_params(fnParam);
N = P.N;

fnOutput = [newpath,filesep,'output1_n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d') '.bin'];
output1=Langmuir_read_output(fnOutput);


omega_off=P.omegae+2*pi*300000;

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


%clear('output1','params')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      plot over time for all k to check the steady state

%kk_ss=[756 876 1100];
kk_ss=(1:length(k));
[aa, bb]=size(kk_ss);
EE_ss=zeros(bb,P.TT_res);
nn_ss=zeros(bb,P.TT_res);

taxis=(1:P.TT_res)*P.Tstep*P.res*1000;
for ii=1:P.QW


fnNN=[newpath,filesep,'nn' num2str(P.SEED,'%03d'), num2str(ii,'%03d'), '_n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d') '.bin'];
fnEE=[newpath,filesep,'EE' num2str(P.SEED,'%03d'), num2str(ii,'%03d'), '_n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d') '.bin'];

[nnComplex,EEComplex] = Langmuir_read_data(fnNN,fnEE,useSinglePrec);
nn2=abs(nnComplex).^2;
EE2=abs(EEComplex).^2;


  for jj=1:bb
      
    %for ttt=1:TT_res
        %baseInd1 = (ttt-1)*N*2+(kk_ss(jj)-1)*2+1;
        %baseInd2 = baseInd1+1;  
        %buffer1(1,ttt)=abs(EE(baseInd1)+1j*EE(baseInd2)).^2;
        %buffer2(1,ttt)=abs(nn(baseInd1)+1j*nn(baseInd2)).^2;
        %{
        baseInd = (ttt-1)*N + (kk_ss(jj)-1) + 1;
        buffer1(1,ttt) = EE2(baseInd);
        buffer2(1,ttt) = nn2(baseInd);
        %}
    %end

    ttt = 1:P.TT_res;
    baseInd = (ttt-1).*N + kk_ss(jj);
    
    EE_ss(jj,:)=EE_ss(jj,:) +EE2(baseInd); %+buffer1(1,:);    
    nn_ss(jj,:)=nn_ss(jj,:) +nn2(baseInd); %+buffer2(1,:);
    %figure(2*jj-1);
    %hold on;
    %plot(taxis,buffer1(1,:),'b');
    %figure(2*jj);
    %hold on;'E_kt'
    %plot(taxis,buffer2(1,:),'b');
    if ~mod(jj,500), disp(['jj/bb: ',num2str(jj/bb*100,'%0.2f'),'% complete']), end
  end


%floor(ii/QW*100)

end
EE_ss=EE_ss/P.QW;
nn_ss=nn_ss/P.QW;

 

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
title([instantiation,'  |E|^2 (dB)'])
colorbar;
if doWriteResults
 saveas(gcf,[newpath,filesep,'results',filesep,'n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d'),filesep, 'E_kt'],'fig');
end
figure;
imagesc(taxis,k,10*log10(abs(nn_ss)))
xlabel('Time (ms)');
ylabel('Wave number');
title([instantiation,'|n|^2 (dB)'])
colorbar;
if doWriteResults
 saveas(gcf,[newpath,filesep,'results',filesep,'n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d'),filesep,'n_kt'],'fig');
end

figure(100)
X1=linspace(0,1/3*sqrt(P.me/P.mi)/P.lambdaD,100);
Y1=(X1*P.lambdaD).^2;
% figure(101);
semilogy(X1,Y1,'k','LineWidth',3)
hold on
% semilogy(-X1,Y1,'k','LineWidth',3);
X2=X1(end)+linspace(0,k(end),1000);
Y2=X2*P.lambdaD*sqrt(P.me/P.mi); Y2=Y2-Y2(1)+Y1(end);
semilogy(X2,Y2,'k','LineWidth',3);
% semilogy(-X2,Y2,'k','LineWidth',3);
Y2=(X2*P.lambdaD).^4*P.mi/P.me; Y2=Y2-Y2(1)+Y1(end);
semilogy(X2,Y2,'k','LineWidth',3)
% semilogy(-X2,Y2,'k','LineWidth',3);
Y2=Y1(end)*ones(1,length(X1));
semilogy(X1,Y2,'k','LineWidth',3)
% semilogy(-X1,Y2,'k','LineWidth',3);
Y2=linspace(0,1,length(X1))*Y1(end);
X2=ones(1,length(X1))*X1(end);
semilogy(X2,Y2,'k','LineWidth',3)
% semilogy(-X2,Y2,'k','LineWidth',3);
hold off
title([instantiation])

resolution=1.0; %in ms
resolution=length(taxis)/(taxis(end)/resolution);
[EE_max, kk_max]=max(EE_ss);
W=.5*P.epsilon0*(EE_max(resolution:resolution:end))/(P.n0*P.Kb*P.Te);
figure(100)
hold on
semilogy(abs(k(kk_max(resolution:resolution:end))),W)
semilogy(abs(k(kk_max(resolution:resolution:end))),W,'b*')
hold off
xlabel('Wave number (m^-^1)')
ylabel('Electrostatic to thermal energy ration')
if doWriteResults
 saveas(gcf,[newpath,filesep,'results',filesep,'n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d'),filesep,'regiems'],'fig');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%resolution=0.5; %in ms
%resolution=length(taxis)/(taxis(end)/resolution);
%line=zeros(1,length(k));
%W2=zeros(1,length(k));

%for iio=0:18
%    line=EE_ss(:,iio*resolution+10);
%    W2=.5*epsilon0*line/(n0*Kb*Te);
    
%figure;
%X1=linspace(0,1/3*sqrt(me/mi)/lambdaD,100);
%Y1=(X1*lambdaD).^2;
%% figure(101);
%semilogy(X1,Y1,'k','LineWidth',3);
%hold on;
%% semilogy(-X1,Y1,'k','LineWidth',3);
%X2=X1(end)+linspace(0,k(end),1000);
%Y2=X2*lambdaD*sqrt(me/mi); Y2=Y2-Y2(1)+Y1(end);
%semilogy(X2,Y2,'k','LineWidth',3);
%% semilogy(-X2,Y2,'k','LineWidth',3);
%Y2=(X2*lambdaD).^4*mi/me; Y2=Y2-Y2(1)+Y1(end);
%semilogy(X2,Y2,'k','LineWidth',3);
%% semilogy(-X2,Y2,'k','LineWidth',3);
%Y2=Y1(end)*ones(1,length(X1));
%semilogy(X1,Y2,'k','LineWidth',3);
%% semilogy(-X1,Y2,'k','LineWidth',3);
%Y2=linspace(0,1,length(X1))*Y1(end);
%X2=ones(1,length(X1))*X1(end);
%semilogy(X2,Y2,'k','LineWidth',3);
%% semilogy(-X2,Y2,'k','LineWidth',3);
%semilogy(abs(k),(W2),'r*');
%Xlabel('Wave number (m^-^1)');
%ylabel('Electrostatic to thermal energy ration');
%title(['Time: ', num2str(taxis(iio*resolution+10)),' (ms)']);
%hold off;
%end
%close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   average over realizations to calculate the spectrum for the last 0.5ms, for all k

    Taverage=(.5e-3);
    OO=80; %how many periods?
    ACF_length = floor(Taverage/(P.Tstep*P.res));

    ACF_E=zeros(OO,N,2*ACF_length-1);
    ACF_n=zeros(OO,N,2*ACF_length-1);
    Spec_E=zeros(N,2*ACF_length-1);
    Spec_n=zeros(N,2*ACF_length-1);
    taxis=(1:P.TT_res)*P.Tstep*P.res*1000;
    
for ii=1:P.QW
%   fnNN=[newpath,filesep,'nn' num2str(SEED+ii-1) '_n' num2str(Nnbeam) '_v' num2str(Nvbeam) '.bin'];
%   fnEE=[newpath,filesep,'EE' num2str(SEED+ii-1) '_n' num2str(Nnbeam) '_v' num2str(Nvbeam) '.bin'];
%   [nnComplex,EEComplex] = Langmuir_read_data(fnNN,fnEE,useSinglePrec);
   nnComplex = reshape(nnComplex,1,1,[]);    EEComplex = reshape(EEComplex,1,1,[]);
    
    TT_res=length(taxis);
%     TT_res=TT_res-ACF_length*30;

for oo=1:OO
 % for ii=1:QW

for kk=1:N
    %floor(kk/N*100)
    %{
    for ttt=TT_res-ACF_length:TT_res-1
        outInd = ttt-(TT_res-ACF_length)+1;
        baseInd = (ttt-1).*N + kk;
        AA(outInd) = EE2(baseInd);
        BB(outInd) = nn2(baseInd);
            %AA(outInd)=EE((ttt-1)*N*2+(kk-1)*2+1)+1j*EE((ttt-1)*N*2+(kk-1)*2+2);
            %BB(outInd)=nn((ttt-1)*N*2+(kk-1)*2+1)+1j*nn((ttt-1)*N*2+(kk-1)*2+2);
    end
    %}
    ttt=TT_res-ACF_length:TT_res-1;
    %outInd = ttt-(TT_res-ACF_length)+1;
    baseInd = (ttt-1).* N + kk;
    
    ACF_E(oo,kk,:)=ACF_E(oo,kk,:) +xcorr(EEComplex(baseInd));%+xcorr(AA);
    %Spec_E(kk,:)=fftshift(fft(ACF_E(kk,:)))/(ACF_length*2-1);
    ACF_n(oo,kk,:)=ACF_n(oo,kk,:)+ xcorr(nnComplex(baseInd));%+xcorr(BB);
    %Spec_n(kk,:)=fftshift(fft(ACF_n(kk,:)))/(ACF_length*2-1);
    
end
TT_res=TT_res-ACF_length;

end


%floor(ii/QW*100)
if ~mod(ii,1), display(['ii/QW: ',num2str(ii/P.QW*100,'%0.3f'),'% complete']),end
end
ACF_E=ACF_E/P.QW/(ACF_length-1)*2;
ACF_n=ACF_n/P.QW/(ACF_length-1)*2;



% thermal_pl=ACF_E(1,:,ACF_length);
% thermal_il=ACF_n(1,:,ACF_length);
% figure;
% plot(k,10*log10(thermal_pl));
% hold on;
% plot(k,10*log10(E_thermal_k_squared),'r');
% hold off;
% xlabel('Wave number (m^-^1)');
% ylabel('Intensity (dB)');
% title([instantiation,'|E|^2 thermal level'])
% figure;
% plot(k,10*log10(thermal_il));
% hold on;
% plot(k,10*log10(n_thermal_k_squared),'r');
% hold off;
% xlabel('Wave number (m^-^1)');
% ylabel('Intensity (dB)');
% title([instantiation,'|n|^2 thermal level'])


faxis=linspace(-1/(P.res*P.Tstep)/2,1/(P.res*P.Tstep)/2,ACF_length*2-1);
TT_res=length(taxis);
% TT_res=TT_res-ACF_length*30;
for oo=1:OO
  
for kk=1:N
    %floor(kk/N*100)
    Spec_E(kk,:)=fftshift(fft(ACF_E(oo,kk,:)))/(ACF_length*2-1)/(faxis(2)-faxis(1));
    Spec_n(kk,:)=fftshift(fft(ACF_n(oo,kk,:)))/(ACF_length*2-1)/(faxis(2)-faxis(1));
end


figure;
imagesc(-faxis,k,10*log10(abs(Spec_E)));
caxis([-80 -40]);
colorbar;
hold on;
plot((omegaL-omegae)/2/pi,k,'k');
xlabel('Frequency (Hz)');
ylabel('Wave number (m^-^1)');
title([instantiation,' "Spec_E"  Time: ', num2str(taxis(TT_res-ACF_length)), '--', num2str(taxis(TT_res)), ' (ms)']);
if noCaxis, else caxis([-80 -40]), end %#ok<UNRCH>
if doWriteResults
 nameWr=['specE--' num2str(oo)];
 saveas(gcf,[newpath,filesep,'results',filesep,'n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d'),filesep, nameWr],'fig');
end

figure;
imagesc(-faxis,k,10*log10(abs(Spec_n)));
caxis([65 110]);
colorbar;
hold on;
plot(abs(Cs*k/2/pi),k,'k');
xlabel('Frequency (Hz)');
ylabel('Wave number (m^-^1)');
title([instantiation,' "Spec_n" Time: ', num2str(taxis(TT_res-ACF_length)), '--', num2str(taxis(TT_res)), ' (ms)']);
if noCaxis, else caxis([80 120]), end %#ok<UNRCH>
if doWriteResults
 nameWr=['specn--' num2str(oo)];
 saveas(gcf,[newpath,filesep,'results',filesep,'n' num2str(Nnbeam,'%03d') '_v' num2str(Nvbeam,'%03d'),filesep, nameWr],'fig');
end

TT_res=TT_res-ACF_length;
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   DO ion- and pl-spectra calculations from model
