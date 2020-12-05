%% This script uses the theodorsen theory for the flutter of a flat plate
% Author: Giampiero Bindolino, Nicola Fonzi, Vittorio Cavalieri
%%% Structural definition
clear all;close all;clc
index=1;
velocities = sqrt(1.4*287*273)*[0.1 0.2 0.3 0.357 0.364];
theo = cell2struct(cell(4, 1),{'U','t','h','alpha'},1);
plot_time = 0;
plot_fft = 0;
plot_eig = 1;
if plot_eig
    plot_time = 0;
    plot_fft = 0;
    Ma_0 = 0.1;
    Ma_end = 0.42;
    velocities = sqrt(1.4*287*273)*linspace(Ma_0,Ma_end,20);
    ww1 = [];
    ww2 = [];
end
for U = velocities
    c = 1;
    b = c/2;
    xf = c/4;
    xac= c/4;
    e = (xf-xac)/c; %distanza fra centro aerodinamico e asse elastico, normalizzato per la corda
    T = 273;
    Re = 4e6;
    muref = 1.716e-5;
    Tref = 273.15;
    SutherlandC = 110.4;
    visc_dyn= muref*(T/Tref)^1.5*((Tref+SutherlandC)/(T+SutherlandC));
    rho(index) = Re*visc_dyn/U/c;
    Ma(index) = U/sqrt(1.4*287*T);
    X = 0.25;
    ra = 0.5;
    wa = 45;
    w_segnato = 0.3185;
    wh = w_segnato*wa;
    mu = 100;
    M = mu*pi*rho(index)*b^2;
    I = ra^2*M*b^2;
    CI =  0;
    CM =  0;
    KI = wa^2*I;
    KM = wh^2*M;
    Stat = X*M*b;
    MS = [M -Stat ; -Stat I];
    CS = [CM 0 ; 0 CI];
    KS = [KM 0 ; 0 KI];
    AS = [zeros(2,2) eye(2) ; -inv(MS)*KS  -inv(MS)*CS];
    %%% Aerodynamic definition
    S=c;
    Pdyn = 0.5*rho(index)*U*U;
    BUS = .5*rho(index)*U^2*c*[ 0; 0; 2*pi ; 2*pi*(xf - xac)];
    MA = rho(index)*pi*c^2/4*[ -1               -(xf-.5*c) ; ...
        -(xf-.5*c) -(c^2/32+(xf-.5*c)^2)];
    CA = rho(index)*pi*c^2/4*U*[ 0  1 ;  ...
        0 -(.75*c-xf)];
    KA = zeros(2,2);
    MAE = MS-MA;
    CAE = CS-CA;
    KAE = KS-KA;
    A1 = .165;
    b1 = .0455;
    A2 = .335;
    b2 = .30;
    AAS = [0    1 ; -b1*b2*(2*U/c)^2  -(b1+b2)*(2*U/c)];
    BAS = [0;1];
    CAS = [b1*b2/2*(U/b)^2 (A1*b1+A2*b2)*(U/b)];
    DAS = [.5];
    AQS= [0   1   -1/U  (c*.75-xf)/U];
    LA = [eye(2)      zeros(2,2)  zeros(2,2); ...
        zeros(2,2)     MAE      zeros(2,2); ...
        zeros(2,2)  zeros(2,2)  eye(2)];
    RA = [zeros(2,2)  eye(2)  zeros(2,2); ...
        -KAE       -CAE   zeros(2,2) ; ...
        zeros(2,4)         AAS];
    RAA = [BUS*DAS*AQS BUS*CAS ;  BAS*AQS  zeros(2,2)];
    AT=inv(LA)*(RA+RAA);
    ZI=-inv(AAS)*BAS*5*pi/180;
    BT=inv(LA)*[0; 0; 0; KI ; 0; 0];DT=[ 0 0 .5*2*pi]';
    CT=[ 1 0 0 0 0 0 ; 0 1 0 0 0 0 ; AQS*pi 2*pi*b1*b2/2*(U/b)^2 2*pi*(A1*b1+A2*b2)*(U/b)];
    sa=ss(AT,BT,CT,DT);
    Fs=10000;
    [Y,T,X]=initial(sa,[0; 5/180*pi; 0 ;0; ZI],0:(1/Fs):10);
    h=-Y(:,1);
    alpha = Y(:,2);
    time = T;
    theo(index).U = velocities(index);
    theo(index).t = time;
    theo(index).h = h;
    theo(index).alpha = alpha;
    if plot_time
        figure
        plot(time,h)
        figure
        plot(time,alpha)
    end
    i0 = 4501;
    time = time(i0:end);
    h = h(i0:end);
    alpha = alpha(i0:end);
    L = length(h);
    H = fft(h);
    H = abs(H/L);
    H = H(1:floor(L/2)+1);
    H = 2*H;
    FreqVect = Fs*(0:floor(L/2))/L;
    if plot_fft
        figure
        plot(FreqVect,H)
        xlim([0,20]);
    end
    [pks,locs] = findpeaks(H,FreqVect);
    f_alpha(index)=locs(min([2,length(locs)]));
    f_h(index)=locs(1);
    L = length(alpha);
    A = fft(alpha);
    A = abs(A/L);
    A = A(1:floor(L/2)+1);
    A = 2*A;
    FreqVect = Fs*(0:floor(L/2))/L;
    if plot_fft
        figure
        plot(FreqVect,A)
        xlim([0,20]);
    end
    [pks,locs] = findpeaks(A,FreqVect);
    [~,ii] = min(abs(locs-7));
    f_alpha(index)=locs(ii);
    if plot_eig
        lambda = eig(AT);
        if any(real(lambda)>0)
            fprintf('Warning: eigenvalue has positive real part @ Ma=%.3f\n',Ma(end));
        end
        ww = abs(imag(lambda)/45);
        [~,i1] = min(abs(ww-0.3));
        [~,i2] = min(abs(ww-1.1));
        ww1 = [ww1; ww(i1)];
        ww2 = [ww2; ww(i2)];
    end
    index=index+1;
end
figure
%plot(velocities,f_alpha/wa*2*pi)
plot(Ma,f_alpha/wa*2*pi,'LineWidth',2)
hold on
%plot(velocities,f_h/wa*2*pi)
plot(Ma,f_h/wa*2*pi,'LineWidth',2)
ylim([0 1.2])
if plot_eig
    figure
    plot(Ma,ww2,'--','LineWidth',2)
    hold on
    plot(Ma,ww1,'--','LineWidth',2)
    ylim([0 1.2])
end