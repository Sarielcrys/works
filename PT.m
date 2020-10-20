% PT-symmetric of double rings/one bus waveguide system
% Source code by CHEN Nuo
% 2020/10/13
% -------------------------------------------------------------------------------------------------------------------------------------------
%%
clear all
clc
close all

phase_detuning = linspace(-pi/12,pi/12,1001);                               % scanning round-trip phase detuning

% basic parameters
radius = 10e-6;
c = 3e8;
lambda0 = 1550e-9;

n = lambda0*680/2/pi/radius;                                                % m = 680, n = 2.0969
omega0 = 2*pi*c/lambda0;
L = 2*pi*radius;                                                            % assuming 2 rings have equal radius
T = L/(c/n);

% position 1 (bus waveguide/ring 1)
k1 = 0.2;
r1 = sqrt(1-k1^2);

% ring 1
a1 = 0.98;                                                                  % assuming ring 1 is negtive
alpha1 = -2*log(a1)/L;

% position 2/ring 2 (coulping between 2 rings)
k2 = linspace(0.001,0.3,1001);
r2 = sqrt(1-k2.^2);
k2fixed = k2(500);
r2fixed = sqrt(1-k2fixed^2);
a2 = linspace(0.4,2.0,1001);
%a2fixed = a2(339);                                                         % condition of PT-symmetic
a2fixed = a2(401);
alpha2 = -2*log(a2)/L;
alpha2fixed = -2*log(a2fixed)/L;

%%
trans1 = zeros(length(phase_detuning),length(a2));

for ii = 1:length(phase_detuning)                                           % kappa2 fixed, scanning a2 (i.e. alpha2，loss/gain)
    for jj = 1:length(a2)
        tao1(jj) = (r2fixed-a2(jj)*exp(-1i*phase_detuning(ii)))/(1-a2(jj)*r2fixed*exp(-1i*phase_detuning(ii)));
        TR = power(abs((r1-a1*tao1(jj)*exp(-1i*phase_detuning(ii)))/(1-a1*r1*tao1(jj)*exp(-1i*phase_detuning(ii)))),2);
        trans1(jj,ii) = TR;
    end
end

%%
for ii = 1:length(phase_detuning)                                           % adjust resolution

    for jj = 1:length(a2)
        if trans1(jj,ii) >= 100
            trans1(jj,ii) = 100;
        end
        if trans1(jj,ii) <= 1e-2
            trans1(jj,ii) = 1e-2;
        end
    end
end

%%
figure(1)
mesh(phase_detuning,a2,10*log10(trans1))
title('Transmission spectrum as the function of gain(k2 = 1.032)')
colormap jet
colorbar
set(gca,'Xtick',(-pi/12:pi/24:pi/12))
set(gca,'XtickLabel',{'-π/12','-π/24','0','π/24','π/12'})
xlabel('Round-trip phase detuning Δφrt (rad)')
set(gca,'Ytick',(0.95:0.05:1.2))
set(gca,'YtickLabel',{'0.95','1.00','1.05','1.10','1.15','1.20'})
ylabel('Round-trip attenuation factor a2')
zlabel('Transmission (dB)')

%%
trans2 = trans1;
for ii = 1:length(phase_detuning)
    for jj = 1:length(r2)
        tao2(jj) = (r2(jj)-a2fixed*exp(-1i*phase_detuning(ii)))/(1-a2fixed*r2(jj)*exp(-1i*phase_detuning(ii)));
        TR = power(abs((r1-a1*tao2(jj)*exp(-1i*phase_detuning(ii)))/(1-a1*r1*tao2(jj)*exp(-1i*phase_detuning(ii)))),2);
        trans2(jj,ii) = TR;
    end
end

%%
for ii = 1:length(phase_detuning)
    for jj = 1:length(r2)
        if trans2(jj,ii) >= 10
            trans2(jj,ii) = 10;
        end
        if trans2(jj,ii) <= 1e-1
            trans2(jj,ii) = 1e-1;
        end
    end
end

%%
figure(2)
mesh(phase_detuning,k2,10*log10(trans2))
title('Transmission spectrum as the function of κ(a2 = 1.041,loss basicly equals gain)')
%title('Transmission spectrum as the function of κ(a2 = 1.079)')
%title('Transmission spectrum as the function of κ(a2 = 1.027)')
colormap jet
colorbar
set(gca,'Xtick',(-pi/12:pi/24:pi/12))
set(gca,'XtickLabel',{'-π/12','-π/24','0','π/24','π/12'})
xlabel('Round-trip phase detuning Δφrt (rad)')
set(gca,'Ytick',(0:0.05:0.2))
set(gca,'YtickLabel',{'0','0.05','0.10','0.15','0.20'})
ylabel('κ')
zlabel('Transmission (dB)')

%% energy coupling formalism 
%% scanning loss/gain
% ring 1
y1 = alpha1*c/n;
kappa1 = k1/T;
gamma1 = -2*log(r1)/T;

% ring 2
y2 = alpha2*c/n;
kappa2 = k2fixed/T;
gamma2 = -2*log(r2)/T;

omega1 = omega0-1i*(y1+gamma1+y2)/4+sqrt(16*kappa2^2-(y1+gamma1-y2).^2)/4;
omega2 = omega0-1i*(y1+gamma1+y2)/4-sqrt(16*kappa2^2-(y1+gamma1-y2).^2)/4;

figure(3)
plot(a2,(real(omega1)-omega0)/1e9,'b.')
hold on
plot(a2,(real(omega2)-omega0)/1e9,'r.')
title('\omega real')
xlabel('gain (a2)')
ylabel('\omega - \omega0 (GHz)')

figure(4)
plot(a2,imag(omega1)/1e9,'g.')
hold on
plot(a2,imag(omega2)/1e9,'y.')
title('\omega image')
xlabel('gain (a2)')
ylabel('(GHz)')

%% scanning kappa
% ring 1
y1k = alpha1*c/n;
kappa1k = k1/T;
gamma1k = -2*log(r1)/T;

% ring 2
y2k = alpha2fixed*c/n;
kappa2k = k2/T;
gamma2k = -2*log(r2)/T;

omega1k = omega0-1i*(y1k+gamma1k+y2k)/4+sqrt(16*kappa2k.^2-(y1k+gamma1k-y2k).^2)/4;
omega2k = omega0-1i*(y1k+gamma1k+y2k)/4-sqrt(16*kappa2k.^2-(y1k+gamma1k-y2k).^2)/4;

figure(5)
plot(k2,(real(omega1k)-omega0)/1e9,'b.')
hold on
plot(k2,(real(omega2k)-omega0)/1e9,'r.')
title('\omega real')
xlabel('\kappa')
ylabel('\omega - \omega0 (GHz)')

figure(6)
plot(k2,imag(omega1k)*T,'g.')
hold on
plot(k2,imag(omega2k)*T,'y.')
title('\omega image')
xlabel('\kappa')

%% compare

figure(7)
plot3((real(omega1)-omega0)*T,a2,10*log10(trans1),'b.')
hold on
plot3((real(omega2)-omega0)*T,a2,10*log10(trans1),'r.')
mesh(phase_detuning,a2,10*log10(trans1))
title('两环损耗等于增益时，扫描\gamma')
colormap jet
colorbar

figure(8)
plot3((real(omega1k)-omega0)*T,k2,10*log10(trans2),'b.')
hold on
plot3((real(omega2k)-omega0)*T,k2,10*log10(trans2),'r.')
mesh(phase_detuning,k2,10*log10(trans2))
colormap jet
colorbar

