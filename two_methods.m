% Stimulation of ADMR ( we compare two methods of energy coupling and power coupling formalisms)
% Source code by CHEN Nuo
% 2020/10/09

% 仿真结果中，蓝色线是功率耦合模中直通端的响应，红色线是功率耦合模下载端响应
% 黄色线是时域耦合模中直通端响应，绿色线是时域耦合模中下载端响应
% ------------------------------------------------------------------------------------------------------------------------------------------
clear all
clc
close all

% base parameters
c = 3e8;                                                                    % speed of light
n = 2.05;                                                                   % effective group refractive index
R = 8e-6;                                                                   % radius of microring
lambda = 1550e-9;                                                           % resonance wavelength

% calculation
L_rt = 2*pi*R;                                                              % round-trip length
T_rt = L_rt/c/n;                                                            % round-trip time
phase_rt = (-pi:pi/100:3*pi);                                               % scanning round-trip phase
omega = phase_rt/T_rt;                                                      % relationship between phase and angular frequency
omega0 = 2*pi*c/(lambda*n);                                                 % resonance angular frequency

% ring(assuming lossless)
a = 1;                                                                      % round-trip field attenuation factor
y0 = 0;                                                                     % intrinsic decay rate gamma0

% coupler 1
k1 = sqrt(0.08);                                                            % field coupling coeffient kappa1
u1 = sqrt(k1^2/T_rt);                                                       % energy coupling coeffient mu1
y1 = u1^2/2;                                                                % coupling decay rate gamma1
t1 = sqrt(1-(k1)^2);                                                        % power transmission coeffient tau1

% coupler 2
k2 = sqrt(0.1);                                                             % same as above
u2 = sqrt(k2^2/T_rt);
y2 = u2^2/2;
t2 = sqrt(1-(k2)^2);                                                       

y = y0+y1+y2;                                                               % total energy decay
p0 = 1/(t1*t2*a);
z0 = t1/(t2*a);

% response
Through_power = zeros(1,length(phase_rt));                                  % matrix of through port using power coupling formalism
Drop_power = zeros(1,length(phase_rt));                                     % matrix of drop port using power coupling formalism
Through_energy = zeros(1,length(omega));                                    % matrix of through port using energy coupling formalism
Drop_energy = zeros(1,length(omega));                                       % matrix of drop port using energy coupling formalism
phase_responseD = zeros(1,length(omega));
phase_responseT = zeros(1,length(omega));

for ii = 1:length(phase_rt)

    Ht1 = (t1-t2*a*exp(-1i*phase_rt(ii)))/(1-t1*t2*a*exp(-1i*phase_rt(ii))); % transfer functions
    Hd1 = -(k1*k2*sqrt(a)*exp(-1i*phase_rt(ii)/2))/(1-t1*t2*a*exp(-1i*phase_rt(ii)));    
    
    Hd2 = -(u1*u2)/((1i*phase_rt(ii)/T_rt)+y);
    Ht2 = (1i*(phase_rt(ii)/T_rt)+y-u1^2)/(1i*(phase_rt(ii)/T_rt)+y);

    Tt1 = (abs(Ht1))^2;                                                     % transmission spectrums
    Td1 = (abs(Hd1))^2;

    Tt2 = (abs(Ht2))^2;
    Td2 = (abs(Hd2))^2;
    
    Pd = angle(Hd1);                                                        % phase responses
    Pt = angle(Ht1);
    
    Through_power(ii) = Tt1;
    Drop_power(ii) = Td1;
    Through_energy(ii) = Tt2;
    Drop_energy(ii) = Td2;
    phase_responseD(ii) = Pd ;
    phase_responseT(ii) = Pt;
end

% dessin
figure(1)                                                                 
box on;
plot(phase_rt,10*log10(Drop_power),'r');
hold on;
plot(phase_rt,10*log10(Through_power),'b');
hold on;
plot(phase_rt,10*log10(Through_energy),'y--');
hold on;
plot(phase_rt,10*log10(Drop_energy),'g--');

set(gca,'XTick',(-pi:pi:3*pi));                                         
set(gca,'XtickLabel',{'-π','0','π','2π','3π'});
xlabel('Round-trip phase detuning Δφrt (rad)');

ylabel('Transmission(dB)');

title('Comparing of two methods');                                          
legend({'D power','T power','D energy','T energy'},'location','southwest');

figure(2)
box on;
plot(phase_rt,phase_responseD);
hold on;
plot(phase_rt,phase_responseT);

set(gca,'XTick',(-pi:pi:3*pi));                                         
set(gca,'XtickLabel',{'-π','0','π','2π','3π'});
xlabel('Round-trip phase detuning Δφrt (rad)');
set(gca,'YTick',(-pi:pi/2:2*pi));                                         
set(gca,'YtickLabel',{'-π','-0.5π','0','0.5π','π','1.5π','2π'});
ylabel('phase responses (rad)');

title('Phase Response at Through and Drop Ports');
legend({'Drop port','Through port'},'location','north');
