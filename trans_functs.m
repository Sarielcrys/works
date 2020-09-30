clear all;
close all;

phase_rt = (-pi:pi/100:pi);                                                 % ������������

% coupler 1 
k1 = sqrt(0.2);                                                             % �������ϵ��k1
t1 = sqrt(1-(k1)^2);                                                        % ������������£����崫��ϵ��t1
M1 = [t1,-1i*k1; -1i*k1,t1];                                                % �������

% ring
R = 8*1e-6;                                                                 % ���廷ǻ�뾶
L_rt = 2*pi*R;                                                              % �󻷳̳���
a = 0.98;                                                                   % ����͸��ϵ��a

% coupler 2
k2 = sqrt(0.2);                                                             % �������ϵ��k2
t2 = sqrt(1-(k2)^2);                                                        % ������������£����崫��ϵ��t2
M2 = [t2,-1i*k2; -1i*k2,t2];                                                % �������

% coupler 2 for critical coulping
% t2 = t1/a;                                                                  % �������ϵ��k2,�ٽ����������
% k2 = sqrt(1-(k1)^2);                                                        % ������������£����崫��ϵ��t2
% M2 = [t2 -1i*k2; -1i*k2 t2];                                                % �������

% response
Through = zeros(1,length(phase_rt));                                        % ֱͨ�˾���
Drop = zeros(1,length(phase_rt));                                           % ���ض˾���


for ii = 1:length(phase_rt)
    P = [0,sqrt(a)*exp(-1i*phase_rt(ii)/2); 1/(sqrt(a)*exp(-1i*phase_rt(ii)/2)),0]; % ΢���Ĵ������
    H = M1*P*M2;
    
    Ht = (t1-t2*a*exp(-1i*phase_rt(ii)))/(1-t1*t2*a*exp(-1i*phase_rt(ii))); %���ݺ�����
    Hd = -(k1*k2*sqrt(a)*exp(-1i*phase_rt(ii)/2))/(1-t1*t2*a*exp(-1i*phase_rt(ii)));    

    Tt = (abs(Ht))^2;                                                       % ���ʹ�����Ӧ
    Td = (abs(Hd))^2;
                                                                            % �ܴ������
    Through(ii) = Tt;
    Drop(ii) = Td;
end

figure(1);                                                                  % ��ͼ
box on;

plot(phase_rt,10*log10(Drop),'r');

hold on;

plot(phase_rt,10*log10(Through),'b');

set(gca,'XTick',(-pi:pi:pi));                                               % ����x����
set(gca,'XtickLabel',{'-��','0','��'});
xlabel('������λʧг��(����rt)');

set(gca,'YTick',[-30,-25,-20,-15,-10,-5,0,]);                               % ����y����
set(gca,'YtickLabel',{'-30','-25','-20','-15','-10','-5','0'});
ylabel('������(dB)');

title('ADMR�Ĺ��ʹ�����Ӧ,��1=��2,a=0.98');                                                 % ����ͼ������
legend('Drop Port','Through Port');                                         % ����ͼ��

