clear all;
close all;

phase_rt = (-0.3*pi:pi/200:0.3*pi);                                         % ������������

% coupler 1 
t1 = 0.95;                                                                  % ������������£����崫��ϵ��t1
k1 = sqrt(1-t1^2);                                                          % �������ϵ��k1

% ring
R = 8*1e-6;                                                                 % ���廷ǻ�뾶
L_rt = 2*pi*R;                                                              % �󻷳̳���
a = (0.92:0.03:0.98);                                                       % ����͸��ϵ��a

% response
Through = zeros(1,length(phase_rt));
alpha = zeros(1,length(a));

%��ͼ
figure(1);  
box on;

for ee = 1:length(alpha)  
    for ii = 1:length(phase_rt)
        Ht = (t1-a(ee)*exp(-1i*phase_rt(ii)))/(1-t1*a(ee)*exp(-1i*phase_rt(ii)));   % ���ݺ�����

        Tt = (abs(Ht))^2;                                                   % ���ʹ�����Ӧ
        Through(ii) = Tt;            
    end
    plot(phase_rt,Through);
    hold on;
end

xlabel('������λʧг��(����rt)');
set(gca,'XTick',(-0.3*pi:0.15*pi:0.3*pi));
set(gca,'XtickLabel',{'-0.3��','-0.15��','0','0.15��','0.3��'});
ylabel('͸����(Tap)');

title('APMR͸����');                                                         % ����ͼ������
legend({'overcoupling,t=0.92','critical coupling,t=a=0.95','undercoupling,t=0.98'},'location','southeast');


