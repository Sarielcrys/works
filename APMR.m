% Power spectral response of APMR
% Source code by Nuo Chen
% 2020/10/02

clear all;
close all;

phase_rt = (-0.3*pi:pi/200:0.3*pi);                                         % 遍历环程相移

% coupler 1 
t1 = 0.95;                                                                  % 耦合无损耗情况下，定义传输系数t1
k1 = sqrt(1-t1^2);                                                          % 定义耦合系数k1

% ring
R = 8*1e-6;                                                                 % 定义环腔半径
L_rt = 2*pi*R;                                                              % 求环程长度
a = (0.92:0.03:0.98);                                                       % 环程透过系数a

% response
Through = zeros(1,length(phase_rt));
alpha = zeros(1,length(a));

%绘图
figure(1);  
box on;

for ee = 1:length(alpha)  
    for ii = 1:length(phase_rt)
        Ht = (t1-a(ee)*exp(-1i*phase_rt(ii)))/(1-t1*a(ee)*exp(-1i*phase_rt(ii)));   % 传递函数法

        Tt = (abs(Ht))^2;                                                   % 功率光谱响应
        Through(ii) = Tt;            
    end
    plot(phase_rt,Through);
    hold on;
end

xlabel('环程相位失谐量(Δφrt)');
set(gca,'XTick',(-0.3*pi:0.15*pi:0.3*pi));
set(gca,'XtickLabel',{'-0.3π','-0.15π','0','0.15π','0.3π'});
ylabel('透射谱(Tap)');

title('APMR透射谱');                                                         % 设置图表名称
legend({'overcoupling,t=0.92','critical coupling,t=a=0.95','undercoupling,t=0.98'},'location','southeast');


