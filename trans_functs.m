clear all;
close all;

phase_rt = (-pi:pi/100:3*pi);                                               % 遍历环程相移

% coupler 1 
k1 = sqrt(0.2);                                                             % 定义耦合系数k1
t1 = sqrt(1-(k1)^2);                                                        % 耦合无损耗情况下，定义传输系数t1
M1 = [t1,-1i*k1; -1i*k1,t1];                                                % 传输矩阵

% ring
R = 8*1e-6;                                                                 % 定义环腔半径
L_rt = 2*pi*R;                                                              % 求环程长度
a = 0.98;                                                                   % 环程透过系数a

% coupler 2
k2 = sqrt(0.2);                                                             % 定义耦合系数k2
t2 = sqrt(1-(k2)^2);                                                        % 耦合无损耗情况下，定义传输系数t2
M2 = [t2,-1i*k2; -1i*k2,t2];                                                % 传输矩阵

% coupler 2 for critical coulping
% t2 = t1/a;                                                                  % 定义耦合系数k2,临界耦合条件下
% k2 = sqrt(1-(k1)^2);                                                        % 耦合无损耗情况下，定义传输系数t2
% M2 = [t2 -1i*k2; -1i*k2 t2];                                                % 传输矩阵

% response
Through = zeros(1,length(phase_rt));                                        % 直通端矩阵
Drop = zeros(1,length(phase_rt));                                           % 下载端矩阵


for ii = 1:length(phase_rt)
    P = [0,sqrt(a)*exp(-1i*phase_rt(ii)/2); 1/(sqrt(a)*exp(-1i*phase_rt(ii)/2)),0]; % 微环的传输矩阵
    H = M1*P*M2;
    
    Ht = (t1-t2*a*exp(-1i*phase_rt(ii)))/(1-t1*t2*a*exp(-1i*phase_rt(ii))); %传递函数法
    Hd = -(k1*k2*sqrt(a)*exp(-1i*phase_rt(ii)/2))/(1-t1*t2*a*exp(-1i*phase_rt(ii)));    

    Tt = (abs(Ht))^2;                                                       % 功率光谱响应
    Td = (abs(Hd))^2;
                                                                            % 总传输矩阵
    Through(ii) = Tt;
    Drop(ii) = Td;
end

figure(1);                                                                  % 绘图
box on;

plot(phase_rt,10*log10(Drop),'r');

hold on;

plot(phase_rt,10*log10(Through),'b');

set(gca,'XTick',(-pi:pi:3*pi));                                              % 设置x坐标
set(gca,'XtickLabel',{'-π','0','π','2π','3π'});
xlabel('环程相位失谐量(Δφrt)');

set(gca,'YTick',[-30,-25,-20,-15,-10,-5,0,]);                               % 设置y坐标
set(gca,'YtickLabel',{'-30','-25','-20','-15','-10','-5','0'});
ylabel('传输率(dB)');

title('ADMR的功率光谱响应,κ1=κ2,a=0.98');                                                 % 设置图表名称
legend('Drop Port','Through Port');                                         % 设置图例

