clear
clc
close all
%初始化
h=2*pi/700;
x0=0:h:2e5*pi/700;
y0=[0.1;0.1];%最后一项是cos(w*t)，当t=0时必须为1
%这里一口气求出一大堆点不太现实，内存有限，所以采取分批计算，逐步绘图来实现
figure(1)
hold on
for k=1:100
    [y1,Output]=ODE_RK4_hyh(x0,h,y0,[0.1,0.35,1.4]);
    Lx=y1(1,:);
    Ly=y1(2,:);
    %采用频闪采样法计算
    %由于w变为了1.4，这里还是z=0平面，采样间隔改为 (3*pi/2:(2*pi/1):2000)/1.4
    yP_List=y1(:,376:500:end);%注，这里376=3/2/1.4*350+1

    %内存有点不够了，分段进行下面的计算
    x0=x0+2e5*pi/700;%时间整体向后移动
    y0=y1(:,end);%把最后一刻的状态当做下一个时间段初始状态
    %绘图
    %1庞加莱截面
    scatter(yP_List(1,1:end),yP_List(2,1:end),0.8,'k','MarkerEdgeAlpha',0.6)
    %disp(k)
end

function [F,Output]=Fdydx(x,y,Input)
%形式为Y'=F(x,Y)的方程，参见数值分析求解常系数微分方程相关知识
%高次用列向量表示，F=[dy(1);dy(2)];y(1)为函数，y(2)为函数导数
d=Input(1);
r=Input(2);
w=Input(3);
dy(1)=y(2);
dy(2)=-y(1)^3+y(1)-d*y(2)+r*cos(w*x);
F=[dy(1);dy(2)];
Output=[];
end

function [y,Output]=ODE_RK4_hyh(x,h,y0,Input)
%4阶RK方法
%h间隔为常数的算法
y=zeros(size(y0,1),size(x,2));
y(:,1)=y0;
for ii=1:length(x)-1
    yn=y(:,ii);
    xn=x(ii);
    [K1,~]=Fdydx(xn    ,yn       ,Input);
    [K2,~]=Fdydx(xn+h/2,yn+h/2*K1,Input);
    [K3,~]=Fdydx(xn+h/2,yn+h/2*K2,Input);
    [K4,~]=Fdydx(xn+h  ,yn+h*K3  ,Input);
    y(:,ii+1)=yn+h/6*(K1+2*K2+2*K3+K4);
end
Output=[];
end
