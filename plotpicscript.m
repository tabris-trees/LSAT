for ii = 1:1

% 动画
    figure;
    sets = ['setup-7/set',num2str(ii)];
    pic_num = 1;
    l = 20000;
    time_ftp = 200;
    for i = 0:time_ftp:time_ftp*100 
        phaseplot(sets,num2str(i),l); % 画图函数
        drawnow; % 画gif动图的时候这个外面的大for循环就存在，直接把以下和for*******以上的代码放到对应的位置即可
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        gifname = strcat(sets,' phase space plot.gif');
        if pic_num == 1
            imwrite(I,map,gifname,'gif','Loopcount',inf,'DelayTime',0.2);
        else
            imwrite(I,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
        end
        pic_num = i/time_ftp + 1;
    end

% 最终
%     figure;
%     sets = ['set',num2str(ii)];
%     phaseplot(sets,num2str(500));
%     picname2 = [sets,' phase-500'];
%     saveas(gcf,picname2,'fig');
%     saveas(gcf,picname2,'png');

% 各向异性：
%     figure;
%     sets = ['set',num2str(ii)];
%     taniso(sets,25000);
%     figname1 = ['aniso_',sets];
%     saveas(gcf, figname1, 'fig');
%     saveas(gcf, figname1, 'png');
end

%######################################################################################################
%% 基础制图
clear;
warning('off');
% 验证有电场版本是否正确
% figure;
% templot('E_set1\Txyz.dat',1,2500);
linemarks = {'-','--',':','-.'};
figure(998);
clf;
% for iii = 1:2
% subplot(1,2,iii)
names = [];
for i = 8:8
    % figure(2);
    % clf;
    % names = [];
    % linemark = linemarks{i};
    linemark = '-';
    header = ["set","w1","Nw","q","bk^2","aplpa",...
    "time","kappa","vz0/va","l0","n","l","dx","dt"];
    for ii = 1:5
%         subplot(2, 2, ii);
        name=vparaplot(['setup-',num2str(i),'\set',num2str(ii),'\Txyz.dat'],ii,20000,[3,5],i,header);
        names = [names,name];
%         title(name);
    end
    % legend(names,'Location','northeast');
%     legend(names,'Location','southeast');
    % print(gcf,'V_para-kappa3.jpg','-r600','-djpeg');
end
% line([2160,2160], [0,1.2], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
h = legend(names,'Location','northeast');
title('the U_{||} evolution with time');
% end
print(gcf,'setup-8/upara.jpg','-r600','-djpeg');

%######################################################################################################

%% poincare map
setupindex = '5';
setindex = '1';
pd = ['setup-',setupindex,'\set',setindex,'\poincare_data.dat'];
data2=importdata(pd);
tttimes = data2(:,1);
pic_num = 1;
figure(333);
clf;
set(gcf,'Units','centimeters','Position',[2 2 14 14]);
x1=[];y1=[];
% for flag = 0.025:10:1000
for flag = 20000:20000
    points = tttimes(tttimes<=flag);
    cutindex = length(points);
    x1 = data2(1:cutindex,2);
    y1 = data2(1:cutindex,3);

    textstr = num2str(flag);
    x=mod(0.003*x1,2*pi);
%     x = 0.003*x1;
    y=y1/10;

    plot(x,y,'.k','MarkerSize',1);
    % hold on;
    xlim([0,2*pi]);
%     ylim([-1,2]);
    text(1,-0.4,textstr);
%     drawnow;
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     gifname = ['setup-',setupindex,'/set',setindex,'poincare.gif'];
%     if pic_num == 1
%         imwrite(I,map,gifname,'gif','Loopcount',inf,'DelayTime',0.2);
%         pic_num = 2;
%     else
%         imwrite(I,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
%     end
end
% print(gcf,['setup-',setupindex,'/set',setindex,'poincare.png'],'-r600','-dpng');
%#####################################################################################################

%% the velocitys space
setupindex = '7';
setindex = '1';
figure(76);
set(gcf,'Units','centimeters','Position',[2 2 40 14]);
clf;
sets = ['setup-',setupindex,'/set',setindex];
plottimes = [0,200,400,800,1600,3200,6400,12800,16000,20000];
pic_num = 1;
% for plottime = linspace(0,20000,101)
%     velocityspace(sets,num2str(plottime),0.5);

for plottime = 1:length(plottimes)
    subplot(2,5,plottime);
    velocityspace(sets,num2str(plottimes(plottime)),0.5);

%     drawnow;
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     gifname = ['setup-',setupindex,'/pitch.gif'];
%     if pic_num == 1
%         imwrite(I,map,gifname,'gif','Loopcount',inf,'DelayTime',0.2);
%         pic_num = 2;
%     else
%         imwrite(I,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
%     end
end
print(gcf,['setup-',setupindex,'/set',setindex,'pitchangle.png'],'-r600','-dpng');
