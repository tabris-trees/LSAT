% for ii = 1:3

% % 动画
%     figure;
%     sets = ['setup-2/set',num2str(ii)];
%     pic_num = 1;
%     time_ftp = 50;
%     for i = 0:time_ftp:time_ftp*100 
%         phaseplot(sets,num2str(i)); % 画图函数
%         drawnow; % 画gif动图的时候这个外面的大for循环就存在，直接把以下和for*******以上的代码放到对应的位置即可
%         F=getframe(gcf);
%         I=frame2im(F);
%         [I,map]=rgb2ind(I,256);
%         gifname = strcat(sets,' phase space plot.gif');
%         if pic_num == 1
%             imwrite(I,map,gifname,'gif','Loopcount',inf,'DelayTime',0.2);
%         else
%             imwrite(I,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
%         end
%         pic_num = i/time_ftp + 1;
%     end

% % 最终
% %     figure;
% %     sets = ['set',num2str(ii)];
% %     phaseplot(sets,num2str(500));
% %     picname2 = [sets,' phase-500'];
% %     saveas(gcf,picname2,'fig');
% %     saveas(gcf,picname2,'png');

% % 各向异性：
% %     figure;
% %     sets = ['set',num2str(ii)];
% %     taniso(sets,25000);
% %     figname1 = ['aniso_',sets];
% %     saveas(gcf, figname1, 'fig');
% %     saveas(gcf, figname1, 'png');
% end

%######################################################################################################

% tts = [0,50,100,500,1000,5000,20000];

% % for tt = 1:length(tts)
% for tt = 1:length(tts)
%     figure(tt);
%     phaseplot('setup-4/set4',num2str(tts(tt)));
% end


%######################################################################################################
%% 基础制图
clear;
warning('off');
% 验证有电场版本是否正确
% figure;
% templot('E_set1\Txyz.dat',1,2500);
linemarks = {'-','--',':','-.'};
figure(999);
clf;
% for iii = 1:2
% subplot(1,2,iii)
names = [];
for i = 4:4
    % figure(2);
    % clf;
    % names = [];
    % linemark = linemarks{i};
    linemark = '-';
    for ii = 1:4
        % subplot(2, 2, ii);
        name=vparaplot(['setup-',num2str(i),'\set',num2str(ii),'\Txyz.dat'],ii,5000,[2],i);
        names = [names,name];
        % title(name);
    end
    % legend(names,'Location','northeast');
    % legend(names,'Location','southeast');
    % print(gcf,'V_para-kappa3.jpg','-r600','-djpeg');
end
h = legend(names,'Location','southeast');
title('the U_{||} evolution with time');
% end
print(gcf,'setup-4/upara.jpg','-r600','-djpeg');

%######################################################################################################

%% poincare map
pd = ['setup-',num2str(3),'\set',num2str(1),'\poincare_data.dat'];
data2=importdata(pd);
tttimes = data2(:,1);
pic_num = 1;
figure(333);
x1=[];y1=[];
for flag = 0.025:10:20000
    points = tttimes(tttimes<=flag);
    cutindex = length(points);
    x1 = data2(1:cutindex,2);
    y1 = data2(1:cutindex,3);

    textstr = num2str(flag);
    x=mod(x1,pi);
    y=y1/10;

    scatter(x,y,3,'k',"filled");
    % hold on;
    % xlim([0,2*pi]);
    ylim([-1,2]);
    text(1,1.8,textstr);
    drawnow;
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    gifname = ['setup-',num2str(3),'\perdopoincare.gif'];
    if pic_num == 1
        imwrite(I,map,gifname,'gif','Loopcount',inf,'DelayTime',0.2);
        pic_num = 2;
    else
        imwrite(I,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
    end

end
%#####################################################################################################

%% the velocitys space
clear;
setupindex = '4';
setindex = '1';
figure(str2double(setindex));
set(gcf,'Units','centimeters','Position',[2 2 40 14]);
clf;
sets = ['setup-',setupindex,'/set',setindex];

% plottimes = [0,20,40,80,160,320,640,1280,1600,2000];
% for plottime = 1:length(plottimes)-1
%     subplot(2,6,plottime);
%     velocityspace(sets,num2str(plottimes(plottime)),num2str(plottimes(plottime)));
%     subplot(2,6,plottime+6);
%     velocityspace(sets,num2str(plottimes(plottime+1)),num2str(plottimes(plottime)));

pic_num = 1;
for plottime = linspace(0,5000,101)
    velocityspace(sets,num2str(plottime),num2str(plottime));
    drawnow;
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    gifname = ['setup-',setupindex,'/pitch-4.gif'];
    if pic_num == 1
        imwrite(I,map,gifname,'gif','Loopcount',inf,'DelayTime',0.2);
        pic_num = 2;
    else
        imwrite(I,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
    end
end
% print(gcf,['setup-',setupindex,'/pitch','-',setindex,'.png'],'-r600','-dpng');