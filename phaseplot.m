function phaseplot(sets, filenames)
    filename=strcat(sets,'\zVxyz',filenames,'.dat');
    disp(filename);
    labeltext=strcat('\Omega_{p}t=',filenames);
    data=importdata(filename);
    z=data(:,1);
    vx=data(:,2);
    vy=data(:,3);
    vz=data(:,4);
    uxy=vx.^2+vy.^2;
    
    set(gcf,'Units','centimeters','Position',[10 10 21 7]);

    subplot(1,2,1);
    scatter(z,vx,1,"black","filled");
    ylabel('V_{x}/V_{A}');
    xlabel('Z\Omega_{p}/V_{A}');
    text(100,1.6,labeltext); 
    ylim([-2.5,2.5]);
%     ylim([-2,2]);
    
%         xticks(1600:800:5600);
%         yticks(-2:1:2);
%         set(gcf,'Units','centimeters','Position',[10 10 21 29]);
    
    subplot(1,2,2);
    scatter(z,vz,1,"red","filled");
    ylabel('V_{z}/V_{A}');
    xlabel('Z\Omega_{p}/V_{A}');
    text(100,1.6,labeltext);
    % xlim([0,1000]);
    ylim([-1,3]);
%     axis([3000,4500,-2,2]);
    
    %     xticks(1600:800:5600);
    %     yticks(-0.5:0.5:1.5);
    %     set(gcf,'Units','centimeters','Position',[10 10 17 25]);
    
%         picname=strcat('phase space-',subfilename);
%         saveas(gcf, picname, 'png');
end