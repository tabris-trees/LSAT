function velocityspace(sets, filenames, prefile)
    filename=strcat(sets,'\zVxyz',filenames,'.dat');
    prefilename=strcat(sets,'\zVxyz',prefile,'.dat');
    disp(filename);
    labeltext=strcat('\Omega_{p}t=',filenames);
    data=importdata(filename);
    predata=importdata(prefilename);
    % z=data(:,1);
    vx=data(:,2);
    vy=data(:,3);
    vz=data(:,4);
    prevx=predata(:,2);
    prevy=predata(:,3);
    prevz=predata(:,4);
    uxy=sqrt(vx.^2+vy.^2);
    pre_energy = prevx.^2+prevy.^2+prevz.^2;
    % energy = vx.^2+vy.^2+vz.^2;
    
%     set(gcf,'Units','centimeters','Position',[2 2 18 18]);

%     plot3(vx,vy,vz,'.','Color','k');
%     xlim([-1,1]);
%     ylim([-1,1]);
%     zlim([-0.5,2]);

    scatter(vz,uxy,1,pre_energy);
%     hold on;
%     colormap 'parula';
    colormap(othercolor('StepSeq_25'));
    % caxis([0,2]);
    ylim([0,2]);
    xlim([-0.5,2.5]);
    title(labeltext);

end