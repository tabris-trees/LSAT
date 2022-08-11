function velocityspace(sets, filenames,pointalpha)
    filename=strcat(sets,'\zVxyz',filenames,'.dat');
    disp(filename);
    labeltext=strcat('\Omega_{p}t=',filenames);
    data=importdata(filename);
    z=data(:,1);
    vx=data(:,2);
    vy=data(:,3);
    vz=data(:,4);
    uxy=sqrt(vx.^2+vy.^2);
    energy = vx.^2+vy.^2+vz.^2;
    
%     set(gcf,'Units','centimeters','Position',[2 2 18 18]);

    % scatter3(vx,vy,vz,3,'k','filled','MarkerFaceAlpha',pointalpha);
    % xlim([-1,1]);
    % ylim([-1,1]);
    % zlim([-0.5,2]);

    scatter(vz,uxy,1,energy);
    scatter(vz,uxy,3,'black','filled','MarkerFaceAlpha' ...
        ,pointalpha,'MarkerEdgeAlpha',pointalpha);

%     hold on;
%     colormap 'parula';
%     colormap(othercolor('Bu_10')); 
    caxis([0,1]);
    ylim([0,2]);
    xlim([-0.5,2.5]);
    title(labeltext);

end