function name=templot(txyz,sets,xmax,posi,setup,header)
    name0 = '';
    setfile = importdata(['setup\setup',num2str(setup),'.txt']);
    setdata = setfile.data(sets,:);
    for x = 1:length(posi)
        index = posi(x);
        setsvalue = num2str(setdata(index+1));
        if setsvalue == "Inf"
            setsvalue = "\infty";
        end
        subname = strcat(header(index),"=",setsvalue);
        name0 = strcat(name0,subname,';');
    end

    name = name0;
    
    data2=importdata(txyz);
%     disp(data2);
    txy=data2(:,1);
    tz=data2(:,2);
    tk=data2(:,3);
    vzl=data2(:,4);
    time=(0:length(txy)-1)*0.025;
    tk2 = tk-tk(1);
    
%     lable=0:500:2500;
%     tick=0:20000:100000;
    
%     figure;
    set(gcf,'Units','centimeters','Position',[10 10 11 9]);
    
    p1=plot(time,tz,time,txy);
    % p1c = get(p1,'Color');
    hold on;
    % plot(time,txy,'Color',p1c);
%     xticks(tick);
    % yticks(0:5:50);
%     xticklabels(lable);
    xlim([0,xmax]);
%     ylim([0,50]);
    xlabel('\Omega_{p}t');
    ylabel('T_{||} and T_{\perp}');
    perptext=strcat("T_{\perp}");
    paratext=strcat("T_{||}");
    % text(xmax*0.7,txy((xmax*40)*0.7),perptext,"Interpreter","latex",'Color',p1c);
    % text(xmax ...
    % *0.7,tz(xmax*40*0.7),paratext,"Interpreter","latex","Color",p1c);
%     legend(paratext, perptext, 'Location', 'southeast');

end