function name = vparaplot(txyz,sets,xmax,posi,setup,header)
    name0 = '';
    setfile = importdata(['setup\setup',num2str(setup),'.txt']);
    setdata = setfile.data(sets,:);
    for x = 1:length(posi)
        index = posi(x);
        setsvalue = num2str(setdata(index));
        if setsvalue == "Inf"
            setsvalue = "\infty";
        end
        subname = strcat(header(index),"=",setsvalue);
        name0 = strcat(name0,subname,';');
    end

    name = name0;
    
    data2=importdata(txyz);

    txy=data2(:,1);
    tz=data2(:,2);
    tk=data2(:,3);
    vzl=data2(:,4);
    time=(0:length(txy)-1)*0.025;
    tk2 = tk-tk(1);
    
    set(gcf,'Units','centimeters','Position',[10 10 11 9]);

    plot(time,vzl);
    hold on;
    xlim([0,xmax]);
    % ylim([-1,0]);
    % yticks(0:0.2:1.2);
    % xticklabels(lable);
    xlabel('\Omega_{p}t');
    ylabel('$U_{\parallel}/V_{A}$','Interpreter','latex');
    text2 = strcat("$",name,"$");
    % text(xmax*0.7,vzl((xmax*40)*0.7),text2,"Interpreter","latex",'Color','black');
    

end
