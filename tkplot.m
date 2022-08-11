function name = tkplot(txyz,sets,xmax,posi,setup,header)
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

    txy=data2(:,1);
    tz=data2(:,2);
    tk=data2(:,3);
    vzl=data2(:,4);
    time=(0:length(txy)-1)*0.025;
    tk2 = tk-tk(1);
    
    set(gcf,'Units','centimeters','Position',[10 10 11 9]);

%% kinetic temperature    
    % plot(time,tk);
    % hold on;
    % xlabel('\Omega_{p}t');
    % ylabel('T_{k}/(1/2m_{p}V_{A})');

    % minus the initial kinetic temperature
    % plot(time,tk2,"LineWidth",1);
    plot(time,tk2);
    hold on;
    xlim([0,xmax]);
    xlabel('\Omega_{p}t');
    ylabel('(T_{k}-T_{k0})/(1/2m_{p}V_{A})');
    
end
