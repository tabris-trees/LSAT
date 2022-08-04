function name=taniso(txyz,sets,xmax,posi,setup)
    name0 = '';
    setfile = importdata(['setup\setup',num2str(setup),'.txt']);
    setdata = setfile.data(sets,:);
    header = ["w1","Nw","q","\sum b_k^2","aplpa","time","\kappa"];
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
    
%     subplot(1,2,1);
    plot(time,txy./tz);
    hold on;
    % xticks(tick);
%     yticks(0:2:16);
    % xticklabels(lable);
    xlim([0,xmax]);
    ylim([0,3]);
    xlabel('\Omega_{p}t');
    ylabel('$T_{\perp}/T_{\parallel}$',"Interpreter","latex")
end