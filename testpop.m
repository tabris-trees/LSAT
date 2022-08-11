%% 
% %===================================================
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is a Alfven wave-particle Heating pic code with matlab *
% 
% modified from a fortran program.
% 
% it won't generate the v_para data!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %===================================================

function testpop(setnumber,set,filenumber)
    %input settings
    tp0 = 0.5; % 质子的初始温度？
    mp = 1.0;
    e = 1.0;
    wp = 1.0;
    beta = 0.01;
    va = 10.0;
    vthp = sqrt(beta)*va; % vthp=sqrt(beta)*va
    b0 = 1.0;
    kappa = set(8);
    Nw = set(3);
    %% 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ω1=0.25Ωp，the waves extend from ω1 = 0.25Ωp to ωN = 0.33Ωp
    
    % w1 = set(2)*wp;
    % wN = w1+0.08;
    % wa = w1:(wN-w1)/(Nw-1):wN;
    % % wa = w1;
    % %% 
    % % $$\sum_k{B_k^2/B_0^2} = 0.013,\ \mathrm{and}\ (B_j/B_1)^2 = (\omega_j/\omega_1)^{-q},\ 
    % % q\ \mathrm{is\ chosen\ as\ 1.667}$$
    
    % q = set(4);
    % sum_bk2 = set(5)*b0;
    % eta = (wa./wa(1)).^(-q);
    % sum_eta = sum(eta);
    % bk12 = sum_bk2/sum_eta;
    % bk = sqrt(bk12*eta);
    % % bk = sqrt(0.09)*b0;

    w1 = set(2)*wp;
    wN = set(14)*wp;
    wa = linspace(w1,wN,Nw);
    % wa = w1;
    %% 
    % $$\sum_k{B_k^2/B_0^2} = 0.013,\ \mathrm{and}\ (B_j/B_1)^2 = (\omega_j/\omega_1)^{-q},\ 
    % q\ \mathrm{is\ chosen\ as\ 1.667}$$
    
    q = set(4);
    sum_bk2 = set(5)*b0;
    eta = (wa./wa(1)).^(-q);
    sum_eta = sum(eta);
    bk12 = sum_bk2/sum_eta;
    bk = sqrt(bk12*eta);
    % bk = sqrt(0.09)*b0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aplpa = set(6); % 波的传播方向
    %calculate plasma parameters

    % maxw = max(wa);
    % maxknumber = maxw/va;
    % wavelength = 2*pi/maxknumber;
    
    
    %set simulation domain
    % dx = wavelength/10;
    % l = 1080*va;
%     l = set(12)*va*((2*pi)/w1);
    l = set(12)*va;
    % dx = 20;
    dx = set(13);
%     dx = 24*pi; % 24Π*24Π的模拟尺度
    % dt = 0.025; % 时间取0.025*Ωp^-1
    dt = set(15);
    ts = set(7);
    ng = l/dx; % 网格个数
    % n = 10000; % 153600个质子
    n = set(11);
    disp(n);
    nt = ts/dt;
    
    
    %compute some other values
    qm = e/mp; % 荷质比，使用Boris方法push the particles的时候会用到（其实就是计算回旋频率的时候会用到吧）
    %% Despersion Relation
    % $$\omega = k_zv_A,\ v_A = B_0/(4\pi n_0 m_i)^{1/2}$$
    
    kz = wa/va;
    kx = tan(aplpa)*kz;
    % vaz = va;
    % vax = va*tan(aplpa);

    %%
    %initial
    % x=(1:n)*0;
    % y=x;
    % z=x;
    % vx=x;
    % vy=x;
    % vz=x;

    rng(1);
    % pphi_k = (0.+rand(1,Nw).*30)/360*(2*pi);
    % pphi_k = [0];
    pphi_k = linspace(0,set(16),Nw);
    if Nw == 1
        pphi_k = 0;
    end
    phi_k = zeros([n,Nw]);
    b_k = zeros([n,Nw]);
    wa_k = zeros([n,Nw]);
    kx_k = wa_k;
    kz_k = wa_k;
    for one = 1:n
        phi_k(one,:)=pphi_k;
        b_k(one,:)=bk;
        wa_k(one,:)=wa;
        kx_k(one,:)=kx;
        kz_k(one,:)=kz;
    end
    
%     phifile = fopen('phi.txt','a');
%     fprintf(phifile,'%6.4f\t',pphi_k./(2*pi));
%     fclose(phifile);
    
    E(1:n)=1;
    EE(1:n)=2;
%     Wa(1:n)=wa;

    % bb0(1:n)=b0;
    
    
    %load
    v0z = set(9)*va;
    v0 = [0,0,v0z];
    l0 = set(10)*l; % 粒子放置的初始位置（z方向总长度的分数倍）
    [r,vv] = dom(l0,n,v0,vthp,kappa,va);
    x = r(1,:); y = r(2,:); z = r(3,:);
    vx = vv(1,:); vy = vv(2,:); vz = vv(3,:);



    %%
    %cycle start
    tic;
    filesaver = ['setup-',num2str(filenumber),'/set',num2str(setnumber)];
    mkdir(filesaver);
    file2=fopen([filesaver,'/Txyz.dat'],'w');
    % filep = fopen([filesaver,'/poincare_data.dat'],'w');
%     file2=fopen([filesaver,'-Txyz.dat'],'w');
    for i=1:nt
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.Average the speed of the z-direction
        vzsum=sum(vz);
        vzl=vzsum/n;
        vzla=vzl/va;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2.Average the temperature
        %prepare the container
        sum_vx_pzi=(1:ng)*0;
        sum_vy_pzi=sum_vx_pzi;
        sum_vz_pzi=sum_vx_pzi;
        k_pzi=sum_vx_pzi;
        vxa=sum_vx_pzi;
        vya=sum_vx_pzi;
        vza=sum_vx_pzi;
        vxs_th_pp=sum_vx_pzi;
        vys_th_pp=sum_vx_pzi;
        vzs_th_pp=sum_vx_pzi;
        vxs_th_sumzi=sum_vx_pzi;
        vys_th_sumzi=sum_vx_pzi;
        vzs_th_sumzi=sum_vx_pzi;
        vxys_th_azi=sum_vx_pzi;
        vzs_th_azi=sum_vx_pzi;
        %count and compute the avrange
        pep=z/dx; 
        ziep=ceil(pep); % 位于第几个网格
        inform=[ziep;vx;vy;vz];
        for p=1:n
            zi=inform(1,p);
            sum_vx_pzi(zi)=sum_vx_pzi(zi)+inform(2,p);
            sum_vy_pzi(zi)=sum_vy_pzi(zi)+inform(3,p);
            sum_vz_pzi(zi)=sum_vz_pzi(zi)+inform(4,p);
            k_pzi(zi)=k_pzi(zi)+1;
        end
        
        for iz=1:ng
            if k_pzi(iz)==0
                k_pzi(iz)=1;
            end
            vxa(iz)=sum_vx_pzi(iz)/k_pzi(iz);
            vya(iz)=sum_vy_pzi(iz)/k_pzi(iz);
            vza(iz)=sum_vz_pzi(iz)/k_pzi(iz);
        end
        
        for p=1:n
            zi=inform(1,p);
            vxs_th_pp(zi)=(inform(2,p)-vxa(zi))^2;
            vxs_th_sumzi(zi)=vxs_th_sumzi(zi)+vxs_th_pp(zi);
            vys_th_pp(zi)=(inform(3,p)-vya(zi))^2;
            vys_th_sumzi(zi)=vys_th_sumzi(zi)+vys_th_pp(zi);
            vzs_th_pp(zi)=(inform(4,p)-vza(zi))^2;
            vzs_th_sumzi(zi)=vzs_th_sumzi(zi)+vzs_th_pp(zi);
        end
        
        for iz=1:ng
            vxys_th_azi(iz)=(vxs_th_sumzi(iz)+vys_th_sumzi(iz))/k_pzi(iz);
            vzs_th_azi(iz)=vzs_th_sumzi(iz)/k_pzi(iz);
        end
        
        vxys_th_a=sum(vxys_th_azi)/ng;
        vzs_th_a=sum(vzs_th_azi)/ng;
        
        txy=(mp/2)*vxys_th_a;
        tz=mp*vzs_th_a;
        tk=(vxys_th_a+vzs_th_a)/(va^2);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3.Out the data
        if i==1 || mod(i*100,nt)==0 % || i*dt==40 || i*dt==600 || i*dt==800 || i*dt==1000 || i*dt==2000 || i*dt==5000
        % if i==1 || i*dt==50 || i*dt==100 || i*dt==500 || i*dt==1000 || i*dt==5000 || i*dt==20000
            time=i*dt;
            if time<1
                time=0;
            end
            ni=num2str(time);
            contin=[filesaver,'/zVxyz',ni,'.dat'];
%             contin2=[filesaver,'/poincare_data',ni,'.dat'];
            file1=fopen(contin,'w');
%             filep=fopen(contin2,'w');
            for p=1:n
                fprintf(file1,'%.5f\t%.5f\t%.5f\t%.5f\n',z(p),vx(p)/va,vy(p)/va,vz(p)/va);
            end
            fclose(file1);
        end
        
        % fprintf(filex, '%.5f\t%.5f\t%.5f\n', vx(1),vy(1),vz(1));
        
        fprintf(file2,'%.5f\t%.5f\t%.5f\t%.5f\n',txy,tz,tk,vzla);    
        
        
        if mod(i*10,nt)==0
            disp(strcat(num2str(i*100/nt),'% of set',num2str(setnumber),'of',num2str(filenumber)));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4.Push the particle,(compute the B-field and velocity of present time)
        % 这一部分还有很大的问题，暂时先按照源程序直接搬过来的：
        
        %prepare the container
    %% 
    % 将所有需要参与磁场计算的参数全部变成n*Nw的二维数组
		x_k=zeros([n,Nw]);
	    z_k=zeros([n,Nw]);
		for k = 1:Nw
			x_k(:,k)=x;
			z_k(:,k)=z;
		end
        
        %compute

    %% Magnetic Field in wave frame
    % $$$$\mathbf{B}_{w}=\sum_{k=1}^{N} B_{k}\left[-\cos (\alpha) \sin \left(\psi_{k}\right) 
    % \mathbf{i}_{x}+\cos \left(\psi_{k}\right) \mathbf{i}_{y}+\sin (\alpha) \sin 
    % \left(\psi_{k}\right) \mathbf{i}_{z}\right]$$$$
    % 
    % $\psi_{k}=k_x x+k_z z+\phi_k \ , \tan(\alpha)=k_x/k_z $，
    % 
    % $$\phi_k$$第k个子波的随机相位差，$N$是波模的个数。
    % 
    % 粒子的运动中不考虑电场的影响（原文应该是在波的参考系中电场是被抵消掉的），运动的整合方式采用Boris算法（具体可以参考particle in cell 
    % simulation的网页）
    % 
    % 
    
    % wave frame
%         psi_k = kz_k.*z_k+phi_k; % 这里稍微有点问题
% 
%         bx=sum(bk.*cos(psi_k),2)';
%         by=-sum(bk.*sin(psi_k),2)';
%         bz(1:n)=b0; 
%         paralle-left-rotate wave


    % lab frame

        psi_k = wa_k*i*dt-(kz_k.*z')+phi_k;
        bx=sum(bk.*cos(psi_k),2)';
        by=-sum(bk.*sin(psi_k),2)';
        bz(1:n)=b0;
        ex=va*by;
        ey=-va*bx;
        ez(1:n)=0;
        


    %% Boris Method to push the particle in magnetic field
    % $$$$v'=v^-+v^-\times t \\\mathrm{t} \equiv-\hat{\mathrm{b}} \tan \frac{\theta}{2}=\frac{q 
    % \mathbf{B}}{m} \frac{\Delta t}{2} \\v^+=v^-+v'\times s \\s=\frac{2t}{1+t^2} 
    % \\\rm{t}\ and\ \rm{s}\ are\ vectors$$$$
   
    

%         tmx=dt*bx*qm/2;
%         tmy=dt*by*qm/2;
%         tmz=dt*bz*qm/2;
%         tm2=tmx.^2+tmy.^2+tmz.^2;
%     
%         sx=(EE.*tmx)./(E+tm2);
%         sy=(EE.*tmy)./(E+tm2);
%         sz=(EE.*tmz)./(E+tm2);
%         
%         vpx=vx+(vy.*tmz-vz.*tmy);
%         vpy=vy+(vz.*tmx-vx.*tmz);
%         vpz=vz+(vx.*tmy-vy.*tmx);
%     
%         vx=vx+(vpy.*sz-vpz.*sy);
%         vy=vy+(vpz.*sx-vpx.*sz);
%         vz=vz+(vpx.*sy-vpy.*sx);
        
        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        % lab-frame

        tm=dt*qm/2;
        s=EE./(E+(bx.^2+by.^2+bz.^2)*tm^2);

        vlasty = vy;

        vxa=vx+tm*ex;
        vya=vy+tm*ey;
        vza=vz+tm*ez;
    
        vxc=vxa+tm*(vya.*bz-vza.*by);
        vyc=vya+tm*(vza.*bx-vxa.*bz);
        vzc=vza+tm*(vxa.*by-vya.*bx);

        vxb=vxa+s.*(tm*(vyc.*bz-vzc.*by));
        vyb=vya+s.*(tm*(vzc.*bx-vxc.*bz));
        vzb=vza+s.*(tm*(vxc.*by-vyc.*bx));

        vx=vxb+tm*ex;
        vy=vyb+tm*ey;
        vz=vzb+tm*ez;
% %         

        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5.Move the particle
        for p=1:n
            z(p)=z(p)+vz(p)*dt;
            if z(p)>l
                z(p)=z(p)-l;
            elseif z(p)<=0
                z(p)=z(p)+l;
            end
                
            % if vy(p) > 0 && vlasty(p) < 0
            %     fprintf(filep,'%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n',i*dt,z(p),vz(p),vx(p),vy(p));
            % end            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Program end
    end
    fclose(file2);
    fclose(filep);
    toc;
    disp(['运行时间: ',num2str(toc)]);
end
%%
function [r,vv]=dom(l,n,v0,vth,kappa,va)
    L(1:n)=l;
%     disp(L);
    s=rand(3,n);
%     disp(s);
    r = zeros(3,n);
    for ir = 1:3
        r(ir,:)=L.*s(ir,:);
    end
    if kappa == Inf
        disp('Maxwellian distribution!')
        vvx = normrnd(v0(1), vth, 1, n);
        vvy = normrnd(v0(2), vth, 1, n);
        vvz = normrnd(v0(3), vth, 1, n);
        vv = [vvx;vvy;vvz];
    elseif kappa == 0
        disp('Uniform distribution!')
        vvx = zeros(1,n);
        vvy = zeros(1,n);
        vvz = linspace(-1,2,n)*va;
        vv = [vvx;vvy;vvz];
    else
        disp('Kappa distribution!')
        vv=rand_kappa3(vth,kappa,n,3,v0);
    end
end
