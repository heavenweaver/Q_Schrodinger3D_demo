%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vb=0.55;                %% Potential barrier height[eV]  1.424 eV - 0.324 eV = 1.1
Vb=1.1;                 %% Potential barrier height[eV]  1.424 eV - 0.324 eV = 1.1
Mass = 0.043;           %% effective mass, constant over all the structure...  mΓ = 0.023mo
Mx=12e-9;               %% map X [m]
My=12e-9;               %% map Y [m]
Mz= 6e-9;               %% map Z [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WLt = 0.5e-9;           %% Wetting Layer thickness [m]
QDh = 3e-9;             %% Quantum dot height [m]
QDd = 12e-9;            %% Quantum dot diameter [m]
%E  = 10e8;              %% Internal Electrical Field (10MV/cm in GaN/AlN interface)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);
z=linspace(-Mz,Mz,Nz)+0.5e-9;

[X,Y,Z]=meshgrid(x,y,z);


a=QDd/2;                      %% side length of the hexagon [m]
h=tan(pi/4) * a*sqrt(2)/2;    %% height of the pyramid with 45deg angle bw side face and base
%h=a/2*sqrt(2)/2;

% Pyramidal base
%idx1=  abs(X)<a*sqrt(2)/2;
%idx2=(tan(pi/4)*X+a>Y) .* (tan(pi/4)*X-a<Y) .* (-tan(pi/4)*X-a<Y) .* (-tan(pi/4)*X+a>Y);
idx1=  abs(X)<a;
idx2=(X+a>Y) .* (X-a<Y) .* (-X-a<Y) .* (-X+a>Y);

idx_QD=idx1.*idx2;
idx_QD(Z < 0) = 0;
idx_QD(Z > QDh) = 0;

% Pyramid facet
idx_QD(Z > -2*h/a/sqrt(2)*X + h) = 0;
idx_QD(Z > +2*h/a/sqrt(2)*X + h) = 0;
idx_QD(Z > -2*h/a/sqrt(2)*Y + h) = 0;
idx_QD(Z > +2*h/a/sqrt(2)*Y + h) = 0;

%idx_QD(Z > +h/a/sqrt(2)*X -h/a*Y + h) = 0;
%idx_QD(Z > -h/a/sqrt(2)*X -h/a*Y + h) = 0;
%idx_QD(Z > +h/a/sqrt(2)*X +h/a*Y + h) = 0;
%idx_QD(Z > -h/a/sqrt(2)*X +h/a*Y + h) = 0;

%idx_QD(Z > +h/a/sqrt(2)*X -h/a*Y + h) = 0;
%idx_QD(Z > -h/a/sqrt(2)*X -h/a*Y + h) = 0;
%idx_QD(Z > +h/a/sqrt(2)*X +h/a*Y + h) = 0;
%idx_QD(Z > -h/a/sqrt(2)*X +h/a*Y + h) = 0;

% Wetting Layer index
idx_WL = (Z < 0) & (Z > -WLt);

% Potential index
idx = idx_QD | idx_WL;
V0=(idx)*0 + (1-idx)*Vb ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0 = Fx*X+V0;        % adding the electric field Fx to the potential in the x-direction
V0 = Fy*Y+V0;        % adding the electric field Fy to the potential in the y-direction
V0 = Fz*Z+V0;        % adding the electric field Fz to the potential in the z-direction

V0=V0-min(min(min(V0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Potential','position',[100 100 1600 800])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1,'fontsize',15)
hold on;grid on;

pcolor(x*1e9,z*1e9,squeeze(V0(round(end/2),:,:))')
plot([-1 1]*10,[0 0],'w--','linewidth',2)

idz=find(z>QDh);
idz=idz(1)-1;
plot([-1 1]*10,[1 1]*z(idz)*1e9,'w--','linewidth',2)
plot([1 1]*0,[-1 1]*10,'b--','linewidth',2)
plot([1 1]*x(end)*1e9,[-1 1]*10,'r--','linewidth',2)

colormap(jet)
colorbar
shading flat
caxis([floor(min(V0(:)))  ceil(max(V0(:)))])

xlim([-1 1]*Mx/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('z (nm)')
title(strcat('Potential (eV): Vxz @y=0nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,4,'fontsize',15)
hold on;grid on;

pcolor(y*1e9,z*1e9,squeeze(V0(:,round(end/2),:))')
plot([-1 1]*10,[0 0],'w--','linewidth',2)
idz=find(z>QDh);
idz=idz(1)-1;
plot([-1 1]*10,[1 1]*z(idz)*1e9,'w--','linewidth',2)
plot([1 1]*0,[-1 1]*10,'b--','linewidth',2)
plot([1 1]*x(end)*1e9,[-1 1]*10,'r--','linewidth',2)

colormap(jet)
colorbar
shading flat
caxis([floor(min(V0(:)))  ceil(max(V0(:)))])

xlim([-1 1]*Mx/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('y (nm)')
ylabel('z (nm)')
title(strcat('Potential (eV): Vyz @x=0nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2,'fontsize',15)
hold on;grid on;

idz=find(z>0);
idz=idz(1);

pcolor(x*1e9,y*1e9,squeeze(V0(:,:,idz)))

colormap(jet)
colorbar
shading flat
caxis([floor(min(V0(:)))  ceil(max(V0(:)))])

xlim([-1 1]*My/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('y (nm)')
title(strcat('Potential (eV): Vxy @z=',num2str(z(idz)*1e9,'%.2f'),'nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,5,'fontsize',15)
hold on;grid on;

idz=find(z>QDh);
idz=idz(1)-1;
pcolor(x*1e9,y*1e9,squeeze(V0(:,:,idz)))

colormap(jet)
colorbar
shading flat
caxis([floor(min(V0(:)))  ceil(max(V0(:)))])

xlim([-1 1]*My/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('y (nm)')
title(strcat('Potential (eV): Vxy @z=',num2str(z(idz)*1e9,'%.2f'),'nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,3,'fontsize',15)
hold on;grid on;
plot(z*1e9,squeeze(V0(round(end/2),round(end/2),:)) ,'b.-')
plot(z*1e9,squeeze(V0(round(end/2),end,:)) ,'r.-')

ylim([floor(min(V0(:)))  ceil(max(V0(:)))])
xlabel('z (nm)')
ylabel('Energy (eV)')
title(strcat('Potential (eV): Vz @y=0nm and \color{blue}x=0nm ; \color{red}x=',num2str(x(end)*1e9),'nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,6,'fontsize',15)
    hold on;grid on;
    view (-38, 30);
    aa=[
	-a   -a   0
	-a   +a   0
	+a   +a   0
	+a   -a   0
	-a   -a   0
    ]*1e9;
    plot3(aa(:,1),aa(:,2),aa(:,3),'k','linewidth',2)

    bb=[
	-(h-QDh)*a/h     -(h-QDh)*a/h     QDh
	-(h-QDh)*a/h     +(h-QDh)*a/h     QDh
	+(h-QDh)*a/h     +(h-QDh)*a/h     QDh
	+(h-QDh)*a/h     -(h-QDh)*a/h     QDh
	-(h-QDh)*a/h     -(h-QDh)*a/h     QDh
    ]*1e9;
    plot3(bb(:,1),bb(:,2),bb(:,3),'k','linewidth',2)

    for j=1:length(aa(:,1))  
  	plot3([aa(j,1) bb(j,1)],[aa(j,2) bb(j,2)],[aa(j,3) bb(j,3)],'k','linewidth',2)
    end
    cc=[
	-2*a   -2*a   0
	-2*a   +2*a   0
	+2*a   +2*a   0
	+2*a   -2*a   0
	-2*a   -2*a   0
       ]*1e9;
    plot3(cc(:,1),cc(:,2),cc(:,3),'k','linewidth',2)

    dd=[
	-2*a   -2*a   -WLt
	-2*a   +2*a   -WLt
	+2*a   +2*a   -WLt
	+2*a   -2*a   -WLt
	-2*a   -2*a   -WLt
       ]*1e9;
    plot3(dd(:,1),dd(:,2),dd(:,3),'k','linewidth',2)

    for j=1:length(cc(:,1))
	plot3([cc(j,1) dd(j,1)],[cc(j,2) dd(j,2)],[cc(j,3) dd(j,3)],'k','linewidth',2)
    end
    hold on

M=max([Mx My]);
xlim([-1 1]*M*1e9)
ylim([-1 1]*M*1e9)
zlim([-0.9 0.9]*M*1e9)

xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title(strcat('Pyramidal Quantum Dot'))