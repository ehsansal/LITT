%--------------------------------------------------------------------------
%**********************************  MONTE CARLO SIMULATION  ******************
%----------------------------------------------------------------------

clear;
clc;
disp('                       ...start Monte Carlo...')
beep5;
tic

%%%%%%%%% parameters
N=50000;
layers=3;       %number of layers
d=[1.5,1.5,1.5];      %layers depth;
% d=3;
Nx=60;
Ny=60;
Nz=60;
Nr=60;
X_tissue=4.5;
Y_tissue=4.5;
R_tissue=0.5*sqrt(X_tissue*X_tissue+Y_tissue*Y_tissue);
applicator_center=[0,0,2.25];
applicator_length=1;
beam_width=0.2;
dx=X_tissue/Nx;
dy=Y_tissue/Ny;
dr=R_tissue/Nr;
dz=sum(d)/Nz;
coszero=1-(1e-6);
cos90=1e-6;
wth=0.001;
chance=0.1;

rand_s=0;
photonout_band=0;
photonout_top=0;
photonout_bottom=0;
photonout_roulette=0;
photonout_tissue=0;
photonout_fresnel=0;
photonout_x=0;
photonout_y=0;
photonout_z=0;
photonout_r=0;
photon_store=0;

% ua=[0.5,0.3,0.5];  %for 1064 NdYAG laser
ua=[0.8,0.8,0.8];     %for 850 Diode Laser
% ua=0.6;
% us=[169,109,169];   %for 1064 NdYAG laser
us=[130,130,130];   %for 850  Diode Laser
% us=108;
% g=[0.952,0.917,0.952];  %for 1064  1064 NdYAG laser
g=[0.955,0.955,0.955]; %for 850  Diode Laser
% g=0.902;
n=[1.33,1.33,1.33];
% n=1.33;
nair=1;
n_tissue=[nair,n,nair];
z0=[];
z1=[];
cos_critical0=[];
cos_critical1=[];
z0(1,1)=0;

%%%%%%%%%%% calculating layers depth into matrix
for l=1:layers
    z1(1,l)=z0(1,l)+d(1,l);
    z0(1,l+1)=z1(1,l);
end

%%%%%%%%%% calculating critical angles

for l=1:layers;
    a1=n_tissue(1,l);
    a2=n_tissue(1,l+1);
    a3=n_tissue(1,l+2);
    if a2>a1
        cos_critical0(1,l)=sqrt(1-(a1*a1)/(a2*a2));
    else
        cos_critical0(1,l)=0;
    end
    if(a2 >a3 )
        cos_critical1(1,l)=sqrt(1-(a3*a3)/(a2*a2));
    else
        cos_critical1(1,l)=0;
    end;
end;


%%%%%%%%% A matrix of final formulation in MC
A_xyz=zeros(Nx,Ny,Nz);
A_rz=zeros(Nr,Nz);

d_c=[];
d_c(1,1)=d(1,1);
for l=1:layers-1;
    d_c(1,l+1)=d_c(1,l)+d(1,l+1);
end

%%%%%%%% starting the program

for t=1:N
    w=1;
    dead=0;
    s=0;
    sleft=0;
    %%%%%%  for horizontal applicator
    z=0;
    y=0;
    x=0;
%             [x,y]=pol2cart(2*pi*rand,0.5*beam_wide);
    x=z+applicator_center(1,1);
    y=y+applicator_center(1,2);
    z=applicator_center(1,3)+applicator_length*(rand-0.5);

    layer=1;
    while z>d_c(1,layer)
        layer=layer+1;
    end
    cosTeta=2*rand-1;       %for g=0;
    sinTeta=sqrt(1-cosTeta*cosTeta);
    psi=2*pi*rand;
    uz=cosTeta;
    uy=sinTeta*sin(psi);
    ux=sinTeta*cos(psi);
    hit=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while dead==0;
        aa=ua(1,layer);
        as=us(1,layer);
        ut=aa+as;
        if sleft==0;
            rand_s=rand_s+1;
            s=-log(rand)/ut;
        else
            s=sleft/ut;
            sleft=0;
        end
        b1=z0(1,layer);
        b2=z1(1,layer);
        if uz>0;
            dl_b=(b2-z)/uz;     % the left step to go
        else
            dl_b=(b1-z)/uz;     % the left step to go
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if uz~=0 && s>dl_b
            sleft=(s-dl_b)*(ut);
            s=dl_b;
            hit=1;
        else
            hit=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if hit==1;
            x=x+s*ux;
            y=y+s*uy;
            z=z+s*uz;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if uz<0
                n1=n(1,layer);
                n0=n_tissue(1,layer);
                r=0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if -uz<=cos_critical0(1,layer),     %%%   without total reflection and uz=cosai
                    r=1;
                else
                    if n0==n1
                        r=0;
                        cosat=-uz;
                    elseif -uz>coszero
                        cosat=-uz;
                        r=(n0-n1)*(n0-n1)/((n0+n1)*(n0+n1));
                    elseif -uz<cos90;
                        cosat=0;
                        r=1;
                    else
                        sinai=sqrt(1-uz*uz);        %sin(alfai)
                        sinat=sinai*n1/n0;      %sin(alfat)
                        if sinat>1;
                            cosat=0;    %cos(alfat)
                            r=1;
                        else
                            cosat=sqrt(1-sinat*sinat);
                            cap=-uz*cosat-sinai*sinat;
                            cam=-uz*cosat+sinai*sinat;
                            sap=sinai*cosat-uz*sinat;
                            sam=sinai*cosat+uz*sinat;
                            r=0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rand > r;
                    %%%%%%%%%%%%%%%%%%%%%%%%%5
                    if layer==1;
                        dead=1;
                        photonout_top=photonout_top+w;
                    else
                        layer=layer-1;
                        ux=ux*n1/n0;
                        uy=uy*n1/n0;
                        uz=-cosat;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    uz=-uz;
                end
                %%%%%%%%%%%%%%%%%%%%
            else
                n1=n(1,layer);
                n2=n_tissue(1,layer+2);
                r=0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if uz<=cos_critical1(1,layer),   %%%   without total reflection and uz=cosai
                    r=1;
                else
                    if n1==n2
                        r=0;
                        cosat=uz;
                    elseif uz>coszero
                        cosat=uz;
                       r=(n2-n1)*(n2-n1)/((n2+n1)*(n2+n1));
                    elseif uz<cos90;
                        cosat=0;
                        r=1;
                    else
                        sinai=sqrt(1-uz*uz);        %sin(alfai)
                        sinat=sinai*n1/n2;      %sin(alfat)
                        if sinat>1;
                            cosat=0;    %cos(alfat)
                            r=1;
                        else
                            cosat=sqrt(1-sinat*sinat);
                            cap=uz*cosat-sinai*sinat;
                            cam=uz*cosat+sinai*sinat;
                            sap=sinai*cosat+uz*sinat;
                            sam=sinai*cosat-uz*sinat;
                            r=0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%
                if rand>r
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    if layer==layers
                        dead=1;
                        photonout_bottom=photonout_bottom+w;
                    else
                        layer=layer+1;
                        ux=ux*n1/n2;
                        uy=uy*n1/n2;
                        uz=cosat;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    uz=-uz;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if hit==0
            aa=ua(1,layer);
            as=us(1,layer);
            ut=aa+as;
            dw=w*aa/ut;
            w=w-dw;
            x=x+s*ux;
            y=y+s*uy;
            z=z+s*uz;
            ix=ceil(x/dx+Nx/2);
            iy=ceil(y/dy+Ny/2);
            iz=ceil(z/dz);
            ii=1;
            ir=ceil(sqrt(x*x+y*y)/dr);
            if ((ix > Nx)||(ix<1)), ii=0;photonout_x=photonout_x+dw;end;
            if ((iy > Ny)||(iy<1)), ii=0 ;photonout_y=photonout_y+dw; end;
            if ((iz > Nz)||(iz<1)), ii=0;photonout_z=photonout_z+dw ; end;
            if ((ir > Nr)||(ir<1)), ii=0 ;photonout_r=photonout_r+dw; end;
            if (ii)
                photon_store=photon_store+dw;
                A_xyz(ix,iy,iz)=A_xyz(ix,iy,iz)+dw;
                A_rz(ir,iz)=A_rz(ir,iz)+dw;
            else
                photonout_tissue=photonout_tissue+dw;
            end;


            %%%%%%%%%%%%%%%%%%%%%% new direction
            g0=g(1,layer);
            if g0==0;
                cosTeta=2*rand-1;
            else
                temp=(1-g0*g0)/(1-g0+2*g0*rand);
                cosTeta=(1+g0*g0-temp*temp)/(2*g0);
            end;
            sinTeta=sqrt(1-cosTeta*cosTeta);
            psi=2*pi*rand;
            cosPsi=cos(psi);
            if psi<pi
                sinPsi=sqrt(1-cosPsi*cosPsi);
            else
                sinPsi=-sqrt(1-cosPsi*cosPsi);
            end

            if abs(uz)>coszero
                ux=sinTeta*cosPsi;
                uy=sinTeta*sinPsi;
                uz=cosTeta*sign(uz);
            else
                temp=sqrt(1-uz*uz);
                ux=ux*cosTeta+sinTeta*(ux*uz*cosPsi-uy*sinPsi)/temp;
                uy=uy*cosTeta+sinTeta*(uy*uz*cosPsi+ux*sinPsi)/temp;
                uz=uz*cosTeta-sinTeta*cosPsi*temp;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rolette
        if w< wth && dead==0
            if w==0
                dead=1;
            elseif rand<chance
                w=w/chance;
            else
                dead=1;
                photonout_roulette=photonout_roulette+w;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
disp('                       ...End of Monte Carlo.')
beep5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
photonout_band=photonout_bottom+photonout_top+photonout_roulette;

A_xyzDensity=[];
A_xyzDensity=A_xyz./(0.01*dx*0.01*dy*0.01*dz*N);
A_rzDensity=[];
for ir=1:Nr
    da=2*pi*(ir-0.5)*0.01*dr*0.01*dr*0.01*dz*N;
    A_rzDensity(ir,:)=A_rz(ir,:)./da;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculating the matrix for comsole Multiphysics
Q= A_xyzDensity;
xplane=zeros(1,Nx);
yplane=zeros(1,Ny);
zplane=zeros(1,Nz);
for i=2:Nx
    xplane(1,i)=X_tissue/Nx+xplane(1,i-1);
    yplane(1,i)=Y_tissue/Nx+yplane(1,i-1);
    zplane(1,i)=sum(d)/Nz+zplane(1,i-1);
end
xyz=zeros(Nx*Ny+3,Nx);
xyz(1,:)=xplane(1,:)/100;
xyz(2,:)=yplane(1,:)/100;
xyz(3,:)=zplane(1,:)/100;

p=0;
for i=1:Nz
    for j=1:Ny
        for k=1:Nx
            xyz(4+p,k)=Q(k,j,i);     %copy xyz(1:3,:) to grid part and xyz(4:end,:) to data file
        end
        p=p+1;
    end
end
%%%%%%%%%%%%%%%%%%%%
DIM=[Nx,Ny,Nz];
dt=3600;
T_bloodi=37;
mwater=0.8;
Ci=ones(1,3)*3600;
Roi=ones(1,3)*1060;
% Roi=1300-300*mwater;  Roi=ones(1,NLayer)*Roi;
%       Ci=1550.3+2639.7*mwater;  Ci=ones(1,NLayer)*Ci;
Ldi=.0557+.5698*mwater;
Ldi=ones(1,layers)*Ldi;
Ld=[];
Ld=ones(DIM);
Ro=[];
Ro=ones(DIM);
DtRoC=[];
DtRoC=ones(DIM);
T_blood=[];
T_blood=ones(DIM).*(T_bloodi+273);
W_blood=[];
W_blood=ones(DIM);
HeatDensity=A_xyzDensity;
for layer=1:layers
    k1=floor(z0(1,layer)/dz)+1;
    k2=floor(z1(1,layer)/dz);
    Ld(:,:,k1:k2)=Ld(:,:,k1:k2).*Ldi(1,layer);
    DtRoC(:,:,k1:k2)=DtRoC(:,:,k1:k2).*(dt./(Roi(1,layer).*Ci(1,layer)));
    Ro(:,:,k1:k2)=Ro(:,:,k1:k2).*Roi(1,layer);
    %     T_M(:,:,k1:k2)=T_M(:,:,k1:k2).*(Ti(1,layer)+273);
    %     W_blood(:,:,k1:k2)=W_blood(:,:,k1:k2).*W_bloodi(1,layer);
end;
Q1=[]; 
Q1=HeatDensity.*DtRoC;

xyz1=zeros(Nx*Ny+3,Nx);
xyz1(1,:)=xplane(1,:);
xyz1(2,:)=yplane(1,:);
xyz1(3,:)=zplane(1,:);

p=0;
for i=1:Nz
    for j=1:Ny
        for k=1:Nx
            xyz1(4+p,k)=Q1(k,j,i);     %copy xyz1(1:3,:) to grid part and xyz(4:end,:) to data file
        end
        p=p+1;
    end
end
xyz2=xyz(4:end,:)/500;
xyz3=xyz(4:end,:)/600;
xyz4=xyz(4:end,:)/1000;

toc