%% =============================================%%
% This code is for TwoPhase 2D Pseudo-time method
% Based on 2018 Ludovic Raess, Thibault Duretz, Yury podladchikov
% Liuhao 2023.7.23


%% =============================================%%

clear

%% Material physics properties


rhofg = 1.0 ;      % fluid density * g
rhosg = 2*rhofg ;  % solid density * g
kmuf0 = 1.0;       % reference permeability
etas0 = 1.0;       % reference bulk viscosity
n     = 3;         % the exponent of porosity
ndim  = 2;         % the number of dimension(2-D)
C = 10;            % Rheological cst. bulk/shear viscosity ration Rheological cst. (bulk to shear viscosity)
phi0  = 0.01;      % reference porosity
phiA  = 2*phi0;    % amplitude of initial porosity perturbation
r     = 2;         % radius of initil porosity perturbation 
lam0  = 1.0;       % standard deviation of initial porosity perturbation
R     = 500;       % Rheological cst. (compaction versus decompaction)
LAMP  = 0.01;      % Effective pressure transition zone
dt    = 1e-5;      % physical time-step

%% model setup
xsize   = 20;
ysize   = 2*xsize;
nx      = 80;        % cell number of x direction
ny      = 2*nx;      % cell numebr of y direction
nt      = 1000;      % number of time steps
tol_nonlin  = 1e-5;  % non-linear tolerance
iterMax = 3e4;       % max non-linear iterations
miu_n   = 1.0;       % numerical compressibility
Vdmp    = 5.0;       % velocity damping for momentum equations
Pfdmp   = 1.4*xsize/nx;% fluid pressure damping for momentum equations
Vr      = 2;         % reduction of PT steps for velocity
Ptr     = 2;         % reduction of PT steps for total pressure
Pfr     = 8;         % reduction of PT steps for fluid pressure
rele    = 1;         % relaxation factor of non-linear viscosity
relk    = 1e-1;      % relaxation factor of non-linear permeability
dt_fact = 1e-3;      % reduction of physical timestep

%% Preprocessing
dx   = xsize/nx;   % grid step in x
dy   = ysize/ny;   % grid step in y
mus  = etas0*phi0/C; %% tricky?
lam  = lam0*sqrt(etas0*kmuf0); % porosity perturbation
rogBG = rhofg*phi0 + rhosg*(1-phi0); % rho*g of background
dtV  = min(dx,dy).^2./mus/(1+miu_n)/2.1/ndim/Vr;  % used for relaxation term. the momentum pseudo time-steps / reduction of PT steps for velocity
xc       = dx/2:dx:xsize-dx/2; % central of cell
yc       = dy/2:dy:ysize-dy/2; % central of cell
[Xc,Yc]  = ndgrid(xc,yc); % design the grid

% Arrays allocation
Pt       = zeros(nx  ,ny  ); % central of cell
Pf       = zeros(nx  ,ny  ); % central of cell
Sxy      = zeros(nx+1,ny+1); % on the node
Vx       = zeros(nx+1,ny  ); % on the grid 
Vy       = zeros(nx  ,ny+1); % on the grid
qDx      = zeros(nx+1,ny  ); % on the grid same with Vx
qDy      = zeros(nx  ,ny+1); % on the grid same with Vy
RPf      = zeros(nx  ,ny  ); % central of cell
Rx       = zeros(nx-1,ny  ); % on the grid same with Vx
Ry       = zeros(nx  ,ny-1); % on the grid same with Vy

%% Initial conditions
radc     = ((Xc - xsize*.5)/lam/4).^2 + ((Yc - ysize*.25)/lam).^2; % ellipse zone
phi      = phi0*ones(nx,ny);
phi(radc<1) = phi(radc<1) + phiA; % distribution of prosity
etas     = mus./(phi*C); % The bulk compaction viscosity
kmuf     = kmuf0*(phi./phi0); % the pporosity  denpendent permeability
phi0bc   = mean(phi(:,end));
qDy(:,[1 end]) = (rhosg-rhofg)*(1-phi0bc)*kmuf0*(phi0bc./phi0).^n;

%% Pseudo Time loop
for it = 1:nt
    phi_o = phi;
    for iter = 1:iterMax
        % Compute the physical properties circle
        rog     = rhofg*phi + rhosg*(1-phi) - rogBG;
        rog_avg = 0.5*(rog(:,1:end-1) + rog(:,2:end));
        etas    = etas*(1-rele) + rele*mus./phi*C.*(1+0.5*(1./R-1).*(1+tanh((Pf-Pt)./LAMP)));
        kmuf    = kmuf*(1-relk) + relk*kmuf0*(phi./phi0).^n;
        divVs   = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
        divqD   = diff(qDx,1,1)/dx + diff(qDy,1,2)/dy;
        phi     = phi_o + (1-phi).*divVs*dt;
        RPt     = divVs + (Pt-Pf)./etas./(1-phi);
        RPf     = RPf*(1-Pfdmp) + divqD - (Pt-Pf)./etas./(1-phi);
        Pt      = Pt - RPt.*2.1*ndim*mus*(1+miu_n)/ny/Ptr;
        % the fluid pressure pseudo-time step
        dtPf    = min(dx,dy).^2./kmuf./2.1./ndim./Pfr; 
        Sxx     = 2*mus*(diff(Vx,1,1)./dx - 1/3*divVs  + miu_n*RPt);
        Syy     = 2*mus*(diff(Vy,1,2)./dy - 1/3*divVs  + miu_n*RPt);
        Sxy(2:end-1,2:end-1) = mus*(diff(Vx(2:end-1,:),1,2)./dy + diff(Vy(:,2:end-1),1,1)./dx);
        
        
        % Compute the relaxation equation
        dtPf2   = dtPf;
        dtPf2(2:end-1,2:end-1) = min(min(min(dtPf(1:end-2,2:end-1), dtPf(3:end  ,2:end-1)), dtPf(2:end-1,2:end-1) ), min(dtPf(2:end-1,1:end-2), dtPf(2:end-1,3:end  )));
        % damp on fluid pressure equation
        Pf      = Pf - RPf.*dtPf2;
        Rx      = Rx*(1-Vdmp/nx) + (diff(Sxx,1,1)/dx + diff(Sxy(2:end-1,:),1,2)/dy - diff(Pt,1,1)/dx);
        Ry      = Ry*(1-Vdmp/ny) + (diff(Sxy(:,2:end-1),1,1)/dx + diff(Syy,1,2)/dy - diff(Pt,1,2)/dy - rog_avg);
        resid   = max([norm(Ry(:))/length(Ry(:)),norm(RPf(:))/length(RPf(:))]); 
        % check the residual
        if(max(resid)<tol_nonlin)
            break
        end 
        
        % calculate new velocity with relaxation term and new Darcy flow
        % k                k-1            update
        Vx(2:end-1,:)  = Vx(2:end-1,:) + Rx.*dtV;
        Vy(:,2:end-1)  = Vy(:,2:end-1) + Ry.*dtV;
        qDx(2:end-1,:) = -0.5*(kmuf(1:end-1,:)+kmuf(2:end,:)).* diff(Pf,1,1)/dx;
        qDy(:,2:end-1) = -0.5*(kmuf(:,1:end-1)+kmuf(:,2:end)).*(diff(Pf,1,2)/dy + (rhofg-rogBG));
    end
    
    dt = dt_fact/(1e-10+max(abs(divVs(:))));
    
    
    % draw the picture
    figure(1),clf,colormap(jet)
    subplot(131),imagesc(xc,yc,phi'),axis xy,axis image;
    colorbar,title(['it=' num2str(it)]),xlabel('\phi')
    
    subplot(132),imagesc(xc,yc,(Pt-Pf)'),axis xy,axis image;
    colorbar,title(['iter=' num2str(iter)]),xlabel('Pt-Pf')
    
    %figure(2),clf,colormap(jet)
    subplot(133),imagesc(xc,yc,log(phi'./phi0)),axis xy,axis image;
    colorbar,title(['iter=' num2str(iter)]),xlabel('log(\phi/\phi0)')

    drawnow
end
