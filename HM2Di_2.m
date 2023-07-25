%% --------------------------------------------------------
% This code is used to calculate two phase flow
%% --------------------------------------------------------

clear


noisy           = 1;
%% Grid setup
ra              = 2;                    %y/x
halfnode        = 42;                   % from the middle to one horizontal side                   
nx              = 2*halfnode-1;         % number of node in x
ny              = 2*halfnode*ra-1;      % number of node in y
nt              = 1000;                  % number of timestep
dt_red          = 1e-3;                 % reduction of explicit timestep for porosity update
dt_fact         = 1;
dt              = 1e-5;                 % initial timestep
xsize           = 20;                   
ysize           = xsize*ra;
dx              = xsize/nx;
dy              = ysize/ny;
% Non-linear solver
CN              = 0.0;                  % the Crank-Nicolson switch from 0 to 0.5
coupled         = 0;                    % optimised(0) or basic(1) solution
tol_nonlin      = 1e-8;                 % nonlinear rolerance
nmax_nonlin     = 100;                  % max number of nonlinear iterations
% Linear solver
tol_linv        = tol_nonlin/1e3;       % linear tolerance momentum
tol_linpt       = tol_nonlin/1e3;       % linear tolerance for total pressure
tol_linpf       = tol_nonlin/1e3;       % linear tolerance for fluid pressure
nPH             = 100;                  % Pre-conditionning calculate times

%% material properties
rhofg           = 1;                    % fluid rho*g
rhosg           = 2*rhofg;              % solid rho*g
kmuf0           = 1;                    % background permeability
etas0           = 1;                    % background bulk viscosity
phi0            = 1e-2;                 % background porosity
n               = 3;                    % exponent of porosity
ndim            = 2;                    % number of dimension
eta2mus         = 10;                   % bulk to shear viscosity ratio
R               = 500;                  % bilinear viscosity reduction
lamp            = 0.01;                 % lambda pe
lam0            = 1;                    % perubation standard devition
phiA            = 2*phi0;               % porosity amplitude
mus0            = etas0*phi0/eta2mus;   % back ground shear viscosity
lam             = lam0*sqrt(etas0*kmuf0);
rogbg           = rhofg.*phi0 + (1-phi0).*rhosg;
%% Initial conditions
Vx              = zeros(nx+1,ny  );
Vy              = zeros(nx  ,ny+1);
Pf              = zeros(nx  ,ny  );
Pt              = zeros(nx  ,ny  );
phi             = phi0*ones(nx  ,ny  );
divV            = zeros(nx  ,ny  );
[XC,YC]         = ndgrid(dx/2:dx:xsize-dx/2,dy/2:dy:ysize-dy/2);
ndif            = 0.3;
arcx            = 0.25;
arcy            = 1.0;
rad             = ((XC-xsize*.5)/ lam*arcx).^2+((YC-ysize*.2)/ lam*arcy).^2;
phi(rad<1)      = phi(rad<1) + phiA;
dtD             = min(dx,dy).^2/(2.1*ndim);
ntD             = fix(ndif/dtD);
for idif = 1:ntD
    phi(2:end-1,2:end-1) = phi(2:end-1,2:end-1) + dtD*(diff(phi(:,2:end-1),2,1)/dx + diff(phi(2:end-1,:),2,2)/dy);
end

%% Numbering Pt and Vx Vy
NumVx = reshape(1:(nx+1)*ny,nx+1,ny  );
NumVy = reshape(1:nx*(ny+1),nx  ,ny+1);
NumPt = reshape(1:nx*ny    ,nx  ,ny  );
NumPf = reshape(1:nx*ny    ,nx  ,ny  );
% max means index after last property
NumVyG= NumVy + max(NumVx(:));
%NumPtG= NumPt + max(NumVyG(:));
%NumPfG= NumPf + max(NumPtG(:));

%% Boundary Conditions of velocities
% Vx 
ibcVxL = NumVx(1      ,:)';
ibcVxR = NumVx(end    ,:)';
ibcVxU = NumVx(2:end-1,1) ;
ibcVxB = NumVx(2:end-1,end);
% Vy 
ibcVyL = NumVyG(1      ,2:end-1)';
ibcVyR = NumVyG(end    ,2:end-1)';
ibcVyU = NumVyG(:      ,1      ) ;
ibcVyB = NumVyG(:      ,end    ) ;

Vx_L  = zeros(1,ny)';
Vx_R  = zeros(1,ny)';
Vy_U  = zeros(nx,1) ;
Vy_B  = zeros(nx,1) ;

%% BC of fluid pressure
Pf_L  = 0*ones(1 ,ny)'; % BC Dirichlet values
Pf_R  = 0*ones(1 ,ny)';
Pf_U  = 0*ones(nx,1 );
Pf_B  = 0*ones(nx,1 );
%% Darcy flux BC 
phi0bc = mean(phi(:,end));
val    = ((rhosg-rhofg)*(1-phi0bc))*kmuf0*(phi0bc./phi0).^n;
qDx_L  = zeros(1,ny)';
qDx_R  = zeros(1,ny)';
qDy_U  = val*ones(nx,1);
qDy_B  = val*ones(nx,1);

% BC type: 1=Dirichlet , 2=Neumann
tbcL   = ones(1,ny)'; tbcL(:) = 2; 
tbcR   = ones(1,ny)'; tbcR(:) = 2; 
tbcU   = ones(nx,1) ; tbcU(:) = 2; 
tbcB   = ones(nx,1) ; tbcB(:) = 2; 

% [W E S N]
ibcPfL = NumPf(1,:)';   ibcPfR = NumPf(end,:)';  ibcPfU = NumPf(:,1);  ibcPfB = NumPf(:,end);
id     = [1,length(ibcPfL),length(ibcPfL)+length(ibcPfR),length(ibcPfL)+length(ibcPfR)+length(ibcPfU),length(ibcPfL)+length(ibcPfR)+length(ibcPfU)+length(ibcPfB)];
tbc    = [tbcL      ;tbcR          ;tbcU           ;tbcB       ];
vBcD   = [Pf_L      ;Pf_R          ;Pf_U           ;Pf_B       ];
vBcN   = [qDx_L     ;qDx_R         ;qDy_U          ;qDy_B      ];
vBc    = zeros(size(tbc)); vBc(tbc==1) = vBcD(tbc==1); vBc(tbc==2) = vBcN(tbc==2);


%% Assemble constant blocks and BC's
iVxC   = NumVx;
iPtL   = ones(size(iVxC));      iPtL(1:end-1,:) = NumPt;
iPtR   = ones(size(iVxC));      iPtR(2:end  ,:) = NumPt;
cPtL   = ones(size(iVxC))/dx;   cPtL([1 end],:) = 0;
cPtR   = -ones(size(iVxC))/dx;  cPtR([1 end],:) = 0;

Idx    = [ iVxC(:); iVxC(:)];
Jdx    = [ iPtL(:); iPtR(:)];
Vdx    = [ cPtL(:); cPtR(:)];

iVyC   = NumVyG;
iPtU   = ones(size(iVyC));      iPtU(:,1:end-1) = NumPt;
iPtB   = ones(size(iVyC));      iPtB(:,2:end  ) = NumPt;
cPtU   = ones(size(iVyC))/dy;   cPtU(:,[1 end]) = 0;% the first and last column is zero
cPtB   = -ones(size(iVyC))/dy;  cPtB(:,[1 end]) = 0;

Idy    = [ iVyC(:); iVyC(:) ]';
Jdy    = [ iPtU(:); iPtB(:) ]';
Vdy    = [ cPtU(:); cPtB(:) ]';


grad   = sparse( [Idx(:); Idy(:)], [Jdx(:); Jdy(:)], [Vdx(:); Vdy(:)], (nx+1)*ny+(ny+1)*nx, nx*ny );

% spy(grad)

div    = -grad';


% Block Coefficients on Div
BCD    = zeros(size(div,1),1);
BCD(NumPt(1   ,:  )) = BCD(NumPt(1    ,:  )) + 1/dx*Vx_L;
BCD(NumPt(end ,:  )) = BCD(NumPt(end  ,:  )) - 1/dx*Vx_R;
BCD(NumPt(:   ,1  )) = BCD(NumPt(:    ,1  )) + 1/dy*Vy_U;
BCD(NumPt(:   ,end)) = BCD(NumPt(:    ,end)) - 1/dy*Vy_B;

%% Linear shear viscosity
mus    = mus0*ones(nx,ny);
mus_axy= 0.25*(mus(1:end-1,1:end-1)+mus(2:end,1:end-1)+mus(1:end-1,2:end)+mus(2:end,2:end));
musv   = zeros(nx+1,ny+1);
musv(2:end-1,2:end-1) = mus_axy;
musv([1,end],:) = musv([2,end-1],:);
musv(:,[1,end]) = musv(:,[2,end-1]);
v = [Vx(:);Vy(:)];
pt= Pt(:);
pf= Pf(:);
%% Timeloop
time = 0;
nfail= 0;
for it = 1:nt
    % start the loop
    iter = 0;
    fprintf('\n**** Nonlin iter %d **** \n', iter );
    resnlv  = 2*tol_nonlin;
    resnlpt = 2*tol_nonlin;
    resnlpf = 2*tol_nonlin;
    phi_o   = phi;
    
    
    dv = zeros(max(NumVyG(:)),1);
    dpt = zeros(max(NumPt(:)),1);
    dpf = zeros(max(NumPf(:)),1);
    
    
    ddv = zeros(max(NumVyG(:)),1);
    ddpt = zeros(max(NumPt(:)),1);
    ddpf = zeros(max(NumPf(:)),1);
    
    rxsc = zeros(max(NumVyG(:)),1);
    rysc = zeros(max(NumPt(:)),1);
    s    = zeros(max(NumPt(:)),1);
    
    % save previous solutions
    divVo = divV;
    
    % Non-linear iterations
    while (resnlv > tol_nonlin || resnlpt > tol_nonlin || resnlpf >tol_nonlin) && iter < nmax_nonlin
        iter = iter + 1;
        fprintf('\n**** Nonlin iter %d **** \n', iter );
        divV = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
        dphidt = (1-phi).*((1-CN)*divV + CN*divVo);
        % Timestep selection
        if (it>1 && iter == 1)
            dt_vol = dt_red/(1e-10 + max(abs(dphidt(:))));
            dt_advS= 0.5*dy/(1e-10 + max(abs(Vy(:))));
            dt     = dt_fact*min(dt_vol , dt_advS);
            
        end
        
        
        % Non linear
        phi = phi_o + dt*dphidt;
        etaphi = mus./phi*eta2mus .* (1+0.5*(1./R-1).*(1+tanh((Pf-Pt)./lamp)));
        kmuf  = kmuf0 .* (phi/phi0).^n;
        kU = zeros(size(Pf));   kU(:,2:end  ) = (kmuf(:,1:end-1) + kmuf(:,2:end))/2;
        kB = zeros(size(Pf));   kB(:,1:end-1) = (kmuf(:,1:end-1) + kmuf(:,2:end))/2;
        kgy   = ((rhofg - rogbg) * (kB-kU))/dy;
        rog   = rhofg.*phi + (1-phi) .* rhosg - rogbg;
        
        
        % Assemble matrix with variable coeff
        %% Block PF
        kL = zeros(size(Pf));   kL(2:end,:  ) = (kmuf(1:end-1,:) + kmuf(2:end,:))/2;
        kR = zeros(size(Pf));   kR(1:end-1,:) = (kmuf(1:end-1,:) + kmuf(2:end,:))/2;
        kU = zeros(size(Pf));   kU(:,2:end  ) = (kmuf(:,1:end-1) + kmuf(:,2:end))/2;
        kB = zeros(size(Pf));   kB(:,1:end-1) = (kmuf(:,1:end-1) + kmuf(:,2:end))/2;
        
        iPf= NumPf;
        iPfL = ones(size(Pf)); iPfL(2:end  ,:) = NumPf(1:end-1,:);
        iPfR = ones(size(Pf)); iPfR(1:end-1,:) = NumPf(2:end  ,:);
        iPfU = ones(size(Pf)); iPfU(:,2:end  ) = NumPf(:,1:end-1);
        iPfB = ones(size(Pf)); iPfB(:,1:end-1) = NumPf(:,2:end  );
        
        cPfC  = (kL+kR)/dx/dx + (kU+kB)/dy/dy + 1./(etaphi.*(1-phi));
        cPfL  = -kL/dx/dx;
        cPfU  = -kU/dy/dy;
        
        % value not used
        %cPfR  = -kR/dx/dx;
        %cPfB  = -kB/dy/dy;
        
        I     = [  iPf(:) ;   iPf(:) ;   iPf(:) ;]';
        J     = [  iPf(:) ;   iPfL(:);   iPfU(:);]';
        V     = [  cPfC(:);   cPfL(:);   cPfU(:);]';
        
        PF    = sparse(I,J,V, nx*ny , nx*ny);
        
        %% PP block
        iPt   = NumPt;
        I     = iPt(:)';
        J     = I;
        V     = ( 1./(etaphi(:).*(1-phi(:))))';
        PP    = sparse(I,J,V);
        
        
        %% PfP block
        iPt = NumPt;  
        I = iPt(:)';  
        J = I;
        V = ( -1./(etaphi(:).*(1-phi(:))) )';
        PfP = sparse(I,J,V);
        
        
        %% VX block
        mus_L = zeros(size(Vx));   mus_L(2:end  ,:) = mus;
        mus_R = zeros(size(Vx));   mus_R(1:end-1,:) = mus;
        mus_U = zeros(size(Vx));   mus_U(2:end-1,1:end-1) = mus_axy;
        mus_B = zeros(size(Vx));   mus_B(2:end-1,2:end) = mus_axy;
        
        iVx   = NumVx;
        iVxL  = NumVx(1:end-1, :     );
        iVxR  = NumVx(2:end  , :     );
        iVxU  = NumVx(2:end-1,1:end-1);
        iVxB  = NumVx(2:end-1,2:end  );
        cVxC  = (mus_L+mus_R)/dx/dx + (mus_U+mus_B)/dy/dy;
        scVx  = max(cVxC(:));                   cVxC([1,end],:) = scVx;
        cVxL  = -mus_L(2:end  , :     )/dx/dx;  cVxL([1,end],:) = 0;
        cVxU  = -mus_U(2:end-1,1:end-1)/dy/dy; 
        Ivx   = [  iVx(:); iVxR(:); iVxB(:) ]';
        Jvx   = [  iVx(:); iVxL(:); iVxU(:) ]';
        Vvx   = [ cVxC(:); cVxL(:); cVxU(:) ]';
        
        %% VY Block
        mus_L = zeros(size(Vy));  mus_L(1:end-1,2:end-1) = mus_axy;
        mus_R = zeros(size(Vy));  mus_R(2:end  ,2:end-1) = mus_axy;
        mus_U = zeros(size(Vy));  mus_U( :     ,2:end  ) = mus;
        mus_B = zeros(size(Vy));  mus_B( :     ,1:end-1) = mus;
        rogy  = zeros(size(Vy));  rogy(  :     ,2:end-1) = 0.5*(rog(:,1:end-1)+rog(:,2:end));
        iVy   = NumVyG;
        iVyL  = NumVyG(1:end-1,2:end-1);
        iVyR  = NumVyG(2:end  ,2:end-1);
        iVyU  = NumVyG( :     ,1:end-1);
        iVyB  = NumVyG( :     ,2:end  );
        cVyC  = (mus_L+mus_R)/dx/dx + (mus_U+mus_B)/dy/dy;
        scVy  = max(cVyC(:));                   
        cVyC(:,[1,end]) = scVy;
        cVyL  = -mus_L(1:end-1,2:end-1)/dx/dx;
        cVyU  = -mus_U( :     ,2:end  )/dy/dy;  
        cVyU(:,[1,end]) = 0;
        Ivy   = [  iVy(:); iVyR(:); iVyB(:) ]';
        Jvy   = [  iVy(:); iVyL(:); iVyU(:) ]';
        Vvy   = [ cVyC(:); cVyL(:); cVyU(:) ]';


        %% VXY Block
        mus_L = zeros(size(Vy));  mus_L( :, :     ) = musv(1:end-1, :);
        mus_R = zeros(size(Vy));  mus_R( :, :     ) = musv(2:end  , :);
        mus_U = zeros(size(Vy));  mus_U( :,2:end  ) = mus;
        mus_B = zeros(size(Vy));  mus_B( :,1:end-1) = mus;
        iVy   = NumVyG( :    ,2:end-1);
        iVxLU = NumVx(1:end-1,1:end-1);
        iVxRU = NumVx(2:end  ,1:end-1);
        iVxLB = NumVx(1:end-1,2:end  );
        iVxRB = NumVx(2:end  ,2:end  );
        cVxLU = (-mus_L(:,2:end-1) + mus_U(:,2:end-1))/(dx*dy);  cVxLU(1  ,:) = 0;
        cVxRU = ( mus_R(:,2:end-1) - mus_U(:,2:end-1))/(dx*dy);  cVxRU(end,:) = 0;
        cVxLB = ( mus_L(:,2:end-1) - mus_B(:,2:end-1))/(dx*dy);  cVxLB(1  ,:) = 0;
        cVxRB = (-mus_R(:,2:end-1) + mus_B(:,2:end-1))/(dx*dy);  cVxRB(end,:) = 0;
        Ivxy   = [   iVy(:);   iVy(:);   iVy(:);   iVy(:) ]';
        Jvxy   = [ iVxLU(:); iVxRU(:); iVxLB(:); iVxRB(:) ]';
        Vvxy   = [ cVxLU(:); cVxRU(:); cVxLB(:); cVxRB(:) ];

        %% Sparse and Block Coefficient's
        K =  sparse( [Ivx(:); Ivy(:); Ivxy(:)], [Jvx(:); Jvy(:); Jvxy(:)], [Vvx(:); Vvy(:); Vvxy(:)], (nx+1)*ny+(ny+1)*nx, (nx+1)*ny+(ny+1)*nx );

        % BC's on velocity
        BCK                  = zeros(size(K,1),1);
        BCK(NumVx(2    ,:))  = BCK(NumVx(2    ,:))  + mus(1  ,:  )'/dx/dx.*Vx_L;
        BCK(NumVx(end-1,:))  = BCK(NumVx(end-1,:))  + mus(end,:  )'/dx/dx.*Vx_R;
        BCK(NumVyG(:, 2   )) = BCK(NumVyG(:, 2   )) + mus(:  ,1  ) /dy/dy.*Vy_U;
        BCK(NumVyG(:,end-1)) = BCK(NumVyG(:,end-1)) + mus(:  ,end) /dy/dy.*Vy_B;

        BCK(ibcVxB) = BCK(ibcVxB) +  (  mus(1:end-1,end) - musv(2:end-1,end) )/dx/dy.*Vy_B(1:end-1) ;   % VyNW
        BCK(ibcVxB) = BCK(ibcVxB) +  ( -mus(2:end  ,end) + musv(2:end-1,end) )/dx/dy.*Vy_B(2:end  ) ;   % VyNE
        BCK(ibcVxU) = BCK(ibcVxU) +  ( -mus(1:end-1,1  ) + musv(2:end-1, 1 ) )/dx/dy.*Vy_U(1:end-1) ;   % VySW
        BCK(ibcVxU) = BCK(ibcVxU) +  (  mus(2:end  ,1  ) - musv(2:end-1, 1 ) )/dx/dy.*Vy_U(2:end  ) ;   % VySE
        BCK(ibcVyR) = BCK(ibcVyR) + ((  mus(end,1:end-1) - musv(end,2:end-1) )/dx/dy.*Vx_R(1:end-1)')'; % VxSE
        BCK(ibcVyR) = BCK(ibcVyR) + (( -mus(end,2:end  ) + musv(end,2:end-1) )/dx/dy.*Vx_R(2:end  )')'; % VxNE
        BCK(ibcVyL) = BCK(ibcVyL) + (( -mus(1  ,1:end-1) + musv(1  ,2:end-1) )/dx/dy.*Vx_L(1:end-1)')'; % VxSW
        BCK(ibcVyL) = BCK(ibcVyL) + ((  mus(1  ,2:end  ) - musv(1  ,2:end-1) )/dx/dy.*Vx_L(2:end  )')'; % VxNW
        
        BCK(NumVyG(:)) = BCK(NumVyG(:)) - rogy(:);  % Gravity term in rhs
        BCK([ibcVxL; ibcVxR; ibcVyU; ibcVyB]) = [Vx_L*scVx; Vx_R*scVx; Vy_U*scVy; Vy_B*scVy];

        %% BC's on Pf
        BcPf  = zeros(nx*ny,1);
        % Boundary coefficient for Boundary condition
        dD    = zeros(nx*ny,1);
        % second Boundary condition coefficient
        cN    = zeros(nx*ny,1);
        % first Boundary condition coefficient
        cD    = zeros(nx*ny,1);
        % Pf diag
        d0    = spdiags(PF,0);
        kmv   = kmuf(:);
        
        dD(ibcPfL(tbcL==1)) = -kmv(ibcPfL(tbcL==1))/dx/dx;   % Dirichlet
        dD(ibcPfL(tbcL==2)) =  ones(size(ibcPfL(tbcL==2)))/dx; % Neumann
        cD(ibcPfL) = dD(ibcPfL).*vBc(id(1):id(2));
        cN(ibcPfL) = dD(ibcPfL).*vBc(id(1):id(2));
        d0(ibcPfL(tbcL==1))   = d0(ibcPfL(tbcL==1))   - 2.*dD(ibcPfL(tbcL==1));  
        BcPf(ibcPfL(tbcL==1)) = BcPf(ibcPfL(tbcL==1)) - 2.*cD(ibcPfL(tbcL==1));
        BcPf(ibcPfL(tbcL==2)) = BcPf(ibcPfL(tbcL==2)) +    cN(ibcPfL(tbcL==2));
        
        dD(ibcPfR(tbcR==1)) = -kmv(ibcPfR(tbcR==1))/dx/dx;   % Dirichlet
        dD(ibcPfR(tbcR==2)) = -ones(size(ibcPfR(tbcR==2)))/dx; % Neumann
        cD(ibcPfR) = dD(ibcPfR).*vBc(id(2)+1:id(3));
        cN(ibcPfR) = dD(ibcPfR).*vBc(id(2)+1:id(3));
        d0(ibcPfR(tbcR==1))   = d0(ibcPfR(tbcR==1))   - 2.*dD(ibcPfR(tbcR==1));
        BcPf(ibcPfR(tbcR==1)) = BcPf(ibcPfR(tbcR==1)) - 2.*cD(ibcPfR(tbcR==1));
        BcPf(ibcPfR(tbcR==2)) = BcPf(ibcPfR(tbcR==2)) +    cN(ibcPfR(tbcR==2));
        
        dD(ibcPfU(tbcU==1)) = -kmv(ibcPfU(tbcU==1))/dy/dy;   % Dirichlet
        dD(ibcPfU(tbcU==2)) =  ones(size(ibcPfU(tbcU==2)))/dy; % Neumann
        cD(ibcPfU) = dD(ibcPfU).*vBc(id(3)+1:id(4));
        cN(ibcPfU) = dD(ibcPfU).*vBc(id(3)+1:id(4));
        d0(ibcPfU(tbcU==1))   = d0(ibcPfU(tbcU==1))   - 2.*dD(ibcPfU(tbcU==1));
        BcPf(ibcPfU(tbcU==1)) = BcPf(ibcPfU(tbcU==1)) - 2.*cD(ibcPfU(tbcU==1));
        BcPf(ibcPfU(tbcU==2)) = BcPf(ibcPfU(tbcU==2)) +    cN(ibcPfU(tbcU==2));
        
        dD(ibcPfB(tbcB==1)) = -kmv(ibcPfB(tbcB==1))/dy/dy;   % Dirichlet
        dD(ibcPfB(tbcB==2)) = -ones(size(ibcPfB(tbcB==2)))/dy; % Neumann
        cD(ibcPfB) = dD(ibcPfB).*vBc(id(4)+1:id(5));
        cN(ibcPfB) = dD(ibcPfB).*vBc(id(4)+1:id(5));
        d0(ibcPfB(tbcB==1))   = d0(ibcPfB(tbcB==1))   - 2.*dD(ibcPfB(tbcB==1));
        BcPf(ibcPfB(tbcB==1)) = BcPf(ibcPfB(tbcB==1)) - 2.*cD(ibcPfB(tbcB==1));
        BcPf(ibcPfB(tbcB==2)) = BcPf(ibcPfB(tbcB==2)) +    cN(ibcPfB(tbcB==2));
        BcPf = BcPf + kgy(:);
        PF   = spdiags(d0,0,PF);

        %% Build full from lower
        PF  = PF + PF' - diag(diag(PF));
        K   = K  + K'  - diag(diag(K ));
        % Linearised operator
        bv=BCK; bpt=BCD; bpf=BcPf;
        fv  = -(   K*v + grad*pt          - bv  );  resnlv  = norm(fv)/length(fv); % nonlinear residuals
        fpt = -( div*v +   PP*pt + PfP*pf - bpt );  resnlpt = norm(fpt)/length(fpt);
        fpf = -(          PfP*pt +  PF*pf - bpf );  resnlpf = norm(fpf)/length(fpf);
        if iter==1
            resnlv0=resnlv; 
            resnlpt0=resnlpt; 
            resnlpf0=resnlpf; 
        end
        if noisy>=1, fprintf('Chk: NonLin res. ||res.v ||=%2.4e, ||res.v /res.v0 ||=%2.4e\n', resnlv,  resnlv/resnlv0   );
                     fprintf('Chk: NonLin res. ||res.pt||=%2.4e, ||res.pt/res.pt0||=%2.4e\n', resnlpt, resnlpt/resnlpt0 );
                     fprintf('Chk: NonLin res. ||res.pf||=%2.4e, ||res.pf/res.pf0||=%2.4e\n', resnlpf, resnlpf/resnlpf0 ); 
        end
        if (resnlv/resnlv0 < tol_nonlin)
            break; 
        end
        AJ=K  ;  BJ=grad;
        DJ=div;  EJ=PP  ;  FJ=PfP;
                 HJ=PfP ;  IJ=PF ;
        
        
        %% Solver
        if coupled == 1
            f   = [fv;    fpt;   fpf];
            Ms  = [AJ     BJ     0  ; ...
                   DJ     EJ     FJ ; ...
                   0      HJ     IJ ];
               
            scm = spdiags(mean(abs(Ms),2), 0, size(Ms,1), size(Ms,2));
            ds  = scm*( umfpack((scm*Ms*scm),'\',(scm*f)) );
            rJ  = Ms*ds-f;
            dv  = ds(1:length(dv));
            dpt = ds(length(dv) + 1:length(dv) + length(dpt));
            dpf = ds(length(dv) + length(dpt)  + 1:end);
            ddpf = PF\bpf;
            
            if noisy>=1, fprintf('\n Resid Jac = %2.2e, Resid Pf = %2.2e \n\n', norm(rJ) ,norm(PF*ddpf-bpf)); end
        else
            % Pre-conditionning (~Jacobi)
            Esc = EJ - FJ*(spdiags(diag(IJ ),0,size(PF,1),size(PF,2))\HJ);   % Esc = E - F*(diag(I  )\H) = E - F*( (1/diag(I  ))*H)
            Asc = AJ - BJ*(spdiags(diag(Esc),0,size(PF,1),size(PF,2))\DJ);   % Dsc % Asc = A - B*(diag(Esc)\D) = A - B*( (1/diag(Esc))*D)
            [Icc ,pI,sI] = chol(IJ ,'lower','vector');                       % Cholesky factors
            [Ascc,pA,sA] = chol(Asc,'lower','vector');                       % Cholesky factors
            Esci         = spdiags(1./diag(Esc),0,size(EJ,1),size(EJ,2));    % trivial inverse
            rx0=resnlv; 
            ry0=resnlpt; 
            rz0=resnlpf;
            for itPH = 1:nPH
                rx  = -( AJ*dv + BJ*dpt          - fv  );
                ry  = -( DJ*dv + EJ*dpt + FJ*dpf - fpt );
                rz  = -(         HJ*dpt + IJ*dpf - fpf );
        
                s(sI,1) = Icc'\(Icc\rz(sI));
                rysc      = -( FJ*s - ry );
                s         =    Esci*rysc;
                rxsc      = -( BJ*s - rx );
                
                ddv(sA,1) = Ascc'\(Ascc\rxsc(sA));

                s         = -( DJ*ddv - rysc );
                ddpt      =    Esci*s;
                s         = -( HJ*ddpt - rz );
                
                ddpf(sI,1) = Icc'\(Icc\s(sI));
                
                dv  = dv  + ddv; % Updates
                dpt = dpt + ddpt;
                dpf = dpf + ddpf;
                
                if noisy>1, fprintf('  --- iteration %d --- \n',itPH);
                            fprintf('  ||res.v ||=%2.2e\n', norm(rx)/length(rx));
                            fprintf('  ||res.pt||=%2.2e\n', norm(ry)/length(ry));
                            fprintf('  ||res.pf||=%2.2e\n', norm(rz)/length(rz));
                end
                if ((norm(rx)/length(rx)) < tol_linv) && ((norm(ry)/length(ry)) < tol_linpt) && ((norm(rz)/length(rz)) < tol_linpf), break; end
                if ((norm(rx)/length(rx)) > (norm(rx0)/length(rx0)) && norm(rx)/length(rx) < tol_nonlin && (norm(ry)/length(ry)) > (norm(ry0)/length(ry0)) && norm(ry)/length(ry) < tol_nonlin && (norm(rz)/length(rz)) > (norm(rz0)/length(rz0)) && norm(rz)/length(rz) < tol_nonlin)
                    fprintf(' > Linear residuals do no converge further:\n'); break;
                end
                rx0=rx; 
                ry0=ry; 
                rz0=rz; 
                if (itPH==nPH)
                    nfail=nfail+1; 
                end
            end%itPH
            if noisy>=1
                fprintf(' - Linear res. ||res.v ||=%2.2e\n', norm(rx)/length(rx) );
                fprintf(' - Linear res. ||res.pt||=%2.2e\n', norm(ry)/length(ry) );
                fprintf(' - Linear res. ||res.pf||=%2.2e\n', norm(rz)/length(rz) ); 
            end
        end%coupled
        [alpha] = LineSearch(nx,ny,dx,dy,dt,CN,divVo,v,pt,pf,dv,dpt,dpf,rhosg,rhofg,rogbg,phi,phi_o,phi0,kmuf0,n,eta2mus,mus,musv,R,lamp,noisy,resnlv0,resnlpt0,resnlpf0,resnlv,resnlpt,resnlpf,qDx_L,qDx_R,qDy_B,qDy_U);
        v   = v  + alpha*dv; % Solution from Newton iters
        pt  = pt + alpha*dpt;
        pf  = pf + alpha*dpf;
        
        filename3 = [pwd,  '/',  'DI-fig',  num2str(it),'-', num2str(nfail) ]
        figure(100),if iter==1, clf; end
        hold on, plot(iter,log10(resnlv),'-rs', iter,log10(resnlpt),'-+', iter,log10(resnlpf),'-o')
        title(['nfail = ' num2str(nfail)]),grid on%,drawnow
        close(gcf)
        Vx   = reshape(v(NumVx(:)) ,[nx+1,ny  ]);
        Vy   = reshape(v(NumVyG(:)),[nx  ,ny+1]);
        Pt   = reshape(pt          ,[nx  ,ny  ]);
        Pf   = reshape(pf          ,[nx  ,ny  ]);
    end%while nonlin
    time = time + dt;
    fprintf('\n>------------------- TIMESTEP %d  -> dt=%2.2e  -> time=%2.2e ----------------------------------------\n',it,dt,time)
        
    qDx = -(kmuf(1:end-1,:) + kmuf(2:end,:))/2.* diff(Pf,1,1)/dx;
    qDy = -(kmuf(:,1:end-1) + kmuf(:,2:end))/2.*(diff(Pf,1,2)/dy + rhofg-rogbg);
    Pe  = Pt-Pf;
    
    filename1 = [pwd,  '/',  'DI-fig2',  num2str(it) ]
    b=figure(2),set(b,'Visible','off'),clf,colormap(jet), tp ={'Pt';'Vx';'Vy';'Pf';'qDx';'qDy'};
    for ipl=1:size(tp,1), eval(['Aname = ',tp{ipl},';']); subplot(3,2,ipl), imagesc(flipud(Aname')), title(tp{ipl}), colorbar, axis image; end
    %drawnow
    print(gcf,'-djpeg','-r600',filename1);
    close(gcf);
    
    filename2 = [pwd,  '/',  'DI-fig4',  num2str(it) ]
    c = figure(4),set(c,'Visible','off'),clf,colormap(jet), tp ={'phi';'log10(kmuf)';'Pe';'etaphi'};
    for ipl=1:size(tp,1), eval(['Aname = ',tp{ipl},';']); subplot(2,2,ipl), imagesc(flipud(Aname')), title(tp{ipl}), colorbar, axis image; end
    %drawnow
    print(gcf,'-djpeg','-r600',filename2);
    close(gcf);

end


