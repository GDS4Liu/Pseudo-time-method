function [alpha] = LineSearch(nx,ny,dx,dy,dt,CN,divVo,v,pt,pf,dv,dpt,dpf,rhosg,rhofg,rogbg,phi,phi_o,phi0,kmuf0,n,eta2mus,mus,musv,R,lamp,noisy,resnlv0,resnlpt0,resnlpf0,nRMe0,nRCTe0,nRCFe0,qDx_L,qDx_R,qDy_B,qDy_U)
    nsteps = 20;
    amin   = 0.1;
    amax   = 2.0;
    dalpha = (amax-amin)/(nsteps-1);
    alphav = amin:dalpha:amax;
    nVx    = (nx+1)*ny;
    % momentum
    nRMe   = zeros(nsteps,1);
    % total pressure
    nRCTe  = zeros(nsteps,1);
    % fluid pressure
    nRCFe  = zeros(nsteps,1);
    
    v0     = v;
    pt0    = pt;
    pf0    = pf;
    
    for ils = 1:nsteps
        v   = v0 + alphav(ils).*dv;
        pt  = pt0+ alphav(ils).*dpt;
        pf  = pf0+ alphav(ils).*dpf;
        
        Pt  = reshape(pt       ,[nx  ,ny  ]);
        Pf  = reshape(pf       ,[nx  ,ny  ]);
        Vx  = reshape(v(1:nVx) ,[nx+1,ny  ]);
        Vy  = reshape(v(nVx+1:end),[nx  ,ny+1]);
        % NonLinearity 1
        divV  = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
        dphidt = (1-phi).*((1-CN)*divV + CN*divVo );
        % NonLinearity 2
        phi     = phi_o + dt*dphidt;
        kmuf    = kmuf0.*(phi/phi0).^n;
        eta_phi = mus./phi*eta2mus .* (1+0.5*(1./R -1).*(1+tanh((Pf-Pt)./lamp)));
        k_U = zeros(size(Pf));  k_U(:,2:end  ) = (kmuf(:,1:end-1) + kmuf(:,2:end))/2;
        k_B = zeros(size(Pf));  k_B(:,1:end-1) = (kmuf(:,1:end-1) + kmuf(:,2:end))/2;
        
        kgy = ((rhofg - rogbg)*(k_B - k_U))/dy; % gravity
        %%%%%%
        rog   = rhofg.*phi + (1-phi ).*rhosg - rogbg;
        rogy  = zeros(size(Vy));  rogy(:,2:end-1) = 0.5*(rog(:,1:end-1)+rog(:,2:end));
        % Momentum residual
        Vx_e        = zeros(nx+1,ny+2); Vx_e(:,2:end-1) = Vx;
        Vy_e        = zeros(nx+2,ny+1); Vy_e(2:end-1,:) = Vy;
        Vx_e(:,  1) = Vx(:,  1); % Free slip
        Vx_e(:,end) = Vx(:,end); % Free slip
        Vy_e(  1,:) = Vy(  1,:); % Free slip
        Vy_e(end,:) = Vy(end,:); % Free slip
        tau_xx      = 2*mus.*(diff(Vx,1,1)/dx - 1/2*divV);
        tau_yy      = 2*mus.*(diff(Vy,1,2)/dy - 1/2*divV);
        tau_xy      =  musv.*(diff(Vx_e,1,2)/dy + diff(Vy_e,1,1)/dx );
        Res_x       = diff(-Pt + tau_xx,1,1)/dx + diff(tau_xy(2:end-1,:),1,2)/dy;
        Res_y       = diff(-Pt + tau_yy,1,2)/dy + diff(tau_xy(:,2:end-1),1,1)/dx - rogy(:,2:end-1);
        RMe         = [Res_x(:) ; Res_y(:)];
        nRMe(ils)   = norm(RMe)/((nx+1)*ny + (ny+1)*nx);
        % Total continuity residual
        RCTe        = -(diff(Vx,1,1)/dx+diff(Vy,1,2)/dy) - (Pt-Pf)./(eta_phi.*(1-phi));
        nRCTe(ils)  = norm(RCTe(:))/(nx*ny);
        % Fluid continuity residual
        qDx = zeros(nx+1,ny); qDy = zeros(nx,ny+1); % zero fluxes
        qDx(2:end-1,:) = -0.5*(kmuf(1:end-1,:)+kmuf(2:end,:)).*diff(Pf,1,1)/dx;
        qDx([1 end],:) = [qDx_L';qDx_R'];
        qDy(:,2:end-1) = -0.5*(kmuf(:,1:end-1)+kmuf(:,2:end)).*diff(Pf,1,2)/dy;
        qDy(:,[1 end]) = [qDy_U qDy_B];
        divqD          = diff(qDx,1,1)/dx + diff(qDy,1,2)/dy - kgy;
        RCFe           = -divqD + (Pt-Pf)./(eta_phi.*(1-phi));
        nRCFe(ils)     = norm(RCFe(:))/(nx*ny);
    end
    nRMe = [nRMe0; nRMe]; nRCTe = [nRCTe0; nRCTe]; nRCFe = [nRCFe0; nRCFe]; alphav = [0 alphav];
    % [~,ibest] =  min(nRMe);
    % [~,ibest] = min(nRCTe);
    [~,ibest] = min(nRCFe);
    if ibest==1
        LSsuccess=0; 
        alpha=0; 
        if noisy>=1 
            fprintf('No descent found - continuing ...'); 
        end
        nRCFe(1) = 2*max(nRCFe); [~,ibest] = min(nRCFe);
        LSsuccess=1; 
        alpha = alphav(ibest);
    else
        LSsuccess=1; 
        alpha = alphav(ibest); 
    end
    
    if noisy>=1, fprintf(' LS: Selected alpha = %2.2f\n', alpha );
        fprintf(' LS: NonLin res. ||res.v ||=%2.4e, ||res.v /res.v0 ||=%2.4e\n', nRMe (ibest), nRMe (ibest)/resnlv0  );
        fprintf(' LS: NonLin res. ||res.pt||=%2.4e, ||res.pt/res.pt0||=%2.4e\n', nRCTe(ibest), nRCTe(ibest)/resnlpt0 );
        fprintf(' LS: NonLin res. ||res.pf||=%2.4e, ||res.pf/res.pf0||=%2.4e\n', nRCFe(ibest), nRCFe(ibest)/resnlpf0 );
    end
end
% 