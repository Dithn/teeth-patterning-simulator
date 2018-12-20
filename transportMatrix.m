function [Tm] = transportMatrix(dx,Nx,R,phiR,diff,boundary_type,reaction_vec,simulation_number)

  if ismember(simulation_number,[4 5 6])
    Tm = (1/dx^2)*(spdiags(ones(Nx,1),-1,sparse(Nx,Nx))+spdiags(-2*ones(Nx,1),0,sparse(Nx,Nx))+spdiags(ones(Nx,1),1,sparse(Nx,Nx)));
    Tm(1,1) = -1/dx^2; Tm(Nx,Nx) = -1/dx^2;
  else
    ll = (R^2/diff)*(phiR(Nx-1)+phiR(Nx))/2;
    l0 = (R^2/diff)*(phiR(1:Nx-1) + phiR(2:Nx))/2;
    lr = (R^2/diff)*(phiR(1)+phiR(2))/2;
    diag_down =  [l0; 0]./expm1([l0; 0]*dx);
    diag_cent = -[0;l0(2:Nx-1);0]./expm1([0;l0(2:Nx-1);0]*dx) - [0;l0(1:Nx-2);0].*exp([0;l0(1:Nx-2);0]*dx)./expm1([0;l0(1:Nx-2);0]*dx);
    diag_up   =  [0; l0].*exp([0; l0]*dx)./expm1([0; l0]*dx);  
    if boundary_type(1) == 1
      diag_cent(1) = 0;     diag_up(2) = 0;
      reaction_vec(1) = 0;
    else
      diag_cent(1) =  - lr/expm1(lr*dx);
      diag_up(2)   = lr*exp(lr*dx)/expm1(lr*dx);
    end
    if boundary_type(2) == 1
      diag_cent(Nx) = 0;    diag_down(Nx-1) = 0;
      reaction_vec(Nx) = 0;
    else
      diag_cent(Nx)   = (1+dx/2)*phiR(Nx) - ll*exp(ll*dx)/expm1(ll*dx);
      diag_down(Nx-1) = ll/expm1(ll*dx);
    end
    Tm = (1/(R*R*dx))*spdiags(diag_down,-1,sparse(Nx,Nx))+spdiags(diag_cent,0,sparse(Nx,Nx))+spdiags(diag_up,1,sparse(Nx,Nx));
  end
end

