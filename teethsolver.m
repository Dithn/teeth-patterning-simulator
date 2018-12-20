%% TEETH PATTERNING SIMULATOR
% Author:      Monika Twarogowska ; mtwarogowska@gmail.com
% Date:        2018/12/19
% Description: Numerical solution of a teeth patterning model introduced in the article "Modeling Edar expression reveals the hidden dynamics of tooth signaling center patterning" - Alexa Sadier, Monika Twarogowska, Klara Steklikova, Luke Hayden, Anne Lambert, Pascal Schneider, Vincent Laudet, Maria Hovorakova, Vincent Calvez and Sophie Pantalacci.


%  TEETHSOLVER Compute time evolution for a chosen simulation_type of the model with the following variables: 
%                I - inhibitor
%                A - activator
%                M - tissue maturation level
%                S - mesenchyme inhibitors
%                C - cells density
%
%              System of reaction-diffusion equations:
%
%              I_t = DI(t)*I_xx + x*psi(t)*I_x + F(I,A,M) - alpha*C*(C>C_threshold)
%              A_t = DA(t)*A_xx + x*psi(t)*A_x + G(I,A)
%              C_t = DC(t)*c_xx + x*psi(t)*C_x - chi(t)*(C*(1-C/Cmax)*A_x)_x
%              M_t = x*psi(t)*M_x + H1(A,M)
%              S_t = x*psi(t)*S_x + H2(A,S)
%
%              where:
%              * time dependent parameters (due to domain growth scaling) are given by
%               
%              DI(t)  = diffI/R(t)^2; DA = diffA/R(t)^2; DC = diffC/R(t)^2; diffI,diffA,diffC>0
%              chi(t) = chi0/R(t)^2; chi0>0 
%              psi(t) = (dR/dt)/R(t); R(t) = R0 + growth_rate * t; growth_rate>0
%
%              * reaction functions are given by
%              F(I,A)  = -sigma(M) * I + A
%                        | -I - mu*(A + 2*Q_low)     if        A < -Q_low
%              G(I,A)  = | -I + mu*A                 if   -Q_low <= A <= Q_high
%                        | -I - mu*(A - 2*Q_high)    if        Q_high < A
%
%              H1(A,M) = m_rate_natural*(1 - M) + m_rate_enhanced*M*(1-M)*(A_delay>maturation_A_threshold)
%              H2(A,S) = 
%
%              * boundary confitions are set to
%              inhibitor        :       Neumann homegeneous at x = 0 and x = R
%              activator        :       Dirichlet at x = 0 and Neumann at x = R and Neumann homogeneous in case of fusion simulation; 
%              cells density    :       Neumann homegeneous at x = 0 and x = R
%              maturation       :       Dirichlet at x = R;
%              mesenchyme       :       Neumann at x = R
%

clear ; close all; clc

fprintf('Time evlution of teeth patterning. Simulation types: \n');
fprintf('1. Sequential patterning without traveling wave \n');
fprintf('2. Wildtype \n');
fprintf('3. Mutant \n');
fprintf('4. Fusion: small domain & weak chemotaxis \n');
fprintf('5. Fusion: small domain & strong chemotaxis \n');
fprintf('6. Fusion: large domain & strong chemotaxis \n');

prompt = 'Chose simulation number (integer from 1 to 6):';
simulation_number = input(prompt);
while ~ismember(simulation_number,[1 2 3 4 5 6])
  prompt = 'Wrong simulation number. Chose an integer between 1 and 6:';
  simulation_number = input(prompt);
end

%% ==================== Part 1: Parameters and Boundary conditions ====================
[parameters, boundary_type, dirichlet] = parameter_set(simulation_number);

%% Model parameters
diff_I = parameters(1);
diff_A = parameters(2);
diff_C = parameters(3);

sigma_t      = parameters(4);
sigma_b      = parameters(5);
mu           = parameters(6);
m_sigma_star = parameters(7);
Q_low        = parameters(8);
Q_high       = parameters(9);
lambda_g     = parameters(10);
alpha_n      = parameters(11);
alpha_a      = parameters(12);
alpha_m_star = parameters(13);
tau          = parameters(14);
ks           = parameters(15);
s_low        = parameters(16);
s_high       = parameters(17);
m_dist_star  = parameters(18);
gamma        = parameters(19);
chi          = parameters(20);
c_max        = parameters(21);
c_star       = parameters(22);
c_mass       = parameters(23);
alpha        = parameters(24);
L            = parameters(25);
T            = parameters(26);

%% Numerical parameters
% Space
Nx = L*100;
xL = 0;
xR = L;
dx = (xR-xL)/(Nx-1);
x  = xL:dx:xR;
				% time
theta = 0.5;  % IMEX scheme for diffusion
if theta >0.5
  dt = 0.5*dx^2;   
else
  dt = 0.5*dx;
end

Nt   = floor(T/dt) + 1;
Ntau = floor(tau/dt) + 1;


%% ==================== Part 2: Variables initializations  ====================

dtP = 10;                              % iteration interval at which plot is visualized
T_plot = 0:dtP:T;

				% Initial values
if ismember(simulation_number,[1 2 3])
  I  = zeros(Nx,2);  I(:,1)  = 0;                                    I0 = I(:,1);
  A  = zeros(Nx,Nt); A(:,1)  = -1+(2/atan(25))*atan(25*x');          A0 = A(:,1);
  C  = zeros(Nx,2);  C(:,1)  = ones(Nx,1);                           C0 = C(:,1);
  M  = zeros(Nx,2);  M(:,1)  = (pi/2+atan(-15*x'+5))/(pi/2+atan(5)); M0 = M(:,1);
  S  = zeros(Nx,2);  S(:,1)  = s_low + 0.05;                         S0 = S(:,1);
else
  I  = zeros(Nx,2);  I(:,1)  = 0;                                    I0 = I(:,1);
  A  = zeros(Nx,Nt); A(:,1)  = 1+0.1*sin(2*pi*x'/L).^2;              A0 = A(:,1);
  C  = zeros(Nx,2);  C(:,1)  = (c_mass/L)*ones(Nx,1);                C0 = C(:,1);
  M  = zeros(Nx,2);  M(:,1)  = ones(Nx,1);                           M0 = M(:,1);
  S  = zeros(Nx,2);  S(:,1)  = zeros(Nx,1);                          S0 = S(:,1); 
end

Rt = zeros(1,Nt);  Rt(1,1) = xR;
DT = zeros(1,Nt);  DT(1,1) = dt;

dtS = 1;
T_save = 0:dtS:T; Nt_save = length(T_save);
IA = zeros(Nx,Nt_save); IA(:,1) = I0;   
AA = zeros(Nx,Nt_save); AA(:,1) = A0;
CA = zeros(Nx,Nt_save); CA(:,1) = C0;
MA = zeros(Nx,Nt_save); MA(:,1) = M0;
SA = zeros(Nx,Nt_save); SA(:,1) = S0;
RA = zeros(1,Nt_save);  RA(1,1) = xR;

k_plot = 1; k_save = 1;
t = 0; ti = 1;

%% ----------------   Part 2.1: Plot ------------------ 
figure1  = figure('units','normalized','outerposition',[0 0 1 1],'Name','Time evolution on a growing domain','NumberTitle','off');
if ismember(simulation_number,[1 2 3])
  plot_grow(x,L,I0,A0,M0,S0,figure1)
else
  plot_chemo(x,L,I0,A0,C0,figure1)
end
%% ===============================================================
%% ==================== Part 3: Main program  ====================
%% ===============================================================

while t<T

  t  = t+dt;
  ti = ti+1;

  %% ----------------   Part 3.1: Domain growth ------------------ 
  R       = xR + lambda_g*(t);
  phiR    = x'*lambda_g/R;
  Rt(ti)  = R;
  sigma   = sigma_t*(M(:,1) > m_sigma_star) + sigma_b*(M(:,1) <= m_sigma_star); 

  %% ----------------   Part 3.2: Chemoattractant  ------------------ 
  if ismember(simulation_number,[4 5 6])
    qR = (1-0.5 * ([C(2:Nx,1);C(Nx,1)]+C(:,1))/c_max);
    qL = (1-0.5 * (C(:,1)+[C(1,1); C(1:Nx-1,1)])/c_max);
    
    L0 = - exp((chi/diff_C)*qR.*([A(2:Nx,ti-1);dirichlet(2,2)]-A(:,ti-1))/2) - exp(-(chi/diff_C)*qL.*(A(:,ti-1)-[dirichlet(2,1);A(1:Nx-1,ti-1)])/2);
    Lm = [exp((chi/diff_C)*qL(2:Nx).*(A(2:Nx,ti-1)-A(1:Nx-1,ti-1))/2);0];
    Lp = [0;exp(-(chi/diff_C)*qR(1:Nx-1).*(A(2:Nx,ti-1)-A(1:Nx-1,ti-1))/2)];
    L0(1)  = - exp((chi/diff_C)*qR(1).*(A(2,ti-1)-A(1,ti-1))/2);
    L0(Nx) = - exp(-(chi/diff_C)*qL(Nx).*(A(Nx,ti-1)-A(Nx-1,ti-1))/2);
    DD = (diff_C*dt/dx^2)*(spdiags(Lm,-1,sparse(Nx,Nx)) + spdiags(L0,0,sparse(Nx,Nx)) + spdiags(Lp,1,sparse(Nx,Nx)));
    C(:,2) = (speye(Nx,Nx)-(1-theta)*DD) \ ((speye(Nx,Nx)+theta*DD)*C(:,1));    
  else
    C(:,2) = zeros(Nx,1);
  end
  %% ----------------   Part 3.3: Inhibitor ------------------ 
  F  = - sigma.* I(:,1) + A(:,ti-1) - alpha*C(:,2).*(C(:,2)>c_star); 
  flow_matrix = transportMatrix(dx,Nx,R,phiR,diff_I,boundary_type(1,:),F,simulation_number);
  %y=full(flow_matrix);
  %y(1:10,1:10)
  %pause
  I(:,2) = (speye(Nx,Nx) - dt*(1-theta)*diff_I*flow_matrix) \ ((speye(Nx,Nx)+dt*theta*diff_I*flow_matrix)*I(:,1) + dt*F);
  
  %% ----------------   Part 3.4: Activator ------------------ 
  G = - I(:,1) - mu*(A(:,ti-1) + 2*Q_low).*(A(:,ti-1) < -Q_low)...
      + mu*A(:,ti-1).*(A(:,ti-1) >= -Q_low).*(A(:,ti-1) <= Q_high)...
      - mu*(A(:,ti-1) - 2*Q_high).*(A(:,ti-1) > Q_high);
  Y = ks*(-(S(:,1)>s_high)+(S(:,1)<s_low)).*(M(:,1)<m_dist_star);
  react = G+Y;
  flow_matrix = transportMatrix(dx,Nx,R,phiR,diff_A,boundary_type(2,:),react,simulation_number);
  A(:,ti) = (speye(Nx,Nx)-dt*(1-theta)*diff_A*flow_matrix) \ ((speye(Nx,Nx)+dt*theta*diff_A*flow_matrix)*A(:,ti-1)+dt*react);
  A_delayed_vector = activatorTimeDelayed(x,Rt,ti,Ntau,A);
  
  %% ----------------   Part 3.5: Maturation and mesenchyme ------------------   
  if ismember(simulation_number,[1 2 3])
    m_law  = alpha_n*(dirichlet(4,1) - M(:,1)) + alpha_a*(A_delayed_vector>alpha_m_star).*M(:,1).*(dirichlet(4,1)-M(:,1));
    M(:,2) = M(:,1) + (phiR*dt/dx).*([M(2:Nx,1); dirichlet(4,2)]-M(:,1)) + dt*m_law;

    s_law  = gamma*A(:,ti).*(M(:,1)<m_dist_star);
    S(:,2) = S(:,1) + (phiR*dt/dx).*([S(2:Nx,1); S(Nx,1)*(boundary_type(5,2)==0)+dirichlet(5,2)*(boundary_type(5,2)==1)]-S(:,1)) + dt*s_law;
  else
    M(:,2) = ones(Nx,1);
    S(:,2) = zeros(Nx,1);
  end
  %% ----------------  Part 3.7: Post-settings ------------------ 

   if ismember(simulation_number,[1 2 3])
     I(:,1) = I(:,2);
     C(:,1) = C(:,2);
     M(:,1) = M(:,2);
     S(:,1) = S(:,2);
   else
     I(:,1) = 0.5*(I(:,2)+I(Nx:-1:1,2)); 
     A(:,ti) = 0.5*(A(:,ti)+A(Nx:-1:1,ti));
     C(:,1) = 0.5*(C(:,2)+C(Nx:-1:1,2));
   end
  
  %% PLOT
  if t > T_plot(k_plot)
    k_plot = k_plot + 1;
    if ismember(simulation_number,[1 2 3])
      plot_grow(x,R,I(:,2),A(:,ti),M(:,2),S(:,2),figure1);
    else
      plot_chemo(x,L,I(:,2),A(:,ti),C(:,2),figure1);
    end
  end
  
  %% SAVE
  if t > T_save(k_save)
    k_save = k_save + 1;
    IA(:,k_save) = I(:,2);
    AA(:,k_save) = A(:,ti);
    CA(:,k_save) = C(:,2);
    MA(:,k_save) = M(:,2);
    SA(:,k_save) = S(:,2);
    RA(1,k_save) = R;
  end
end








