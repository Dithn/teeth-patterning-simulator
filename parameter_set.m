function [parameters,boundary_type,dirichlet_val] = parameter_set(simulation_number)
%   PARAMETER_SET - provides parameter values for three types of simulations
%                   1) sequential patterning
%                   2) wildtype
%                   3) mutant
%               param_list - vector of all model parameters
%               boundary type - matrix 5x2 with 0-1 values; 0-homogeneous Neumann, 1-Dirichlet; first column is for x=0, second column is for x=R; rows correspond sequentially to I,A,C,M,S 
%               dirichlet_val - matrix 5x2 with values in case of Dirichlet boundary type
  

  %% ==================== Part 1: Boundary conditions ====================
  % parameters_list:
  %        1. diff_I  - diffusivity of inhibitor              
  %        2. diff_A  - diffusivity of activator
  %        3. diff_C  - diffusivity of cells
  %        4. sigma_t - auto-inhibition of inhibitor for unsteady solution (Turing instability, mature zone)
  %        5. sigma_b - auto-inhibition of inhibitor for bistable solution (immature zone)
  %        6. mu - auto activation/inhibition of activator
  %        7. m_sigma_star - maturation threshold changing sigma between sigma_touring and sigma_bist
  %        8. Q_low - separation of kinetics of activator in three intervals
  %        9. Q_high
  %       10. lambda_g - rate of domain growth
  %       11. alpha_n - rate of spontaneous maturation
  %       12. alpha_a - rate of activator enhanced maturation
  %       13. alpha_m_star - activator threshold for enhanced maturation
  %       14. tau - delay time of enhanced maturation
  %       15. ks - rate of activator production by mesenchyme signal
  %       16. S_low - low mesenchyme signal threshold
  %       17. S_high - high mesenchyme signal threshold
  %       18. m_dist_star - maturation threshold for distal activation
  %       19. gamma - rate of production of mesenchyme signal
  %       20. chi - chemo-sensitivity of cells
  %       21. c_max - maximal cells density
  %       22. c_star - cells density threshold for inhibitor degradation (>0 in table)
  %       23. c_mass - total mass of cells
  %       24. alpha - rate of inhibitor degradation by cells
  %       25. L    - initial length of a domain
  %       26. T    - final time of simulation


  parameters_matx  = [ 1      1      1      1    1    1;           
		       0.01   0.01   0.01   0.01 0.01 0.01;
		       1      1      1      0.02 0.02 0.02;
		       0.55   0.55   0.4    0.55 0.55 0.55;
		       20     20     6      0    0    0;
		       0.4    0.4    0.4    0.4  0.4  0.4;
		       0.7    0.7    0.7    0.7  0.7  0.7;
		       0.7    0.7    0.7    0.7  0.7  0.7;
		       0.7    1      1      0.7  0.7  0.7;
		       0.005  0.005  0.005  0    0    0;
		       0.0005 0.0005 0.0005 0    0    0;
		       0.1    0.1    0.1    0    0    0; 
		       0.9    0.9    0.9    0    0    0;
		       100    100    100    0    0    0;
		       10     10     10     0    0    0;
                       1.1    1.1    1.1    0    0    0;
                       1.4    1.4    1.4    0    0    0;
		       0.1    0.1    0.1    0    0    0;
		       0.001  0.001  0.001  0    0    0;
		       1      1      1      0.1  0.5  0.5;
		       1      1      1      3    3    3;
		       1      1      1      2.9  2.9  2.9;
                       1      1      1      3	 3    3;	       
		       1      1      1      0.3  0.3  0.3;
		       1      1      1      2.6  2.6  3
		       800    800    800    150  150  150];

  parameters = parameters_matx(:,simulation_number);
  
  %% ==================== Part 2: Boundary conditions ====================
  %  Boundary conditions are equal for all simulations
  if ismember(simulation_number,[1 2 3])
    boundary_type = [0 0; 1 0; 0 0; 1 1; 1 0];
    dirichlet_val = [0 0; -1 0; 0 0; 1 0; 0 0];
  else
    boundary_type = zeros(5,2);
    dirichlet_val = zeros(5,2);
  end

  

  
