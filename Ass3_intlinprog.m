clear all

% Object function
C = [10, 0, 20, 11, 12, 7, 9, 20, 0, 14, 16, 18]; 

% Supply constraint
Ain = blkdiag(ones(1,4), ones(1,4), ones(1,4));
supply = [25; 55; 35];
% S2 = 60, Q4
supply2 = [25; 60; 35];

% Demand constraint
Aeq = [eye(4), eye(4), eye(4)];
demand = [15; 45; 30; 25];


%% Interior point
options = optimoptions('intlinprog', 'Display', 'iter');

  x = intlinprog(C,12, Ain, supply, Aeq, demand, zeros(12,1), [], options);
  
  Ain*x;
  Aeq*x;
  
  
  x2 = intlinprog(C,12, Ain, supply2, Aeq, demand, zeros(12,1), [], options);
  
  Ain*x2;
  Aeq*x2;
  
%% Dual simplex
  options = optimoptions('intlinprog', 'Display', 'iter');

  x3 = intlinprog(C,12, Ain, supply, Aeq, demand, zeros(12,1), [], options);
  
  Ain*x;
  Aeq*x;
  
  % S2 = 60, Q4
  supply = [25; 60; 35];
  x4 = intlinprog(C,12, Ain, supply2, Aeq, demand, zeros(12,1), [], options);
  
  Ain*x2;
  Aeq*x2;
  
  %% Results
  X = [x, x3, x2, x4]
  
  Used = Ain*X