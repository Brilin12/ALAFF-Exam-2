function [ A ] = Create_Poisson_problem_A( N )

  % Create the archtypical matrix A for an N x N Poisson problem (5-point
  % stencil.
  A = zeros(N^2, N^2);

  for i = 1:N^2
      % Set the diagonal
      A(i,i) = 4;

      % Set the entries of the first sub and super diagonals
      if not (mod(i,N) == 0)
          A(i,i+1) = -1;
      end
      if not (mod(i-1,N) == 0)
          A(i,i-1) = -1;
      end

      % Set the other off-diagonal entries
      if (i+N <= N^2)
          A(i,i+N) = -1;
      end
      if (i-N > 0)
          A(i,i-N) = -1;
      end
  end
  
end





