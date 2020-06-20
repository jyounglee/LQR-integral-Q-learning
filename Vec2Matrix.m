function Y = Vec2Matrix(X, n)
% Transform a vector into a symmetric matrix
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)

   N = n*(n+1)/2;
   A = zeros(Y,1);
   for j = 1:n
     y = y + j-1;
     for i = 1:j
         A(i,j) = X(y+i);
     end
   end
   Y = A + A' - diag(diag(A));
