function Y = Matrix2Vec(X, n)
% Transform a matrix into a column vector
%
% Copyright (C) 2009 Jaeyoung Lee (jyounglee@yonsei.ac.kr)

   y = 0;
   N = n*(n+1)/2;
   p = zeros(N,1);
   for j = 1:n
     y = y + j-1;
     for i = 1:j
         p(y+i) = X(i,j);
     end
   end
   Y=p;
