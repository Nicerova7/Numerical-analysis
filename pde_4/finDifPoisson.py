# Finite differences method to solve Poisson equation.

"""
 d^2u(x,y)/dx^2  + d^2u(x,y)/dy^2 = f(x,y)  

 a <= x <= b
 c <= y <= d

 Boundary Conditions:
 
 u(x,y) = g(x,y)   If x = a   or   x = b  and  c <= y <= d
  
 and

 u(x,y) = g(x,y)   If y = c   or   y = d and   a <= x <= b

INPUT: a, b, c, d
        m,n >= 3
        tol -> tolerance
        N   -> Number of iterations. 

OUTPU: w_i,j approximation to u(xi,yi) to all i = 1, ... , n-1 and
                                       to all j = 1, ... , m-1 

       Or a message that the number of iterations was exceeded

"""
