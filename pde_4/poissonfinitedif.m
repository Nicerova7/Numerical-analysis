% Finite difference to solve Poisson equation

%{
    du/dxx + du/dyy = f(x,y);  Poisson eq (if f(x,y) = 0 will be Laplace eq)
    
    xl <= x <= xr
    yb <= y <= yt 
    

    Initial and Boundary conditions
    u(x,yb) = g1(x)
    u(x,yt) = g2(x)
    u(xl,y) = g3(y)
    u(xr,y) = g4(y)

%}

% We will use w = poissonfinitedif(0,1,1,2,4,4)

function w = poissonfinitedif(xl,xr,yb,yt,M,N)

% modify to use
f = @(x,y) 0;  
g1 = @(x) sin(pi*x);
g2 = @(x) sin(pi*x);
g3 = @(y) 0;
g4 = @(y) 0;

m = M+1; n = N+1;  mn = m*n;

h = (xr-xl)/M;  h2 = h^2;
k = (yt-yb)/N;  k2 = k^2;

% Grid values
x = xl + (0:M)*h;
y = yb + (0:N)*k;

A = zeros(mn,mn);  b = zeros(mn,1);

for i = 2:m-1
    for j = 2:n-1
        A(i + (j-1)*m,i-1 + (j-1)*m) = 1/h2;
        A(i + (j-1)*m,i+1 + (j-1)*m) = 1/h2;
        
        A(i + (j-1)*m,i + (j-1)*m) = -2/h2 - 2/h2;
        
        A(i + (j-1)*m,i + (j-2)*m) = 1/k2;
        A(i + (j-1)*m,i + j*m) = 1/k2;
        
        b(i + (j-1)*m) = f(x(i),y(j));
    end
end

for i = 1:m
    j = 1; A(i + (j-1)*m,i + (j-1)*m) = 1;
           b(i + (j-1)*m) = g1(x(i));
    j = n; A(i + (j-1)*m,i + (j-1)*m) = 1;
           b(i + (j-1)*m) = g2(x(i));
end

for j = 2:n-1
    i = 1; A(i + (j-1)*m,i + (j-1)*m) = 1;
           b(i + (j-1)*m) = g3(y(j));
    i = m; A(i + (j-1)*m,i + (j-1)*m) = 1;
           b(i + (j-1)*m) = g4(y(j));
end

v = A\b;  % b/A
w = reshape(v(1:mn),m,n);
mesh(x,y,w')



