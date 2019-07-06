% Backward difference method to solve heat equation
% with Neumann conditions

%{  
    du(x,t)/dt = D*d^2u(x,t)/dx^2
    
    a < x < b
    t1 < t < t2

    Boundary Conditions

    u(x,0) = f(x) 
    du(a,t)/dx = 0
    du(b,t)/dx = 0

INPUT: a,b
       t1,t2
       D -> Diffusion coefficient
       M -> number of steps in x direction
       N -> number of steps in t direction

%}

% w = heatbdneumann(D,a,b,t1,t2,M,N) 
% We are going to use heatbdneumann(1,0,1,0,1,20,20)
% N = 20 undiconditionally stable

function w = heatbdneumann(D,xl,xr,yb,yt,M,N)

%modify to use
f = @(x) sin(2*pi*x).^2; 

h = (xr-xl)/M; k = (yt-yb)/N; 
m = M+1 ; % add two values ( dirichlet is M-1 )
n = N;

lambda = D*k/(h*h);
a = diag(1+2*lambda*ones(m,1)) + diag(-lambda*ones(m-1,1),1);
a = a + diag(-lambda*ones(m-1,1),-1);  % define A matrix

% Neumann conditions
% Thanks to two extra values (M+1) we can put this coefficients
a(1,:) = [-3 4 -1 zeros(1,m-3)]; 
a(m,:) = [zeros(1,m-3) -1 4 -3];

w(:,1) = f(xl+(0:M)*h)'; % we need to know w0 and wM, two extra values

for j =1:n
    b = w(:,j); 
    b(1) = 0; % thanks to approximate of first derivate 
    b(m) = 0; % and the homogeneus conditions both are zeros.
    w(:,j+1) = a\b; % w/a
end

x = (0:M)*h;
t = (0:n)*k;

mesh(x,t,w')
view(60,30);
axis([xl xr yb yt -1 1]);



