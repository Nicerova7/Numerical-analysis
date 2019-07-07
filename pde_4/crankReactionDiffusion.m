% Crank-Nicolson method to solve reaction-diffusion equation

%{  
    du(x,t)/dt = D*d^2u(x,t)/dx^2 + C*u(x,t)
    
    a < x < b
    t1 < t < t2

    Boundary Conditions

    u(x,0) = f(x) 
    u(a,t) = l(t)
    u(b,t) = r(t)

INPUT: a,b
       t1,t2
       D -> Diffusion coefficient
       M -> number of steps in x direction
       N -> number of steps in t direction

%}

% w = crankReactionDiffusion(C,D,a,b,t1,t2,M,N) 
% We are going to use crankReactionDiffusion(9.5,1,0,1,0,1,20,20)
% N = 20 undiconditionally stable

function w = crankReactionDiffusion(C,D,xl,xr,yb,yt,M,N)

%modify to use
f = @(x) sin(pi*x).^2; % L = 1
r = @(t) 0*t;
l = @(t) 0*t;


h = (xr-xl)/M; k = (yt-yb)/N; 
m = M-1 ; n = N;

lambda = D*k/(h*h);

% define A matrix
a = diag(2-k*C+2*lambda*ones(m,1)) + diag(-lambda*ones(m-1,1),1); % less k*C 
a = a + diag(-lambda*ones(m-1,1),-1); 

% define B matrix
b = diag(2+k*C-2*lambda*ones(m,1)) + diag(lambda*ones(m-1,1),1); % add k*C
b = b + diag(lambda*ones(m-1,1),-1); 

lside = l(yb+(0:n)*k); % define left side
rside = r(yb+(0:n)*k); % define rigth side
w(:,1) = f(xl+(1:m)*h)'; % initial conditions

for j =1:n
    sides = [lside(j) + lside(j+1); zeros(m-2,1); rside(j)+ rside(j+1)];
    w(:,j+1) = a\(b*w(:,j)+lambda*sides); % w/a
end

w = [lside;w;rside];
x = xl+(0:M)*h;
t = yb+(0:N)*k;

mesh(x,t,w')
view(60,30);
axis([xl xr yb yt 0 1.5]);



