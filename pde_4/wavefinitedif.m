% Finite differences to solve wave equation

%{
    du/dtt = c^2 * du/dxx    ; u(x,t) ; wave equation

    a <= x <= b
    0  <=  t

    Initial and boundary conditions
    u(x,0) = f(x)
    du(x,0)/dt = g(x)
      
    u(a,t) = l(t)  
    u(b,t) = r(t)   

%}

% we use wavefinitedif(2,0,1,0,1,20,40)

function w = wavefinitedif(c,xl,xr,yb,yt,M,N)

%modify to use
f = @(x) sin(pi*x);
g = @(x) 0*x;
l = @(x) 0*x;
r = @(x) 0*x;

h = (xr-xl)/M ; k = (yt-yb)/N;
m = M-1; n = N;

lambda = c*k/h;

% define A matrix
a = diag(2 - 2*(lambda^2)*ones(m,1)) + diag((lambda^2)*ones(m-1,1),1);
a = a + diag((lambda^2)*ones(m-1,1),-1);

lside = l(yb + (0:n)*k);
rside = r(yb + (0:n)*k);


gg = g(yb + (1:m)*h)';

w(:,1) = f(xl+(1:m)*h)';
w(:,2) = 1/2*a*w(:,1) + k*gg + (1/2)*(lambda^2)*[lside(1);zeros(m-2,1);rside(1)];

for j = 2:n
   w(:,j+1) = a*w(:,j) - w(:,j-1) + (lambda^2)*[lside(j);zeros(m-2,1);rside(j)];
end

w = [lside;w;rside];
x = (0:m+1)*h;
t = (0:n)*k;

mesh(x,t,w')
view(60,30);
axis([xl xr yb yt -1 2]);