% We are going to use heatfd(0,1,0,1,10,250)

function w = heatfd(xl,xr,yb,yt,M,N)

f = @(x) sin(2*pi*x).^2;
r = @(t) 0*t;
l = @(t) 0*t;

D = 1;
h = (xl-xr)/M; k = (yt-yb)/N; 
m = M-1 ; n = N;

lambda = D*k/(h*h);
a = diag(1-2*lambda*ones(m,1)) + diag(lambda*ones(m-1,1),1);
a = a + diag(lambda*ones(m-1,1),-1);

lside = l(yb+(0:n)*k);
rside = l(yb+(0:n)*k);
w(:,1) = f(xl+(1:m)*h)';

for j =1:n
    w(:,j+1) = a*w(:,j)+lambda*[lside(j);zeros(m-2,1);rside(j)];
end

w = [lside;w;rside];
x = (0:m+1)*h;
t = (0:n)*k;

mesh(x,t,w')
view(60,30);
%axis([xl xr yb yt -1 1]);



