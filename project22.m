%SECTION 6.4, COMPUTER PROJECT 8
%Interpret the results of the following numerical experiment and draw some
%conclusions.

%a. Define p to be the polynomial of degree 20 that interpolates the
%function f(x)=(1+6x^2)^-1 at 21 equally spaced nodes in the interval
%[-1,1]. Include the endpoints as nodes. Print a table of f(x), p(x) and
%f(x)-p(x) at 41 equally spaced points in the interval.

%b.Repeat the experiment using the Chebyshev nodes given by xi =
%cos((i-1)pi/20) for i[1,21]

%c.With 21 equally spaced knots, repeat the experiment using a cubic
%interpolating spline.

%PART A
%This is a lagrange interpolation.
lint = abs(1--1); %length of the interval of interpolation
h = lint/21; %Compute the distance between equispaced nodes

init = -1;
xs = zeros(1,21);
y = zeros(1,21);
coefl = ones(1,21);

%Insert the nodes in the vector nodes. Equispaced. Starting in -1.
for i = 1:21
    xs(i)= init;
    init= init +h;
end

%Compute the coefficients of each of the polynomials li(x)
for i = 1:21
    for j = 1:21
        if i ~= j
            coefl(i)= coefl(i)*1/(xs(i)-xs(j));
        end
    end
end

%Compute the image at each of the nodes
for i = 1:21
    y(i) = 1/(1+6*xs(i)^2);
end


syms p(x)
p(x) =0;

syms l(x)
l(x) = 1;

for j = 1:21
    for i = 1:21
        if i~=j
            l(x)=(x-xs(i))*l(x); %we get lj(x)
        end
    end
    p(x) = p(x) + y(j)*l(x)*coefl(j);
    l(x) = 1;
end

h1 = 2/41;
count = 1;
valf = zeros(1,41);
valp = zeros(1,41);
valdif = zeros(1,41);

for i = -1:h1:1
    valf(count) = 1/(1+6*i^2);
    valp(count) = p(i);
    valdif(count) = valf(count)-valp(count);
    count = count +1;
end

fprintf('For lagrange');
fprintf('\nf(x)\t\tp(x)\t\tf(x)-p(x)')
for i = 1:41
    fprintf('\n%f\t\t%f\t\t%f',valf(i),valp(i),valdif(i))
end



%PART B

%The chebyshev nodes are given: xi = cos[(i-1)pi/20]

nodes = zeros(1,21);
for i = 1:21
    nodes(i)=cos((i-1)*pi/20);
end


%Compute the coefficients of each of the polynomials li(x)
for i = 1:21
    for j = 1:21
        if i ~= j
            coefl(i)= coefl(i)*1/(nodes(i)-nodes(j));
        end
    end
end

%Compute the image at each of the nodes
for i = 1:21
    y(i) = 1/(1+6*nodes(i)^2);
end


syms p(x)
p(x) =0;

syms l(x)
l(x) = 1;

for j = 1:21
    for i = 1:21
        if i~=j
            l(x)=(x-xs(i))*l(x); %we get lj(x)
        end
    end
    p(x) = p(x) + y(j)*l(x)*coefl(j);
    l(x) = 1;
end

h1 = 2/41;
count = 1;
valf = zeros(1,41);
valp = zeros(1,41);
valdif = zeros(1,41);

for i = -1:h1:1
    valf(count) = 1/(1+6*i^2);
    valp(count) = p(i);
    valdif(count) = valf(count)-valp(count);
    count = count +1;
end

fprintf('\nFor Chebyshev');
fprintf('\nf(x)\t\tp(x)\t\tf(x)-p(x)')
for i = 1:41
    fprintf('\n%f\t\t%f\t\t%f',valf(i),valp(i),valdif(i))
end
%PART C

%Lo comento para que vaya mas rapido, pero en principio funciona bien


for i = 1:21
   hold on
   plot(xs(i),y(i),'.--k')
end

u = zeros(1,19);
for i = 1:19
    u(i) = 2*(h+h);
end
vecth = zeros(1,18);
for i = 1:18
    vecth(i)=h;
end
b = zeros(1,20);
for i = 1:20
    b(i) = (6/h)*(y(i+1)-y(i));
end

mmatrix = diag(u)+diag(vecth,1)+diag(vecth,-1);
v = zeros(1,19);
for i = 1:19
    v(i)= b(i+1)-b(i);
end

vtrans = transpose(v);
solutions = linsolve(mmatrix,vtrans);
z = zeros(1,21);
for i = 2:20
    z(i) = solutions(i-1);
end
A = zeros(1,20);
B = zeros(1,20);
C = zeros(1,20);

for i=1:20
    
    A(i) = (1/(6*h))*(z(i+1)-z(i));
    B(i) = z(i)/2;
    C(i) = (-h/6)*z(i+1)-(h/3)*z(i)+(1/h)*(y(i+1)-y(i));
end

for i= 1:20
    syms g(x) %Plot for each spline
	g(x) = y(i)+(x-xs(i))*(C(i)+(x-xs(i))*(B(i)+(x-xs(i))*A(i)));
	x=linspace(xs(i),xs(i+1));
	plot(x,g(x),'b');
    hold on
end
hold off

h1 = 2/41;
count = 1;
valf = zeros(1,41);
valp = zeros(1,41);
valdif = zeros(1,41);
for i = -1:h1:1
    valf(count) = 1/(1+6*i^2);
    valp(count) = g(i);
    valdif(count) = valf(count)-valp(count);
    count = count +1;
end

fprintf('\nFor cubic splines');
fprintf('\nf(x)\t\tp(x)\t\tf(x)-p(x)')
for i = 1:41
    fprintf('\n%f\t\t%f\t\t%f',valf(i),valp(i),valdif(i))
end


