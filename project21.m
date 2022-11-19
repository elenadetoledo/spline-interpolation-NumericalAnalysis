%SECTION 6.4, COMPUTER PROJECT 7

%Insert the values for the points. I´ll use a 2x1 matrix. Each line will
%identify a point and each column one axis. First column will be the x axis
%and the second one will be the y axis.


%I´ll use A for Sx and E for Sy

A = zeros(11,2);
A(1,1) = 1.95;
A(2,1) = 2.5;
A(3,1) = 2.55;
A(4,1) = 2.8;
A(5,1) = 3;
A(6,1) = 3.10;
A(7,1) = 3.45;
A(8,1) = 3.52;
A(9,1) = 3.8;
A(10,1) = 4.3;
A(11,1) = 5.45;

A(1,2) = 1.61;
A(2,2) = 5.75;
A(3,2) = 2.11;
A(4,2) = 4.75;
A(5,2) = 2.10;
A(6,2) = 3.51;
A(7,2) = 4.90;
A(8,2) = 1.60;
A(9,2) = 5.90;
A(10,2) = 0.75;
A(11,2) = 1.10;

E=zeros(11,2);

E(1,1)=0.75;
E(2,1)=1.1;
E(3,1)=1.6;
E(4,1)=1.61;
E(5,1)=2.1;
E(6,1)=2.11;
E(7,1)=3.51;
E(8,1)=4.75;
E(9,1)=4.9;
E(10,1)=5.75;
E(11,1)=5.9;


E(1,2)=4.3;
E(2,2)=5.45;
E(3,2)=3.52;
E(4,2)=1.95;
E(5,2)=3;
E(6,2)=2.55;
E(7,2)= 3.10;
E(8,2)=2.8;
E(9,2)=3.45;
E(10,2)=2.5;
E(11,2)=3.8;


%First, I´ll plot the data points
t= tiledlayout(1,1);
title(t,'Script letter: L')

for (i = 1:11)
   hold on
   plot(A(i,1),A(i,2),'.--k')
end


n = 11 -1 ;

%We need to store the distance between nodes. I'll use a vector called h
%for this purpose

h = zeros(1,n);
for i = 1:n
    h(i) = A(i+1,1)-A(i,1);
end

hy = zeros(1,n);
for i = 1:n
    hy(i) = E(i+1,1)-E(i,1);
end

mh = zeros(1,n-2) %need to insert the hs from h0 to hn-2
for i = 1:n-2
    mh(i) = h(i+1);
end

mhy = zeros(1,n-2) %need to insert the hys from hy0 to hyn-2
for i = 1:n-2
    mhy(i) = hy(i+1);
end

u = zeros(1,n-1);
for i = 1:n-1 %ui = 2(hi + hi-1)
    u(i) = 2*(h(i+1)+h(i));
end

uy = zeros(1,n-1);
for i = 1:n-1 %ui = 2(hi + hi-1)
    uy(i) = 2*(hy(i+1)+hy(i));
end

mmatrix = diag(u)+diag(mh,1)+diag(mh,-1)
mmatrixy = diag(uy)+diag(mhy,1)+diag(mhy,-1)

b = zeros(1,n);
for i = 1:n
    b(i) = (6/h(i))*(A(i+1,2)-A(i,2));
end

by = zeros(1,n);
for i = 1:n
    by(i) = (6/hy(i))*(E(i+1,2)-E(i,2));
end

v = zeros(1,n-1)
for i = 1:n-1
    v(i)= b(i+1)-b(i)
end

vy = zeros(1,n-1)
for i = 1:n-1
    vy(i)= by(i+1)-by(i)
end


vtrans = transpose(v);
vytrans = transpose(vy);
%Now that the system has been created, following the naming convention seen
%in the book, I'll use the function linsolve to solve the system of
%equations
solutions = linsolve(mmatrix,vtrans);
solutionsy =linsolve(mmatrixy,vytrans);

solutionsZ = transpose(solutions);
solutionsZy = transpose(solutionsy);

%Insert 0's in the positions z0 and zn:

z = zeros(1,11);
for i = 2:10
    z(i) = solutions(i-1)
end

zy= zeros(1,11);
for i = 2:10
    zy(i) = solutionsy(i-1)
end



D = zeros(1,n); %A was already used above to denote the matrix that contains the interpolation points
B = zeros(1,n);
C = zeros(1,n);

F = zeros(1,n); %E was already used above to denote the matrix that contains the interpolation points
G = zeros(1,n);
K = zeros(1,n);


for i=1:10
    
    D(i) = (1/(6*h(i)))*(z(i+1)-z(i));
    B(i) = z(i)/2;
    C(i) = (-h(i)/6)*z(i+1)-(h(i)/3)*z(i)+(1/h(i))*(A(i+1,2)-A(i,2));
end

for i=1:10
    
    F(i) = (1/(6*hy(i)))*(zy(i+1)-zy(i));
    G(i) = zy(i)/2;
    K(i) = (-hy(i)/6)*zy(i+1)-(hy(i)/3)*zy(i)+(1/hy(i))*(E(i+1,2)-E(i,2));
end
%Now that we have the coefficients for the splines, we will create and plot
%each spline


for i= 1:10
    syms f(x) %Plot for each spline
	f(x) = A(i,2)+(x-A(i,1))*(C(i)+(x-A(i,1))*(B(i)+(x-A(i,1))*D(i)));
	x=linspace(A(i,1),A(i+1,1));
	plot(x,f(x),'b');
    hold on
end

for i= 1:10
    syms g(x) %Plot for each spline
	g(x) = E(i,2)+(x-E(i,1))*(K(i)+(x-E(i,1))*(G(i)+(x-E(i,1))*F(i)));
	x=linspace(E(i,1),E(i+1,1));
	plot(g(x),x,'g');
    hold on
end

    
hold off


