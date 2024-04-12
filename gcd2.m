clear;clc
syms x
k = 1;
p = (x-4)^(k)*(x-5)^(2*k)*(x-6)^(3*k)
p = coeffs(p, x);
p = double(fliplr(p));
p=p/norm(p,2)
pp=p';
q = polyder(p);  
    
m = length(p) - 1; 
n = length(q) -1 ;
Syl = zeros(m+n-k+1,m+n-2*k+2);

for i=1:n-k+1
    Syl(i:i+m,i)=p;
end

for j=1:m-k+1
    Syl(j:j+n,i+j)=q;
end
Syl
[H]=svd(Syl);
r = sum(H < 1e-10);

Syl = zeros(m+n-r+1,m+n-2*r+2);

for i=1:n-r+1
    Syl(i:i+m,i)=p;
end

for j=1:m-r+1
    Syl(j:j+n,i+j)=q;
end

[~,H,V] = svd(Syl);


%%
%[H] = svd(Syl)
%r2 = sum(H < 1e-10)

%% find u and v
gcd_cofactor = V(:, end);

w = gcd_cofactor(1:n-r+1);
v = gcd_cofactor(n-r+2:end);

c = length(w)-1;
x = length(v)-1;

W = zeros(c+r+1,r+1);

for i=1:r+1
    W(i:i+c,i)=w;
end

q=q';
u=W\q;
roots(u)
z=length(u)-1;

U=zeros(z+x+1,x+1);
for i=1:x+1
    U(i:i+z,i)=u;
end

pp;
v=U\pp;
roots(v)


%%
p=u;
q = polyder(p);

m = length(p) - 1;
n = length(q) -1 ;
Syl = zeros(m+n-k+1,m+n-2*k+2);

for i=1:n-k+1
    Syl(i:i+m,i)=p;
end

for j=1:m-k+1
    Syl(j:j+n,i+j)=q;
end

[H]=svd(Syl);
r = sum(H < 1e-10);

Syl = zeros(m+n-r+1,m+n-2*r+2);

for i=1:n-r+1
    Syl(i:i+m,i)=p;
end

for j=1:m-r+1
    Syl(j:j+n,i+j)=q;
end

[~,H,V] = svd(Syl);

gcd_cofactor = V(:, end);

w = gcd_cofactor(1:n-r+1);
v = gcd_cofactor(n-r+2:end);

c = length(w)-1;
x = length(v)-1;

W = zeros(c+r+1,r+1);

for i=1:r+1
    W(i:i+c,i)=w;
end

q=q';
u=W\q;
roots(u)
z=length(u)-1;

U=zeros(z+x+1,x+1);
for i=1:x+1
    U(i:i+z,i)=u;
end


v=U\p;
roots(v)

%%
p=u;
q = polyder(p);

m = length(p) - 1;
n = length(q) -1 ;
Syl = zeros(m+n-k+1,m+n-2*k+2);

for i=1:n-k+1
    Syl(i:i+m,i)=p;
end

for j=1:m-k+1
    Syl(j:j+n,i+j)=q;
end

[H]=svd(Syl);
r = sum(H < 1e-10);

Syl = zeros(m+n-r+1,m+n-2*r+2);

for i=1:n-r+1
    Syl(i:i+m,i)=p;
end

for j=1:m-r+1
    Syl(j:j+n,i+j)=q;
end

[~,H,V] = svd(Syl);

gcd_cofactor = V(:, end);

w = gcd_cofactor(1:n-r+1);
v = gcd_cofactor(n-r+2:end);

c = length(w)-1;
x = length(v)-1;

W = zeros(c+r+1,r+1);

for i=1:r+1
    W(i:i+c,i)=w;
end

q=q';
u=W\q;
roots(u)
z=length(u)-1;

U=zeros(z+x+1,x+1);
for i=1:x+1
    U(i:i+z,i)=u;
end


v=U\p;
roots(v)

%%
p=u;
q = polyder(p);

m = length(p) - 1;
n = length(q) -1 ;
Syl = zeros(m+n-k+1,m+n-2*k+2);

for i=1:n-k+1
    Syl(i:i+m,i)=p;
end

for j=1:m-k+1
    Syl(j:j+n,i+j)=q;
end

[H]=svd(Syl);
r = sum(H < 1e-10);

Syl = zeros(m+n-r+1,m+n-2*r+2);

for i=1:n-r+1
    Syl(i:i+m,i)=p;
end

for j=1:m-r+1
    Syl(j:j+n,i+j)=q;
end

[~,H,V] = svd(Syl);

gcd_cofactor = V(:, end);

w = gcd_cofactor(1:n-r+1);
v = gcd_cofactor(n-r+2:end);

c = length(w)-1;
x = length(v)-1;

W = zeros(c+r+1,r+1);

for i=1:r+1
    W(i:i+c,i)=w;
end

q=q';
u=W\q;
roots(u)
z=length(u)-1;

U=zeros(z+x+1,x+1);
for i=1:x+1
    U(i:i+z,i)=u;
end


v=U\p;
roots(v)

