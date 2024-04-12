clear;clc
format compact 
P=[1 -2 -15 20 44 -48];
Q=[1 3 -15 -19 30];
m = length(P) - 1
n = length(Q) - 1 
k = 1
Syl = zeros(m+n-k+1,m+n-2*k+2);

for i=1:n-k+1
    Syl(i:i+m,i)=P;
end

for j=1:m-k+1
    Syl(j:j+n,i+j)=Q;
end


[~,H,~] = Ksvd(Syl)

H = diag(H)

r = sum(H < 1e-10)

Syl = zeros(m+n-r+1,m+n-2*r+2)

for i=1:n-r+1
    Syl(i:i+m,i)=P;
end

for j=1:m-r+1
    Syl(j:j+n,i+j)=Q;
end

[U,H,V]=svd(Syl);
gcd_cofactor = V(:, end);

w=gcd_cofactor(1:n-r+1)
v=gcd_cofactor(n-r+2:end)



c=length(w)-1
Con = zeros(c+r+1,r+1)

for i=1:r+1
    Con(i:i+c,i)=w;
end

Q=Q'
u=Con\Q
u=u/u(1)

