[root_estimates,multiplicity_structure]=GCDRoot(P,theta,tol,phi)

%%
if ~exist('theta', 'var')
    theta = 10^-8;
end

if ~exist('tol', 'var')
    tol = 10^-10;
end

if ~exist('phi','var')
    phi = 100;
end

%%
u0 = P;
s = length(P);

%計算u_{m-1}=f
u(1,1:s) = P;
m = length(P)-1;

for j = 2:m+1    
    u0 = u(j-1,1:s); 
    u00 = polyder(u0);

    len_diff = length(u0) - length(u00);
    if len_diff > 0
        u00 = [zeros(1, len_diff), u00];
    elseif len_diff < 0
        u0 = [zeros(1, -len_diff), u0];
    end

    u(j,1:s) = gcd(u0, u00);
end

 u_m_1=u(m,1:s);


%%
p = tol * norm(u_m_1,2);

for m = 1:s
    for l=1:p
    [varsigma,x] = inverseIteration(u_m_1,l);
        if varsigma < theta * norm(u_m_1)
            [u,v,w] = GCDsystem(u_m_1);
            [u0,v0,w0] =initialterate(x,u_m_1);
                
            [H]=Ksvd(v0);

            Conv_uv = convmtx(u',m)*v';
            Conv_uw = convmtx(u',m)*w';
            residual = norm([Conv_uv;Conv_uw] - [u_m_1; u_m_1]);
        end
    end
    tol = max(tol,phi*residual);
    dm = length(H)-1;
end

%%
k=d1;
L = zeros(1, k);

for j = 1:k
    for t = 1:length(d)
        if d(t) >= k-j+1
            L(j) = t;
            break;
        end
    end
end

for m = 1:s
    roots{m} = roots(dm);
end


root_estimates = roots{m};
multiplicity_structure = L(j);
