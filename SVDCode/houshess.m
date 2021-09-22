function [H,Q]=houshess(A)
% [H,Q]=HOUSHESS(A) computes the matrices H and Q such that H=Q¡¯AQ.
[n,m]=size(A);
if n~=m 
    error('Only square matrices'); 
end
Q=eye(n); 
H=A;
for k=1:n-2
[v,beta]=vhouse(H(k+1:n,k)); 
I=eye(k);  
N=zeros(k,n-k); 
m=length(v);
R=eye(m)-beta*(v*v.');
H(k+1:n,k:n)=R*H(k+1:n,k:n);
H(1:n,k+1:n)=H(1:n,k+1:n)*R; 
P=[I, N; N.', R]; 
Q=Q*P;
end 
end
