%This program takes an m-by-n matrix A and computes the SVD decomposition.
%U and V are m-by-m and n-by-n orthonormal matrices such that A=U*SIGMA*V'
%for m-by-n zero matrix SIGMA with diagonal consisting of the singular values of
%A.


function [U,SIGMA,V] = SVDtrunc(A,p)
    [m,n] = size(A);
    dual = 0;   %indicate whether working with A or A'
    if m > n    %standardizes m>n to m<n
        A = A';
        [m,n] = size(A);
        dual = 1;
    end
    [H,Q] = Hessenberg(A'*A);
    [SIGMA2,G] = IterateGivensTrunc(H,50,p); %iterations determine accuracy
    SIGMA = sqrt(abs(SIGMA2));
    SIGMAINV = zeros(n);
    for i = 1:n %construct pseudoinverse of SIGMA
        if SIGMA(i,i) ~= 0
            SIGMAINV(i,i) = 1/SIGMA(i,i);
        end
    end
    V = Q'*G';
    U = A*V*SIGMAINV;
    U = U(1:m,1:m); %trim U to be m-by-m
    SIGMA = SIGMA(1:m,1:n); %trim SIGMA to be m-by-n
    if dual == 1 %convert answer from A' to A if swapped
        [U,V] = deal(V,U);
    end
end