%This function takes an n-by-n Upper Hessenberg matrix A and uses Givens
%rotations to transform it into another Hessenberg matrix via orthonormal
%similarity transformations. In particular, the function iteratively
%multiplies by a Givens rotation on the left for each column (except the
%last), and then at the end multiplies on the right by the transpose of the
%product of Givens matrices. The end product of Givens matrices is called
%Q, and the output upper Hessenberg Matrix is called H.

function [H,Q] = GivensTrunc(A,m)
    [n,~] = size(A);
    m = min(n,m);
    for j = 1:m-1
        alpha = sqrt(A(j,j)^2 + A(j+1,j)^2);
        c = A(j,j)/alpha;
        s = A(j+1,j)/alpha;
        G(:,:,j) = eye(n);
        G(j,j,j) = c;
        G(j+1,j+1,j) = c;
        G(j,j+1,j) = s;
        G(j+1,j,j) = -s;
        A = G(:,:,j)*A;
    end
    Q = G(:,:,1);
    for j = 2:m-1
        Q = G(:,:,j)*Q;
    end
    H=A*Q';
end