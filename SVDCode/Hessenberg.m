%This function takes an n-by-n matrix A and converts it to upper Hessenberg
%form H using Householder reflections. Q is the composition of Householder
%reflections such that H=QAQ'.

function [H,Q] = Hessenberg(A)
	[n,~] = size(A);
    for j = 1:n-2
        R = Householder(A(j+1:n,j));
        P(:,:,j) = blkdiag(eye(j),R);
        A = P(:,:,j)*A*P(:,:,j)';
    end
    Q = P(:,:,1);
    for j = 2:n-2
        Q= P(:,:,j)*Q;
    end
    H = A;
end