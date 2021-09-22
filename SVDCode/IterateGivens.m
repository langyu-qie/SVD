%This function is for iterating the Givens function. H is the resulting
%Upper Hessenberg matrix, and Q is the composition of Givens rotations

function [H,Q] = IterateGivens(A,k)
    [n,~] = size(A);
    if (k==0)
        H=A;
        Q=eye(n);
    else   
        [H,Q] = Givens(A);
        for i = 1:k-1
            [H,G] = Givens(H);
            Q = G*Q;
        end
    end
end