%This function is for iterating the Givens function. H is the resulting
%Upper Hessenberg matrix, and Q is the composition of Givens rotations

function [H,Q] = IterateGivensTrunc(A,k,m)
    [H,Q] = GivensTrunc(A,m);
    for i = 1:k-1
        [H,G] = GivensTrunc(H,m);
        Q = G*Q;
    end
end