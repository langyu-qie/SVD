%This function is for iterating the Givens function. H is the resulting
%Upper Hessenberg matrix, and Q is the composition of Givens rotations

function [H,Q] = IterateGivensGraph(A,k)
    [n,~] = size(A);
    [H,Q] = Givens(A);
    for j = 1:n
        ColumnNormH(j,1) = norm(A(1:n,j));
        ColumnNormH(j,2) = norm(H(1:n,j));
    end
    for i = 1:k-1
        [H,G] = Givens(H);
        Q = G*Q;
        for j = 1:n
            ColumnNormH(j,i+2) = norm(H(1:n,j));
        end
    end
    x = 1:k+1;
    loglog(x,ColumnNormH(1,1:k+1));
    hold on
    for i = 2:n
        loglog(x,ColumnNormH(i,1:k+1));
    end
    legend
    hold off
end