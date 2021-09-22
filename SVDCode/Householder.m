%This function takes a vector x as an input and outputs the Householder
%reflection R such that Rx has the 2-norm of x in the first component and
%zeros in the others. 

function R = Householder(x)
    n = length(x);
    alpha = sign(x(1))*norm(x);
    x = -x;
    x(1) = x(1) - alpha;
    R = eye(n) - (2/(x.'*x))*x*x.';
end
    
    