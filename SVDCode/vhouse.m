function [v,beta]=vhouse(x)
n=length(x); 
x=x/norm(x); 
sigma=x(2:n).'*x(2:n); 
v=[1; x(2:n)]; 
if sigma==0
beta=0; 
else
   xnrm=sqrt(x(1)^2+sigma); 
   if x(1)<=0
      v(1)=x(1)-xnrm; 
   else
      v(1)=-sigma/(x(1)+xnrm); 
   end
   beta=2*v(1)^2/(sigma+v(1)^2); 
   v=v/v(1); 
end
end