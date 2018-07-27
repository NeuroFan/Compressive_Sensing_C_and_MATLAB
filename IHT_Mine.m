function x=IHT_Mine(y,A,k,iteration)
[M,N]=size(A);
xnew=zeros(N,1);
r=y-A*xnew;
for i=1:iteration
   g=A'*(r);
   w=0.001;%(g'*g)/(g'*(A'*A)*g);
   x=xnew+w*g;
   [v,s]=sort(abs(x),'descend');
   T=s(1:k);
   xnew(:)=0; 
   xnew(T)=x(T);
   r=y-A*xnew;
end
x=xnew;
end
