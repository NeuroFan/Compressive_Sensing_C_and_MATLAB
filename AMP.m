function xhat=AMP(y,H,landa,iterAMP,m,n)
    r=y;
    s=zeros(n,1);
    sqM=sqrt(m);

for t=1:iterAMP
    
    maxs=sort(abs(r),'descend');
    sigma=mean(maxs(1:20));%norm(r)/sqM;    
    g=H'*r;
    s=s+g;
    
    %minOFmaxs=mean(maxs(1:10));
    for i=1:n
           s(i)=sign(s(i))*max(abs(s(i))-sigma,0);
    end
    %for i=1:n %thresholding
     %   a=middle(i);
      %  s(i)=sign(a)*max(abs(a)-sigma,0);
    %end
     b=sum(abs(s)>0)/m;
     r=y-H*s+b.*r;
end
    xhat=s;
end  