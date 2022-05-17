function x=cst(phi,n)

if(mod(n,2)==0)
    x=1./tan(phi);
else
    x=1./sin(phi);
end
end