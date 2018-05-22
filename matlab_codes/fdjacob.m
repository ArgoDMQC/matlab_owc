function J = fdjacob(x,y,W_i,ubrk)

epsilon = 0.1;  %??????

m=length(x);
n=length(ubrk);

residual = nlbpfun(ubrk);

for i=1:n
    
    ubrk1 = ubrk;
    ubrk1(i) = ubrk1(i) + epsilon;
    residual1 = nlbpfun(ubrk1);
    
    df = (residual1 - residual)/epsilon;
    
    J(:,i) = df;

end