function J = fdjacob(x,y,W_i,ubrk)

% Cecile Cabanes, 2019 (CC) changes of LMA.m &fdjacob.m to get similar convergence and
% solutions as lsqnonlin.m with the option 'levenberg-marquardt'
% (optimization toolbox function)

epsilon = sqrt(eps);  % CC
%epsilon=0.1;       %??????
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