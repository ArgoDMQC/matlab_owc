
function [ubrk, resnorm, residual] = LMA(ubrk)


global A breaks nbr1 ubrk_g
global xf yf W_i xblim

nbr = length(ubrk);
increment = 1e-10;
lambda =  diag(increment*ones(1,nbr));
e1=1; 
max_cyc = 50;
cycles = 0;
while e1>1e-15 
    %Jacobian of least squares fit for P
    J = fdjacob(xf,yf,W_i,ubrk);
    residual = nlbpfun(ubrk);
    previous_residual = residual;
    cntr=0;
    while sum(previous_residual.^2)<=sum(residual.^2)
        previous_residual = residual;
        N = J'*J + lambda;
        deltaU = inv(N)*-J'*residual;
        %Test Size of U's and distance between breaks
        new_ubrk = ubrk+deltaU';
        mb = length (new_ubrk);
        fnumer = zeros (size (new_ubrk) );
        fnumer(1) =  exp (new_ubrk(1));
        fdenom = 0.0;
        for i = 2:mb 
            fnumer(i) = fnumer(i-1) + exp (new_ubrk(i));
        end  % for numerator computation.
        fdenom = 1.0 + fnumer(mb);
        ftem =  (xblim(2) - xblim(1) ) / fdenom;
        deltas = ftem*exp(new_ubrk);
        %If distance between breaks is too small stop process
        % and do not except further changes to U's.  Also stop if
        % lamda get's so large matlab gives Inf values
        if ~isempty(find(deltas<0.001)) || ~isempty(find(~isfinite(deltaU)))
            break
        end

        residual = nlbpfun(ubrk+deltaU');
        lambda = lambda*10;
        cntr=cntr+1;
    end
    if ~isempty(find(deltas<0.001)) || ~isempty(find(~isfinite(deltaU)))
        break
    end
    lambda = lambda/10;
    e1 = norm(deltaU);
    ubrk = ubrk + deltaU';
    cycles = cycles +1;
    if cycles > 50
        disp(['Iteration count has exceeded preset maximum of ' ...
            num2str(max_cyc) ' cycles'])
        break
    end
end

resnorm = sum(residual.^2);


