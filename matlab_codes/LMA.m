
function [ubrk, resnorm, residual] = LMA(ubrk)


global A breaks nbr1 ubrk_g
global xf yf W_i xblim
% Cecile Cabanes, 2019 (CC) changes of LMA.m & fdjacob.m to get similar convergence and
% solutions as lsqnonlin.m with the option 'levenberg-marquardt' (matlab
% optimization toolbox function)

nbr = length(ubrk);
%increment = 1e-10;
increment = 1e-2;    % CC: change initial value of increment
lambda =  diag(increment*ones(1,nbr));
e1=1;
max_cyc = 50;
cycles = 0;
done = false;
tolFun = 1e-6;
tolX = 1e-6;
while ~done
    %while e1>1e-15 % CC: new convergence criteria
    % Jacobian of least squares fit for P
    J = fdjacob(xf,yf,W_i,ubrk);
    residual = nlbpfun(ubrk);
    previous_residual = residual;
    cntr=0;
    while sum(previous_residual.^2)<=sum(residual.^2)
        % previous_residual = residual; % CC: only lambda is updated in the inner while loop
        N = J'*J + lambda;
        % deltaU = inv(N)*-J'*residual;
        deltaU = inv(N)*-J'*previous_residual; % CC: only lambda is updated in the inner while loop
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
            %disp(['LMA.m : Iteration stopped because distance between breaks is to small'])
            break
        end
        residual = nlbpfun(ubrk+deltaU');
        
        if sum(previous_residual.^2)>sum(residual.^2)&&cntr==0 % CC : if the condition is met, we exit the loop before multiplying lambda by 10
            break
        end
        lambda = lambda*10;
        cntr=cntr+1;
        
        if ~isempty(find(~isfinite(lambda))) % Massimo Paciaroni: stop the algorithm if lambda get's very large
            %disp(['LMA.m : Iteration stopped because lambda get''s infinite'])
            break
        end
        
    end
    if ~isempty(find(deltas<0.001)) || ~isempty(find(~isfinite(deltaU)))||~isempty(find(~isfinite(lambda)))
        break
    end
    lambda = lambda/10;
    e1 = norm(deltaU);
    ubrk = ubrk + deltaU';
    
    
    % CC: new convergence criteria
    if e1< tolX*(sqrt(eps)+norm(ubrk))
        done=true;
    end
    if abs(sum(residual.^2) - sum(previous_residual.^2)) <= tolFun*sum(previous_residual.^2)
        done=true;
    end
    
    cycles = cycles+1;
    
    if cycles > 50
        disp(['LMA.m : Iteration count has exceeded preset maximum of ' ...
            num2str(max_cyc) ' cycles'])
        break
    end
    
end

resnorm = sum(residual.^2);


