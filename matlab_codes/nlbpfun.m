
function residual = nlbpfun (ubrk_i)

% function residual = nlbpfun (ubrk)
%
% Ned Campion, September 2005
% Breck Owens, rewritten Dec 2006 to include some fixed break points

global A breaks nbr1 ubrk_g
global xf yf W_i xblim

debug = 0;

     % The passed arguments are assumed in the following order:
     % lstest = leastsq ('llbpfun', [1 1 1], options, '', x, y)
     % lsqnonlin(@nlbpfun, b0, lb, ub)
     
     % use ubrk rather than breaks for convergence reasons

% if some break points are fixed, load up fixed ones and then add in
% the ones that we are trying to fit
if nbr1 > 1
    ubrk(1:nbr1-1) = ubrk_g(1:nbr1-1);
    ubrk(nbr1:length(ubrk_g)) = ubrk_i;
else
    ubrk = ubrk_i;
end

mb = length (ubrk);
fnumer = zeros (size (ubrk) );
fnumer(1) =  exp (ubrk(1));
fdenom = 0.0;
for i = 2:mb 
    fnumer(i) = fnumer(i-1) + exp (ubrk(i));
end  % for numerator computation.
fdenom = 1.0 + fnumer(mb);

ftem =  (xblim(2) - xblim(1) ) / fdenom;
      % If you want to watch the convergence, remove the ';' on breaks = .
breaks = xblim(1) + ftem.*fnumer;

if debug
    disp(['nlbpfun: ubrk = ' num2str(ubrk)])
    disp(['nlbpfun: breaks = ' num2str(breaks)])
end


if ~isempty(find(diff(breaks)==0))
    jj=find(diff(breaks)==0);
    breaks(jj+1)=breaks(jj+1) + 0.00001;
end



% get piece-wise linear parameters and residuals
[A, residual] = brk_pt_fit (xf, yf, W_i, breaks);

return
