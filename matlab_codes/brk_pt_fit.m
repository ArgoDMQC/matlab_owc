function [A, residual] = brk_pt_fit (x, y, W_i, b);

% function [A, fit] = brk_pt_fit(x, y, W_i, b);
%
% routine to get least-squares estimates for piecewise linear fit
% with break points at prescribed points fits
%
%  y = A(1) + A(2)*x+eps                     for x(1) <= x <= b(1)
%  y = A(1) + A(2)*b(1) + A(3)*(x-b(1))+eps  for b(1) <= x <= b(2)
%  ...                                            ...
%  y = A(1) + A(2)*(b(2)-b(1)) + ...
%                         A(m+2)*(x-b(m)+eps for b(m) <= x <= x(n)
%  where x = vector of observed indepedent variables [n]
%        y = vector of observed dependent variables  [n]
%        W_i = inverse of weights for fit [n]
%              if W_i is empty, then equal weighting for each point
%        b = vector of break points as values of x   [m]
%
%        A = fitting coefficients                    [m+1]
%        yg = fitted estimate of y
%
% Tne inverse of the weighting matrix is passed to this routine through a
% global variable W_i
% It returnes the E matrix relating the observations to the fit parameters
% of the linear fits through the global variable E
%
% Ned Campion, September 2005
% Breck Owens, June 2006

% pass E back to fit_cond for debugging plots

%
if nargin < 4
    % no break point specified, do linear fit
    b = [];
end

m = length(b);
n = length(x);
n1 = length(y);


% check to see if break points are included
if n1 ~= n
    fit = 999;
    A = zeros(m+2,1);
    disp(' *** brk_pt_fit: input vectors not consistent.');
    return
end
% set up vectors to make sure they are the same
x = reshape(x,n,1);
y = reshape(y,n,1);

     % add  1st point as break point, to make indexing pretty
btem = [x(1) b ];
      
     % form matrix
     % include an intercept as well as the trends between each break 
     % point
E = zeros(n,m+2);
E(:,1) = ones(n,1);
ixb = sorter (btem, x);
for j = 1:m+1
    ib = find (ixb == j);% pointer to x values greater than break point j
    E(ib,j+1) = x(ib) - btem(j);
    ii = find (ixb > j); % pointer for break points less than the one just
                       % below the x value
    if (~isempty (ii) )
        E(ii,j+1) = btem(j+1) - btem(j);
    end  % if lower triangular.
end  % for coefficients.

     % get least squares estimate of the A's
if ~isempty(W_i)
    B = E'*W_i*E;
else
    B = E'*E;
end

if det(B) == 0
    fit = 999;
    A = zeros(m+2,1);
    residual = y;
    disp('Error in brk_pt_fit determinant of matrix is 0')
    return
end


% calculate fit parameters
if ~isempty(W_i)
    A = B \ E'*W_i*y;
else
    A = B \ E'*y;
end
% calculate fit estimates
residual = y - E*A;

return;
 
function pointer = sorter(msites, sites)
 
[ignored,index] = sort([msites(:).' sites(:).']);
pointer = find(index>length(msites))-(1:length(sites));

