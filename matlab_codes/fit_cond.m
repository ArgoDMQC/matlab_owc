
function [xfit, condslope, condslope_err, time_deriv, time_deriv_err, sta_mean, ...
    sta_rms, NDF, fit_coef, fit_breaks] = fit_cond(x,y,n_err,lvcov,varargin)

%-----------------------------------------------------------------------
% [xfit, condslope, condslope_err, time_deriv, time_deriv_err, sta_mean, ...
%    sta_rms, NDF, fit_coef, fit_breaks] = fit_cond(x,y,n_err,lvcov,varargin)
%
%
% INPUT:
% x,y = observations
% n_err = rms a priori error estimate for each data point
% lvcov = matrix expressing vertical and horizontal covariances (SEE WJO)
%
% Additional parameters that can be passed to the routine.
%        The data is passed as a pair of arguments, the name of the
%        variable and its value
%
% initial_breaks    1st guess for break points
%                   [Default equally spaced between 2nd and next to last
%                   points]
% max_no_breaks     maximum number of break points to be evaluated.
%                   Valid range includes -1 and 0, where
%                   -1 = mean offset only
%                    0 = mean offset and linear trend
%                   [Default = 4];
% number_breaks     specify one or more number of break points to be
%                   fitted
%                    0 = linear fit with no break points
%                   -1 = mean offset
% breaks            evaluate fit for a specified number of break points,
%
%                   if number_breaks > length(breaks), then we are
%                   specifying some breaks, but will compute the
%                   other breaks.
% nloops            specify the number of loops to compute error
%                   of piecewise fit
%                   [Default = 200]
%
% OUTPUT:
% xfit              profile number (unique of x)
% condslope         fit estimate for each profile
% condslope_err     estimated rms error for each profile
% time_deriv        estimated change per profile
% time_deriv_err    estimated rms error in change per profile
% sta_mean          mean difference between estimate and actual values
%                   averaged over each profile
% sta_rms           rms difference between estimates and actual values
%                   averaged over each profile
% NDF               The effective number of independent observations when
%                   the off-diagonal coariance (lvcov) is taken into
%                   account
%
%------------------------------------------------------------------------
%
% Ned Campion, November 2005
% modified by Breck Owens October, 2008
%
%
% To decide which fit is optimal, we will use the small sample variation of the
% Akaike Information Criterion to choose the fit.   Having chosen the
% number of parameters, we then use the F-test to see if the reduction in
% variance relative to the original variance is statistically significant.
%
% We have also used the correlation matrix for the horizontal and vertical
% scales to estimate the number of effective degrees of freedom for the
% fits and for estimating the uncertainties of the fits
% This function implements a non-linear fit of a piecewise linear fit.  The
% methodology is described in:
% Jones, R.H. and I. Dey, 1995, Determining one or more change points.
% Chemistry and Physics of Lipids, 76, 1-6.
%
%
% Cecile Cabanes, 2017: force the fit to an offset only if NDF <13. and display
% a warning : to track change see change config 129

if nargin < 4
    disp('FIT_COND inputs must have at least 4 arguments')
    return
end

if exist('lsqnonlin')==0
disp('Using LMA.m instead of optimization toolbox functions')
end

% variables for the non-liner least squares fitting routine
global A breaks nbr1 ubrk_g
global xf yf W_i xblim

%% default values for possible variables passed to fit_cond
max_brk_dflt = 4; % default maximum number of break points
max_brk_in = [];
max_brk = [];
nbr1 = -1; % 1st break to consider
brk_init = []; % initial guess for the break points
setbreaks = 0;

nloops = 200;  %Number of loops for the profile fit error

%% parameters for lsqnonlin
MaxFunEvals = 1000; % maximum function calls
TolFun = 1e-6; % convergence criteria for lsqnonlin

%% form x and y variable for the fit
% get all profile numbers passed to this routine
xfit = unique(x);
nfit = length(xfit);
% exclude bad points rescale variables to make sure everything is O(1)
% before making the fits
% remove bad points from series
good = find(isfinite(y) & isfinite(x));
x = x(good);
y = y(good);
n_err = n_err(good);
lvcov = lvcov(good,good);
npts = length(x);
if npts == 0 % we failed to find a single good point
    disp('WARNING fit_cond called with no valid data points')
    xp = [];
    condslope      =NaN;
    condslope_err  =NaN;
    time_deriv     =NaN;
    time_deriv_err =NaN;
    sta_mean       =NaN;
    sta_rms        =NaN;
    NDF            =[];
    fit_coef       =[];
    fit_breaks     =[];
    return
end

%% condition the series so that the fit is well behaved
% first sort by the independent variable
[x,i] = sort(x);
y = y(i);
n_err = n_err(i);
lvcov = lvcov(i,i);
x = reshape(x,npts,1);
y = reshape(y,npts,1);
n_err = reshape(n_err,npts,1);
% scale x to go from -1 to 1
x0 = (x(npts)+x(1))/2.;
if x(1) ~= x(npts)
    xscale = (x(npts)-x(1))/2.;
else
    xscale = 1;
end
% remove mean of y and scale y by the standard deviation
y0 = mean(y);
yscale = std(y,1);
if yscale == 0
    yscale = 1;
end
% x and y used for fitting routines
xf = (x-x0)/xscale;
yf = (y-y0)/yscale;
n_err =  n_err/yscale;
xfit = (xfit-x0)/xscale;
% get profile times that we will use as independent variable to get
% error statistics, xp could be different than xfit if there is a profile
% with no good data
[xp,ip] = unique(xf);
nprof = length(xp);

err_var = (n_err).^2; % convert errors from rms to variance

W_i = diag(mean(err_var)./err_var); % weights for weighted Least-squares
% W_i = W_i*lvcov; % include off diagonal terms

% use correlation matrix to compute number of degrees of freedom
NDF = sum(ones(npts,1) ./ (lvcov*ones(npts,1)));

% Residual sum of squares for initial series
RSS0 = sum(((yf).^2) ./err_var);

if length(xp) > 3
    % find 2nd and next to last profile and use them as limits for the break points
    xblim = [xp(2) xp(nprof-1)];
else
    xblim = [1 1]; % should never be used in this case, ie too few profiles
end
%% parse input to find if parameters are to be set
nvarargin = length(varargin);
if rem(nvarargin,2) ~= 0
    disp('FIT_COND - Input must be in the form:')
    disp('x, y, error, parameter, value, ...')
    help fit_cond
    return
end
if nvarargin > 0
    % parse input
    for n=1:nvarargin/2
        parm = varargin(2*(n-1)+1);
        value = varargin(2*n);
        if ~iscellstr(parm)
            disp('FIT_COND - Input must be in the form:')
            disp('flt, parameter, value, ...')
            disp('where parameter is a string')
            help fit_cond
            return
        end
        % find parm in list
        param = lower(char(parm));
        switch param
            case {'initial_breaks'},
                % initial guess for break points included as input
                brk_init = value{:};
                % convert to rescaled units
                brk_init = (brk_init-x0)/xscale;
                brk_init = (brk_init-xblim(1))/diff(xblim);

            case {'max_no_breaks'},
                if ~isempty(value{:})
                    max_brk_in = value{:};
                    nbr1 = -1;
                end
            case {'number_breaks'},
                pbrk = value{:};
                nbr1 = pbrk;
                max_brk_in = pbrk;
            case {'nloops'},
                nloops = value{:};
            case { 'breaks'}
                if ~isempty(value{:})
                    breaks=value{:};
                    breaks = (breaks-x0)/xscale;
                    nbr = length(breaks);
                    nbr0=nbr;  % fix ccabanes 6/10/2020
                    setbreaks=1;
                end
            otherwise,
                disp(['FIT_COND: Parameter ' param ' not found in parameter list'])
        end % end switch
    end
end
% setbreaks is a flag that has been set to indicate that the break points
% have been prescribed


%% initialize variable for search over number of break points
b_pts = ones(max_brk_in,max_brk_in+1)*NaN;% break points, add one to allow for 1st point
b_A = ones(max_brk_in+2,max_brk_in+1)*NaN;% parameters of piece-wise linear fit
RSS = ones(1,max_brk_in+2)*NaN; % residual sum of squares
AIC = ones(1,max_brk_in+2)*NaN; % AICc test to choose optimal fit

if setbreaks % we have set breaks, now check to see range of break points to be tested
    if isempty(max_brk_in) % we have only specified break points
        max_brk_in = nbr;
        nbr1 = nbr;
    elseif max_brk_in > nbr % we have specified some breaks that are to be fixed
        nbr1 = nbr+1; % do fit possible breaks as those specified up to the maximum number of breaks
        % get fit with specified breaks
        [A, residual] = brk_pt_fit (xf, yf, W_i, breaks);
        b_pts(1:nbr,nbr+1) = breaks';
        b_A(1:nbr+2,nbr+2) = A(1:nbr+2);
        RSS(nbr+2) = sum(residual.^2 ./ err_var);
        p = 2*(nbr+1); % number of parameters
        AIC(nbr+2) = NDF*log(RSS(nbr+2)/npts) + NDF*(NDF+p)./(NDF-p-2);
    else
        nbr1 = nbr; % we have specified same numbr of breaks as specified or made an error entering the break points
    end
    max_brk = max_brk_in;
    pbrk = nbr1:max_brk;
else % no break points entered
    if isempty(max_brk_in); % no max break points set
        max_brk_in = max_brk_dflt; % set maximum break points to default
    end
    max_brk = max_brk_in;
    pbrk = nbr1:max_brk;
end

% if too few profiles, limit fit to offset only %    
if nprof < 6
    disp(['WARNING: Only have ' num2str(nprof) ' good profiles, will estimate offset only'])
    pbrk = -1;
end

% if NDF < 13
%    disp(['WARNING: Only have ' num2str(NDF) ' degree of freedom, will estimate offset only'])  % change config 129
%     pbrk = -1;
% end
% if NDF < 2*(max_brk+2)+1   % if NDF is low, AIC criterium will not be valid anymore, there is a maximum number of breakpoints that can be tried
%    if NDF>2*(nbr1+2)+1
%    pbrk =[nbr1:floor((NDF-1)/2 -2)]; 
%    disp(['WARNING: Only have ' num2str(NDF) ' degree of freedom, there is a maximum number of breakpoints that can be tried (' num2str(max(pbrk)) ')'])  % change config 129
%    else
% 	   if setbreaks==1 %
% 	   pbrk = nbr;
% 	   max_brk =nbr;
% 	   nbr1=nbr;
% 	   disp(['WARNING: Only have ' num2str(NDF) ' degree of freedom, will estimate fit with fixed breakpoints only'])  % change config 129
% 	   else
% 	   pbrk = -1;
% 	   disp(['WARNING: Only have ' num2str(NDF) ' degree of freedom, will estimate offset only'])  % change config 129
% 	   end
%    
%    end
% end
if NDF < 2*(max_brk+2)+1 & ~setbreaks  % ccabanes 28/09/2020:  display a warning if NDF is low
     disp([' WARNING : Only have ' num2str(NDF) ' degree of freedom, the fit selected by the software might not be the ''best one''.']);
     disp('  You can change the small spatial and temporal scales in the configuration file, which are used to estimate NDF. You can also define your own breakpoints in the set_calseries.m file.')
     disp(' ')
end


%% Now evaluate range of fits
for nbr = pbrk
    if nbr == -1
        % Offset only, since this is an error weighted average it will not
        % necessarily be zero for yf
        E = ones(npts,1);
        b_A(1,1) = inv(E'*W_i*E) * E'*W_i*yf;
        residual = yf-E*b_A(1,1);
        RSS(1) = sum(residual.^2 ./err_var);
        AIC(1) = NDF*log(RSS(1)/npts)+NDF*(NDF+1)/(NDF-3);
    elseif nbr == 0
        % Linear Fit, no break points
        [A, residual] = brk_pt_fit (xf, yf, W_i);
        b_A(1:2,2) = A(1:2);
        RSS(2) = sum(residual.^2 ./ err_var);
        AIC(2) = NDF*log(RSS(2)/npts)+NDF*(NDF+2)/(NDF-4);
    else
        % includes break points
        nbr2 = length(brk_init);
        if length(brk_init) >= nbr % there are enough initial guesses for break points
            [m n]= size(brk_init);
            if (m>n) brk_init = brk_init'; end
            b_guess = brk_init(1:nbr);
        else
            % give 1st guess for breaks as even distributed between the 2nd and
            % next to last point
            b_guess = -1+2*[1:nbr]/(nbr+1);
        end
        % convert breaks to ubr which insures that the breaks remain in order
        b_g = [-1 b_guess];
        ubrk_g = [];
        for n=1:nbr
            ubrk_g(n) = log((b_g(n+1)-b_g(n))/(1-b_g(nbr+1)));
        end
        if setbreaks
            %if nbr1 == max_brk;  % break points are set get linear-lsq fit to offset and slopes
            if  nbr0 == max_brk % correction ccabanes 06/10/2020
                % break points are static
                [A, residual] = brk_pt_fit (xf, yf, W_i, breaks);
            else % we are supposed to fit over a limited number of breaks
                try
                    [ubrk, resnorm, residual] = lsqnonlin (@nlbpfun, ubrk_g(nbr1:nbr), [],[], ...
                        optimset('DISPLAY','off','MaxFunEvals',MaxFunEvals, 'TolFun',TolFun) );
                    % need to shuffle returned breaks to include ones that are set.
                catch % if error in lsqnonlin get last iteration
                    %#### need to fix the following
                     %[ubrk resnorm residual] = LMA(ubrk_g);
                    [ubrk resnorm residual] = LMA(ubrk_g(nbr1:nbr));  %ccabanes fix
                end
                ubrk = [ubrk_g(1:nbr1-1) ubrk];
            end
        else
            % get non-linear least-squares fit for break points
            try
                [ubrk, resnorm, residual] = lsqnonlin (@nlbpfun, ubrk_g, [],[], ...
                    optimset('DISPLAY','off','MaxFunEvals',MaxFunEvals, 'TolFun',TolFun) );
            catch % if error in lsqnonlin get last iteration
                [ubrk resnorm residual] = LMA(ubrk_g);
            end
        end % end if setbreaks

        b_pts(1:nbr,nbr+1) = breaks';
        b_A(1:nbr+2,nbr+2) = A(1:nbr+2);
        RSS(nbr+2) = sum(residual.^2 ./ err_var);
        p = 2*(nbr+1); % number of parameters
        AIC(nbr+2) = NDF*log(RSS(nbr+2)/npts) + NDF*(NDF+p)./(NDF-p-2);
    end % end if nbr == -1
end % end for nbr

%% select best fit
%if setbreaks & nbr1 == max_brk
if setbreaks && nbr0 == max_brk %fix ccabanes 6/10/2020
    best = pbrk+2;
else%% Decide which fit to use (offset, linear or piecewise fit)
    %    good = find(AIC > 0); % since the maximum likelyhood (residual variance)
    % is scaled by the a priori error at each point, log(RSS/npts) can be
    % negative
    if nbr1 > 1 % we have specified at least one break, but searched over more breaks
        pbrk = (nbr1-1):max_brk; % add in the case for the specified breaks
    end
    good = pbrk+2;
    [AIC_b,best] = min(AIC(good)); % use Akaike Information Criteria for low number of points to choose fit
    best = good(best);
end

if best > 2 % best fit includes break points
    np = (best-1)*2; % number of parameters for fit with break points
else
    np = best;% number of parameters without break point
end

%if setbreaks & nbr1 == max_brk
if setbreaks && nbr0 == max_brk %fix ccabanes 6/10/2020
    comment = ['Fit evaluated using '];
else
    comment = ['Best model found with '];
end
if best>2
    comment = [comment num2str(best-2) ' break points'];
elseif best==2
    comment = [comment 'linear fit'];
else
    comment = [comment 'offset value only'];
end
disp(comment)

% m = best-2;
if best > 2  % there were break points
    breaks = b_pts(1:best-2,best-1)';
else  % a offsset or linear fit was chosen
    breaks = [];
end

%% get paramaeter values for best fit and reconstruct E matrix
A = b_A(1:best,best);

btem = [xf(1) breaks];
E = zeros(npts,best);
E(:,1) = ones(npts,1);
ixb = sorter (btem, xf);
if best > 1
    for j = 1:best-1
        ib = find (ixb == j);% pointer to x values greater than break point j
        E(ib,j+1) = xf(ib) - btem(j);
        ii = find (ixb > j); % pointer for break points less than the one just
        % below the x value
        if (~isempty (ii) )
            E(ii,j+1) = btem(j+1) - btem(j);
        end  % if lower triangular.

    end  % for coefficients.
end

%% reconstruct fit to get parameter and uncertainty estimates

% get uncertainties in fit parameters
B_i = inv(E'*W_i*E);
P = B_i*E'*W_i*diag(err_var)*W_i*E*B_i;
P1 = diag(P); % error variance of fit parameters assuming diagonal error matrix only
P = B_i*E'*W_i*diag(err_var)*lvcov*W_i*E*B_i;
P2 = diag(P); % error variance of fit parameters including off-diagnoal error covariances

% reduce E matrix to have only one value per profile
btem = [xfit(1) breaks];
E = zeros(length(xfit),best);
E(:,1) = ones(length(xfit),1);
ixb = sorter (btem, xfit);
if best>=2
    % recompute A(1) in case xfit(1) is not equal to xf(1) (ccabanes,28/09/2020)
    A(1)=A(1)+A(2)*(xfit(1)-xf(1));
    for j = 1:best-1
        ib = find (ixb == j);% pointer to x values greater than break point j
        E(ib,j+1) = xfit(ib) - btem(j);
        ii = find (ixb > j); % pointer for break points less than the one just
        % below the x value
        if (~isempty (ii) )
            E(ii,j+1) = btem(j+1) - btem(j);
        end  % if lower triangular.

    end  % for coefficients.
end

% fit values
yfit = E*A;

% factor to increase Monte Carlo error estimate taking into account the
% fact that we have correlated noise
P3 = (E*P2)./(E*P1);

% save values of fit while we do Monte Carlo simulation for error estimates
real_A = A; real_E = E; real_breaks = breaks; real_xf = xf; real_yf = yf;

%% Calculate uncertainty in fit
clear ubrk_g
if best==1
    err = ones(nfit,1)*P2(1); % uncertainty in offset term
else
    err = 0;
    for i = 1:nloops
        yf = real_yf + n_err.*randn(size(yf));
        if best ==2
            [A, residual] = brk_pt_fit (xf, yf, W_i);
            % E for linear case is already calculated
            % recompute A(1) in case xfit(1) is not equal to xf(1) (ccabanes,28/09/2020)
            A(1)=A(1)+A(2)*(xfit(1)-xf(1));
        elseif setbreaks
            [A, residual] = brk_pt_fit (xf, yf, W_i, breaks);
            % E stays fixed if breaks are specified
            % recompute A(1) in case xfit(1) is not equal to xf(1) (ccabanes,28/09/2020)
             A(1)=A(1)+A(2)*(xfit(1)-xf(1));
        else
            % give initial guess as the fitted break points to speed up the
            % calculation
            nbr = length(real_breaks);
            b_g = [-1 real_breaks];
            for n=1:nbr
                ubrk_g(n) = log((b_g(n+1)-b_g(n))/(1-b_g(nbr+1)));
            end
            try
                [ubrk resnorm residual] = lsqnonlin (@nlbpfun, ubrk_g, [],[], ...
                    optimset('DISPLAY','off','MaxFunEvals',MaxFunEvals, 'TolFun',TolFun) );
            catch
                [ubrk resnorm residual] = LMA(ubrk_g);
            end
            % get E matrix with new break points
            btem = [xfit(1) breaks];
            E = zeros(length(xfit),best);
            E(:,1) = ones(length(xfit),1);
            ixb = sorter (btem, xfit);
            % recompute A(1) in case xfit(1) is not equal to xf(1) (ccabanes,28/09/2020)
            A(1)=A(1)+A(2)*(xfit(1)-xf(1));
            for j = 1:best-1
                ib = find (ixb == j);% pointer to x values greater than break point j
                E(ib,j+1) = xfit(ib) - btem(j);
                ii = find (ixb > j); % pointer for break points less than the one just
                % below the x value
                if (~isempty (ii) )
                    E(ii,j+1) = btem(j+1) - btem(j);
                end  % if lower triangular.

            end  % for coefficients.
        end % end if best == 2
        err = err + (yfit - E*A).^2;
    end % end for nloops
    err = err / nloops;
    % rescale error to reflect the decrease do to the off diagonal covariances
    err = err.*P3;
end % end if best == 1

A = real_A; E = real_E; breaks = real_breaks; xf = real_xf; yf = real_yf;


%% get residual statistics for each profile
W=diag(W_i);
sta_mean = NaN*ones(1,nfit);
sta_rms = NaN*ones(1,nfit);
ip1 = [0 ; ip];
for n=1:nfit
    if ~isempty(find(xp==xfit(n)))
        index = find(xp==xfit(n));
        sta_mean(n) = sum( W(ip1(index)+1:ip1(index+1)) .* ...
            (yf(ip1(index)+1:ip1(index+1))) ) / ...
            sum(W(ip1(index)+1:ip1(index+1)));
        sta_rms(n)  = sum( W(ip1(index)+1:ip1(index+1)) .* ...
            (sta_mean(n)-yf(ip1(index)+1:ip1(index+1))).^2 ) / ...
            sum(W(ip1(index)+1:ip1(index+1)));
        sta_rms(n) = sqrt(sta_rms(n));
    end
end


%% convert back to original units
xp = xp*xscale+x0;
xfit = xfit*xscale+x0;
% pass the coefficients and breaks back
break_pts = breaks*xscale+x0;
A(1) = A(1)*yscale+y0;
P(1) = P(1)*yscale;
if length(A) > 1
    A(2:best)= A(2:best)*yscale/xscale;
    P(2:best)= P(2:best)*yscale/xscale;
end

yfit = yfit*yscale+y0;
err = sqrt(err)*yscale;
yerr = err;
n_err = n_err*yscale;
sta_mean = sta_mean*yscale+y0;
sta_rms = sta_rms*yscale;
% get time-derivatives and time-derivative errors
ixb = sorter([min(xfit) break_pts],xfit);
if best == 1
    time_deriv = zeros(nfit,1);
    time_deriv_err = ones(nfit,1)*NaN;
elseif best == 2
    time_deriv = A(2)*ones(nfit,1);
    time_deriv_err = P2(2)*ones(nfit,1);
else
    for j=1:best -1
        ib = find (ixb == j);
        time_deriv(ib)=A(j+1);
        time_deriv_err(ib)=P2(j+1);
    end
end

condslope = yfit';
condslope_err = yerr';
% return fit parameters
fit_coef = A;
fit_breaks = break_pts;

return

%% convenient function to fill out arrays for piecewise-fit matrix
function pointer = sorter(msites, sites)
% used to find the interval in piecewise linear fit
[ignored,index] = sort([msites(:).' sites(:).']);
pointer = find(index>length(msites))-(1:length(sites));


