
function [S_h, P_h] = interp_climatology(S, Theta, P, S_f, Theta_f, P_f)

% [S_h, P_h] = interp_climatology(S, Theta, P, S_f, Theta_f, P_f)
% routine to interpolate climatological station salinity and pressure
% values onto float theta's
%
% INPUT:
%           S       salinity for climatology station data
%           Theta   potential temperature for climatology
%           P       Pressure for climatology
%                   these values should be passed as [N_level, N_station]
%                   where N_level is the number of levels and
%                   N_station is the number of stations
%
%           S_f     salinity for float
%           Theta_f potential temperature for float
%           P_f     Pressure for float
%                   these values are assumed to be a vector of values for
%                   the float [length = N_fl]
%
% OUTPUT:
%           S_h     interpolated salinity value on float theta surface
%           P_h     interpolated pressure value on float theta surface
%                   these values will be returned in matrices of 
%                   size [N_fl, N_station]
%                   
% Breck Owens, July 2005
%


[N_level, N_stat] = size(S);
N_fl = length(S_f);

S_h = ones(N_fl,N_stat)*NaN;
P_h = ones(N_fl,N_stat)*NaN;

% make sure that the station data has finite values and discard possible
% bad values in the middle of a profile and truncate the number of levels
% to the maximum level that has data
max_level = 0;
for n=1:N_stat
    ok = find(isfinite(S(:,n)) & isfinite(Theta(:,n)) & isfinite(P(:,n)));
    if ~isempty(ok)
        lok = length(ok);
        S(1:lok,n) = S(ok,n);
        Theta(1:lok,n) = Theta(ok,n);
        P(1:lok,n) = P(ok,n);
        max_level = max(max_level,lok);
    end
end
if max_level > 0
    S = S(1:max_level,:);
    Theta = Theta(1:max_level,:);
    P = P(1:max_level,:);
    N_level = max_level;
else
    disp('Error in INTERP_CLIMATOLOGY, No good climatological data found')
    return
end

ok = find(isfinite(S_f) & isfinite(Theta_f) & isfinite(P_f));
nok = length(ok); % number of finite float observations
for lok = 1:nok
    nl = ok(lok);
    DP = P - P_f(nl);
    DS = S - S_f(nl);
    DT = Theta-Theta_f(nl);

    [val,ia] = min(abs(DP),[],1); % will give vector of level for each station that has P closest to float P_f
    clear val
    
    for n=1:N_stat % loop through all climatology stations
        tst = sign(DT(:,n))*sign(DT(ia(n),n)); % tst = -1 when DT is different from DT for closest P
        S1 = [];
        P1 = [];
        i1 = min(find(tst(ia(n):N_level) < 0)); % look for a theta match below the float pressure
        % 1st point to change sign after closest pressure point, the 
        % interpolated value for DT = 0 will be between level i1 and i1-1
        i2 = max(find(tst(1:ia(n)) < 0)); % look for a theta match above the float pressure
        % last point to change sign before closest pressure point, the
        % interpolated value for DT = 0 will be between level i2 and i2+1
        if ~isempty(i1)
            % there is a theta value at a deeper level
            i1 = i1+ia(n)-1;
            wt = DT(i1,n)/(DT(i1,n)-DT(i1-1,n));
            P1 = wt*DP(i1-1,n)+(1-wt)*DP(i1,n);
            S1 = wt*DS(i1-1,n)+(1-wt)*DS(i1,n);
        else
            P1 = [];
            S1 = [];
        end
        if ~isempty(i2)
            % there is a theta value at a lower level
            wt = DT(i2,n)/(DT(i2,n)-DT(i2+1,n));
            P1 = [P1 wt*DP(i2+1,n)+(1-wt)*DP(i2,n)];
            S1 = [S1 wt*DS(i2+1,n)+(1-wt)*DS(i2,n)];
        end
        if ~isempty(P1)
            % if there are two nearby values of theta, choose the closest one
            [P_h(nl,n),k] = min(abs(P1));
            P_h(nl,n) = P1(k)+P_f(nl);
            S_h(nl,n) = S1(k)+S_f(nl);
        end
    end
end
