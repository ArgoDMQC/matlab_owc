function cmap = make_colorbar(colorbar_limits,color_level_boundaries,colortable)

% eg colorbar_limits = [-.1 .1];
% eg color_level_boundaries = [ -.06 -.04 -.02 -.01 -.005 .005 .01 .02 .04 .06];
% eg colortable = [
%         0         0    0.7000
%         0         0    0.9375
%         0    0.3125    1.0000
%         0    0.6875    1.0000
%    0.8000    1.0000    1.0000
%    1.0000    1.0000    1.0000
%    1.0000    1.0000    0.6000
%    1.0000    0.8750         0
%    1.0000    0.6000         0
%    1.0000    0.4000         0
%    1.0000    0.1250         0
%    ]

%In order to get nice blocks of color in the colorbar that
%correspond to the intervals between contours, define the
%color level boundaries and the rgb values of colours that
%fall between them. The color level boundaries don't have to be
%the same as the contours. ie you can have several contours between
%a color level boundary.

%this program builds a colormap by dividing the colorbar into evenly
%spaced blocks, and then repeating each defined colour the correct
%number of times so that the final colorbar changes color at the
%defined color levels.

%The size of the block is chosen with a Highest Common Factor routine.
%The routine uses a scaling, so that decimal color level boundaries
% appear as 'integers'. The default scaling is 1e6, so that
% boundaries with any reasonable number of decimal places produce
% the expected behaviour.


if nargin < 3
    colortable = jet(length(color_level_boundaries)+1);
end

%we need to know how many total blocks of color we require between colorbar_limits.
%then we'll work out how often to repeat each one.

%First check there is one more colour than colour level.
%we expect colors to be defined between each colour level, and above and below the
%first and last colour levels.

if length(colortable) ~= length(color_level_boundaries)+1
    disp('error')
    return
end

%size of colorblock is minimum gap between contours, or between
%contours and colorbar_limits.

%first find countours that fall within colorbar_limits

iok = find(color_level_boundaries > colorbar_limits(1) & color_level_boundaries < colorbar_limits(2));
use_levels =  [colorbar_limits(1) color_level_boundaries(iok) colorbar_limits(2)];
iok = [iok iok(end)+1]; %include one extra colour at the end; this is now the set of colours in the colorbar

diflev = diff(use_levels);
blocksize = hcf(diflev);  %highest common factor of diflev; use hcf default scaling.

clear cmap
numc = 0;

for k = 1:length(use_levels)-1
    numblock = round((use_levels(k+1)-use_levels(k))/blocksize);
    for kk = 1:numblock
        numc = numc+1;
        cmap(numc,:) = colortable(iok(k),:);
    end
end

return

function h = hcf(list,factor)

%find the highest common factor of a list of numbers.
%
%This was originally written to work with non-integer
%contour intervals, so I allow a user scaling to convert
%'reals', eg 0.01 0.02 and so on, into 'integers'

%algorithm from PDK; 7 Mar 2005; coded by BAK.

if nargin == 1
    factor = 1000000;
end

list = list(:)';  %force to a row vector
list = round(list*factor); %apply scaling and convert to integers
list = unique(list);  %eliminate duplicates and sort ascending

if list(1) < 1
    h = nan;
    return
    %I only deal with sensible lists of positive numbers !
end


%Now we have real cases to deal with

while length(list) > 1
    list = unique([list(1) list(2:end)-list(1)]); %subtract smallest, sort and eliminate duplicates
    if list(1) == 1  %hcf must be unity, so bale out
        list = 1;
    end
end

h = list/factor;

return


    
