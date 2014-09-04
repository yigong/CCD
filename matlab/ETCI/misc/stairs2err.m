function h = stairs2err(varargin)
% function h = stairs2err(Y,E)
% function h = stairs2err(X,Y,E)
% function h = stairs2err(X,Y,L,U)
% function h = stairs2err(...,'linespec')
% function h = stairs2err(...,'linespec','key',value)
%
% A stairs plot that centers each bin on the x value, instead of going to the right of it.
% Also adds one more horizontal line at the right edge of the stairs, for the last bin.
% For use with [n,x] = hist(...).
%
% Error bars are vertical lines through the horizontal center of each "step" of stairs.
%
% Key names and values should match those of plot or stairs. e.g. 'linewidth', 2
%
% Returns h, the handle for the plot object.

v = varargin;
[x,y,lowerbar,upperbar,linespec,key,val] = inputHandling(nargin,v);

holdState = get(gca,'NextPlot');    %either 'next' or 'add'
hold on;

xoffset = 0.5 * mode(x(2:end) - x(1:end-1));    
%we want to shift the graph by half of the most frequent / obvious delta-x

h = zeros(1,length(x)+1);
if isempty(key)
    %make stairs plot
    h(1) = stairs(x - xoffset,y,linespec);
    %finish stairs plot - last stair
    h(2) = plot([x(end)-xoffset; x(end)+xoffset], [y(end); y(end)], linespec);
    
    %make error bar lines
    for i=1:length(x)
        htmp = plot([x(i); x(i)], [y(i)-lowerbar(i); y(i)+upperbar(i)], linespec);
        h(2+i) = htmp;
    end
else
    %make stairs plot
    h(1) = stairs(x - xoffset, y, linespec, key, val);
    %finish stairs plot - last stair
    h(2) = plot([x(end)-xoffset; x(end+xoffset)], [y(end); y(end)], linespec, key, val);
    
    %make error bar lines
    for i=1:length(x)
        h(2+i) = plot([x(i); x(i)], [y(i)-lowerbar; y(i)+upperbar], linespec, key, val);
    end
end

%return hold to previous state
set(gca,'NextPlot',holdState);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,lowerbar,upperbar,linespec,key,val] = inputHandling(nargin,v)

%default
linespec = 'k';
key = '';
val = [];

%check which arguments are strings
stringFlags = false(1,nargin);
for i=1:nargin
    stringFlags(i) = ischar(v{i});
end

%how many data arguments are there? data = before any string arguments given
%--this section must define x, y, lowerbar, upperbar
if ~isempty(stringFlags)
    nDataArgs = find(stringFlags,1) - 1;
else
    nDataArgs = nargin;
end
switch nDataArgs
    case 1
        error('Not enoguh input arguments...')
    case 2
        % Y, E
        y = v{1};
        lowerbar = v{2};
        upperbar = lowerbar;
        x = 1:length(y);
        if ~all(size(y==lowerbar))
            error('Data arguments not the same size...')
        end
        
    case 3
        % X, Y, E
        x = v{1};
        y = v{2};
        lowerbar = v{3};
        upperbar = lowerbar;
        if ~all(size(y==lowerbar)) || ~all(size(x==y))
            error('Data arguments not the same size...')
        end
    case 4
        % X, Y, L, U
        x = v{1};
        y = v{2};
        lowerbar = v{3};
        upperbar = v{4};
        if ~all(size(y==lowerbar)) || ~all(size(x==y)) || ~all(size(y==upperbar))
            error('Data arguments not the same size...')
        end
    otherwise
        error('Too many data arguments...')
end

%what about linespec or key/values?
if ~isempty(stringFlags)
    linespec = v{find(stringFlags,1)};
    if nargin > find(stringFlags,1)
        %there is an additional key/value
        key = v{find(stringFlags,1)+1};
        val = v{find(stringFlags,1)+2};
    end
    if nargin > find(stringFlags,1) + 2
        error('Sorry, stairs2err only handles one key/value at this time. Go code another.')
    end
end