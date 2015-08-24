function h = stairs2(varargin)
% function h = stairs2(...)
%
% A stairs plot that centers each bin on the x value, instead of going to the right of it.
% For use with [n,x] = hist(...).
%
% Input arguments are the same as stairs. First argument must be the vector of x values.
%
% Returns h, the handle for the plot object.

x = varargin{1};    %the vector of x values

xoffset = 0.5 * mode(x(2:end) - x(1:end-1));    
%we want to shift the graph by half of the most frequent / obvious delta-x

%now, we need to call stairs with all the original input arguments. hrm.
functionString = 'stairs(varargin{1} - xoffset';
for i=2:nargin
    functionString = [functionString, ', varargin{',num2str(i),'}'];
end
functionString = [functionString, ');'];
%now functionString contains all the varargin arguments.

h = eval(functionString);