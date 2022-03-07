function [tf] = endsWith(str, suffix)
% [tf] = endswith(str,suffix)
%
% Logical returns true if the string ends in the specified suffix
%
% Input
% -----
% str: string in question.
% suffix: suffix of interest.
% 
%
% Output
% ------
% tf: true or false.
%
%
% Notes
% -----
% This file is from matlabtools.googlecode.com
% Modified by BPS 4 Sept 2012
str = {str};
n = length(suffix);
tf = zeros(length(str),1);

for i=1:length(str)
    foo = char(str(i));
    if length(foo) < n
        tf(i) =  false;
    else
        tf(i) = strcmp(foo(end-n+1:end), suffix);
    end
end

end