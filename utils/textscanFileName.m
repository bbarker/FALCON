function [C, position] = textscanFileName(fileName, formatSpec, varargin)
%
% A simple wrapper for textscan that takes a fileName instead of
% a fileID as its primary argument.
%
% Brandon Barker    9/29/2013
evalStr = '[C, position] = textscan(fileID, formatSpec';
for i = 1:length(varargin)
    evalStr = [evalStr ', varargin{' num2str(i) '}'];
end
evalStr = [evalStr ');'];
fileID = fopen(fileName, 'r');
eval(evalStr);
fclose(fileID);