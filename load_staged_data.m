function [A varargout]=load_staged_data(graphname)
% LOAD_STAGED_DATA Loads a .mat file from a staged directory
%
% load_staged_data is a helper function to load a graph provided with the
% regardless of the current working directory.  
%
% This requires setting the global directory: 
%
%   
%
% Example:
%   A = load_staged_data('cond-mat-2005-fix-cc');
%


% David F. Gleich
% Copyright, Purdue University, 2014

% History
% 2011-10-02: Initial coding based on load_gaimc_graph
% 2013-01-28: Added option to load coordinates

path = '/home/dgleich/research/2014/nexpokit-staged-data';
data=load(fullfile(path,[graphname '.mat']));
A = data.P;


