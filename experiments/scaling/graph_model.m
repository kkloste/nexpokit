function A = graph_model(type,n,varargin)
% GRAPH_MODEL Create the output of a graph model.
%
% A = graph_model(type,n,'parameter',value,'parameter',value,...)
% creates a graph using the model type with n vertices and where
% the parameters are specified in a list.
%
% Current types are:
%   'forest-fire','preferential-attachment','copying'
% 
% Parameters for each type are:
%   'forest-fire':
%     'prob': probability of following an edge (required)
%     'initial': size of initial clique (optional, default=2)
%
%   'copying':
%     'prob': probability of making an error copying an edge (required)
%     'initial': size of initial clique (optional, default=2)
%
%   'preferential-attachment':
%     'initial': size of the initial clique (optional, default=2)
%
% This function calls the C++ code in the same directory and reads the
% output from the file matlab.smat
%
% Example:
%   A = graph_model('forest-fire',n,'prob',0.5);
%   A = graph_model('forest-fire',n,'prob',0.5,'initial',5);
%

% History
% :2011-09-28: Initial coding

opts = struct(varargin{:});

mydir = fileparts(mfilename('fullpath'));
%output = fullfile(mydir,'matlab.smat');
output = tempname;
args = [];

switch type
    case {'forestfire','forest-fire','forest_fire'}
        prog = 'forest_fire_model';
        prob = opts.prob;
        initial = 2;
        if isfield(opts,'initial')
            initial = opts.initial;
        end
        args = [args sprintf(' -n %i',n)];
        args = [args sprintf(' -m %i', initial)];
        args = [args sprintf(' -p %i', prob)];
        args = [args sprintf(' -o %s',output)];
        
    case {'copying','copyingmodel'}
        prog = 'copying_model';
        
        prob = opts.prob;
        initial = 2;
        if isfield(opts,'initial')
            initial = opts.initial;
        end
        args = [args sprintf(' -n %i',n)];
        args = [args sprintf(' -m %i', initial)];
        args = [args sprintf(' -p %i', prob)];
        args = [args sprintf(' -o %s',output)];
        
    case {'pa','preferential_attachement'}
        prog = 'preferential_attachment_model';
        
        initial = 2;
        if isfield(opts,'initial')
            initial = opts.initial;
        end
        args = [args sprintf(' -n %i',n)];
        args = [args sprintf(' -k %i', initial)];
        args = [args sprintf(' -o %s',output)];
        
    otherwise
        error('graph_model:invalidArg','unknown type "%s"',type);
end

if exist(output,'file')
    delete(output);
end

cmd = sprintf('%s %s',fullfile(mydir,prog),args);
system(cmd);

A = fast_smat_read_mex(output);

delete(output);
