classdef SineFit < handle
    % 3 or 4 parameter sine fitting class
    
    properties
        %data
        Y           % matrix of sample data, rows of data in columns of channels 
                    % data may be real or complex
        t           % time, may be scalar sample period(s) or vectors of channel time series
                    % there can be a single column of time or columns of channel times
        W           % scalar or matrix of angular frequencied in rad/s
                    %   scalar: known frequency or initial guess for all channels
                    %   matrix: columns of known or initial guesses per4 channel
                    %   initial guesses may contain 2 or more rows of w - will force 4-parameter fit
                    %      using the mean of the column of w values   
        H = @cos    % sine function handle.  may be @sin or (default) @cos
        
        % fitter options
        TolX = 1e-4 % 4=parameter fit stopping criteria
        MaxIter = 500  % max number of 4-P fit iteations

        % fit results            
        params      % cell array fit parameters rows of {Offs, A, phi, w} in columns of channels
        Y_est       % fitted sine curves, same format as Y
        Y_resid     % Y - Y_est
        rmserr      % Root-Mean-Squared of the residual
        iter        % number of iterations of 4-parameter fit
        exitflag    % 1: converged, 0: MaxIter exceeded
               
    end
    
    %% Constructor
    methods
        function obj = SineFit(varargin)
            % Constructor for SineFit Class
            %   USEAGE:    obj = SineFit(<name-value pairs>)
            %       name-value pairs:
            %           'Y':  matrix of sample data, rows of data in columns of channels 
            %           't': time, may be scalar sample period(s) or vectors of channel time series
            %                   there can be a single column of time or columns of channel times
            %           'W (or w)': scalar or matrix of angular frequencied in rad/s 
            %           'TolX': 4-parameter fitter stopping criteria
            %           'MaxIter': 4-parameter fit maximum number of iterations before stopping
    
            % handle the input arguments
            i = 1;
            while i <= nargin
               switch varargin{i}
                   case {'Y','y'}
                       i = i+1;
                       obj.Y = varargin{i};
                   case {'T','t'}
                       i = i+1;
                       obj.t = varargin{i};
                   case {'W','w'}
                       i = i+1;
                       obj.W = varargin{i};
                   case 'TolX'
                       i = i+1;
                       obj.TolX = varargin(i);
                   otherwise                      
                       warning('Unrecognized parameter %s',varargin{i});
               end
            end
            i = i+1;
            end
    end
        
%%-------------------------------------------------------------------------
% Public methds found in external method files
    methods (Access = public)
        %[obj.params,obj.Y_est,obj.Y_resid,obj.rmserr,obj.iter,obj.exitflag]=fitter(obj,P4)
        fitter(obj,varargin)
        plot(obj,varargin)
    end
    
end    