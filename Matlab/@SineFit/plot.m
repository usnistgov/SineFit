function plot(obj,channel,varargin)
% plot override functions
%   Plots on the existing plot or subplot
%   Useage: plot(channel,<options>
%

% handle the options
i = 1;
while i <= nargin-2
    switch varargin{i}
        case{'y','Y'}
            % plot the original signal;            
            t = getChannelTime(obj,channel);
            plot(t,obj.Y(:,channel),'-k')
        case{'y_est','Y_est','Y_Est'}
            % plot the estimated signal
            t = getChannelTime(obj,channel);
            plot(t,obj.Y_est(:,channel),'-r')
        case{'y_resid','Y_resid','Y_Resid'}
            % plot the estimated signal
            t = getChannelTime(obj,channel);
            plot(t,obj.Y_resid(:,channel),'.b')
            
        otherwise
            warning('unknown argument %s',varargin{i})
    end
    i = i+1;
end

end

%%=========================================================================
% Local Functions
function [t] = getChannelTime(obj,channel)
    [row,col] = size(obj.t);
    
    % find out if therre is a per-channel time
    if col == 1, channel = 1; end
    
    nS = size(obj.Y,1);     % length of the data
    if (row == 1)
        % obj.t is a scalar for all channels
        t = (0:(nS-1)).'*obj.t(:,channel);
    else
        t = obj.t(:,channel);
    end
        
end