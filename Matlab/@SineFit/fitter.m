function varargout = fitter(obj,varargin)
% IEEE std 1057 Standard for Digitizing Recorders
% Fit to sinewave data.  Performs either a three-parameter fit to data with 
% a known frequency or a four-parameter fir if the frequency input is an initial
% first guess
%
% Usage:
% params = sinefit(y,t,w)
% params = sinefit(y,t,w,TolX)
% params = sinefit(y,t,w,TolX,MaxIter)
% [params,y_est,y_err,rmserr,iter,exitFlag] = sinefit(...)
% [...] = sinefit(H,y,t,w,...)
%
%INPUTS:
% y             Input data 
% t (scalar)    Sampling interval for linearly spaced data
% t (vector)    Time Vector, can be unequally spaced
% w (scalar)    Angular frequency in rad/sec, either known or initial guess
% w (vector)    Will use 4-parameter fit to find a frequency
%               between the two values. w must be length < 2
% H             Optional fit function handle e.g. @sin. Default H = @cos.
%
% TolX          For 4-parameter fit, the termination tolarance
%                Default Tolx = 1e-4
%                Forces a 4-parameter fit even if w is a scaler
% MaxIter       For 4-parameter fit, the manimum number of iterations 
%               Default MaxIter = 500
%
%OUTPUTS:
% params        Cell Array of fit parameters {Offs,A,phi,w}
%               - Offs  the offset
%               - A     the peak amplitude
%               - phi   the phase in rad
%               -w      the angular frequency in rad/s
% y_est          best fit time series
% y_err          residual time series
% rmserr        Root-Meam_Squared of the residual
% iter          number of iterations of 4-parameter fit
% exitflag      1: converged, 0: MaxIter exceeded
%
% Notes:
%   This can handle a matrix of data (and times).  The rows should be the
%   samples and the columns should be the channels.  If time is a vector,
%   one column will work for all channels or there can be a time row for
%   each  channel.  Same for other parameters, for example, if there are
%   two angular frequencies, then there can be a single column of two rows
%   or a column of pairs for each channel.
%
%   The assumption is that there are more samples per channels then there
%   are channels. If this is not the case, then an error will be thrown and
%   NAN will be returned
%

P4 = false;     % Should a 4-parameter fit be performed?
if nargin>1
    if varargin{1} == 1 
        P4 = true;
    end
end%-----------------




%------------------
% get the data
y = obj.Y;
t = obj.t;
W = obj.W;

% optional arguments
H = obj.H;  % H, a function handle


% First step is to check that the channel data is columns of samples.  If
% there are more columns than rows, throw an error.
[rows,cols]=size(y);   
if cols>rows
    disp('y argument has are more columns than rows.  It is assumed that there will be more samples than columns of data channels.')
    disp('data must be rows of samples and columns of channels')
    varargout(1:nargout) = {nan};
    return
end

% ensure the data channels are in columns 
% nS = number of Samples, nC = number of Channels
if rows == 1, y = y(:); nS=cols; nC = rows; else, nS = rows; nC = cols; end 

% ensure the time scalar or time vector channels are in columns
% nT = number of time vectors or scalars, ncT is the number of channels of time info
[rowt,colt]=size(t);
    nT = rowt;
    ncT = colt;
if (rowt == 1 && colt ~= nC)
    t = t(:);
    nT = colt;
    ncT = 1;
elseif (rowt ~= nS && colt == nC)
    disp('Error: the time vector must be the same length as the data vectors')
    varargout(1:nargout) = {nan};
    return
end
%     
    
% ensure the frequency info channels are in columns   
% nW = number of frequency rows, ncW = number of frequency channels
[roww,colw]=size(W);
    nW = roww;
    ncW = colw;
if (roww == 1 && colw ~= nC)
    W = W(:);
    nW = colw;
    ncW = roww;
end
if nW > 1, P4 = true; end

% pre-allocate the results
params = cell(4,nC);
exitFlag = ones(1,nC);
iter = zeros(1,nC);
% loop through each of the channels and perform the fit
onevec = ones(nS,1);
colY = 1;   % index for th3e sample channel
colT = 1;   % index for the time (scalar or vector)
colW = 1;   % index for the frequency (scalar or vector)

while colY <= nC
    
    % time vector, either they are equal length or t is a scaler
    tvec = t(:,colT);
    if not(nS==nT)
        if nT==1
            tvec = (0:(nS-1)).'*t(:,colT);
            fprintf('Channel %d fitted using linearly equally spaced time vector\n',colY)
        else
            disp('Error: the data and time vectors lengths are unequal.')
            varargout(1:nargout) = {nan};
            return
        end
    end
      
    if length(W(:,colW))>1 % indicates the need to find best w between two values
        w = mean(W(:,colW));        
    else
        w = W(1,colW);
    end
    
    wt = w*tvec;    
    % 3-parameter sinefit - performed for both 3 and 4-parameter fitting
    cosvec=H(wt);
    sinvec=H(wt-pi/2);
    D0=[cosvec sinvec onevec];
    [Q,R] = qr(D0,0);
    x0 = R\(Q.'*y(:,colY));
    
    %4 - parameter sinefit
    if P4
        fprintf('Performing 4-Parameter Fit on channel %d\n',colY)
        x0=[x0;0];
        success = 0;
        while success==0
            iter(1,colY) = iter(1,colY)+1;
            x0_old=x0;
            w=w+real(x0_old(4));
            wt = w*tvec;
            cosvec=H(wt);
            sinvec=H(wt-pi/2);
            D0=[cosvec sinvec onevec -x0(1)*tvec.*sinvec+x0(2)*tvec.*cosvec];
            if isreal(y(:,colY))
                [Q,R] = qr(D0,0);
                x0 = R\(Q.'*y(:,colY));
            else
                x0=lscov(D0,y(:,colY));
            end
            err = max(abs((x0(1:3)-x0_old(1:3))));            
            if err < obj.TolX || iter(1,colY) > obj.MaxIter
                success = 1;
            end            
        end
        x0(4)=real(w+x0(4));
        if iter(1,colY) > obj.MaxIter
            exitFlag(1,colY) = 0; % no success or incomplete fit
        else
            exitFlag(1,colY) = 1;
        end

    else  % not P4
        fprintf('Performing 3-Parameter Fit on channel %d\n',colY)
        
    end % if P4
    [params(1:3,colY)]=getParams(x0);
    params{4,colY}=w;
    colY = colY+1;
    if colT > 1, colT = colT+1; end
    if colW > 1, colW = colW+1; end      
end
    
% set up the output
    obj.params=params;


% y_est output
y_est = zeros(nS,nC);
colY = 1;
colT = 1;
while colY <= nC
    [Offs,A,pha,w] = params{:,colY};
    
    if length(A)==1
        y_est(:,colY) = Offs + A*H(w*t(:,colT)+pha);
    else
        y_est(:,colY) = Offs(1) + A(1)*H(w*t(:,colT)+pha(1))+...
            1i*(Offs(2) + A(2)*H(w*t(:,colT)+pha(2)-pi/2));
    end
    
    colY = colY+1;
    if colT > 1, colT = colT+1; end
end
obj.Y_est = y_est;


% y_err output (residual)
y_err = zeros(nS,nC);
colY = 1;
while colY <= nC
    y_err(:,colY) = y(:,colY)-y_est(:,colY);
    colY = colY+1;
end
obj.Y_resid = y_err;

rmserr = zeros(1,nC);
colY = 1;
while colY <= nC
    rmserr(1,colY) = norm(y_err(:,colY))/sqrt(length(y_err(:,colY)));
    colY = colY+1;
end
obj.rmserr = rmserr;

obj.iter=iter;
obj.exitflag=exitFlag;
end

%%-------------------------------------------------------------------------
% local functions
function [params] =  getParams(x0)
A0 = x0(1);
B0 = x0(2);
if isreal(x0)
    A = (A0^2+B0^2)^.5;
    phi = atan(-B0/A0);
    Offs = x0(3);
    if A0<0,phi=phi+pi;end
else
    B0=-B0*1i;
    A0_old = A0;
    A0 = real(A0)+1i*imag(B0);
    B0 = real(B0)+1i*imag(A0_old);
    
    A(1)=abs(A0);
    A(2)=abs(B0);
    
    phi(1) = angle(A0);
    phi(2) = angle(B0);
    Offs(1)=real(x0(3));
    Offs(2)=imag(x0(3));
end
params = {Offs,A,phi};

end


