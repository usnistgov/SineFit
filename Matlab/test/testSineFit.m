classdef testSineFit < matlab.unittest.TestCase
    % Class-based unit testing for the Freq_Cal sine fitter
    %
    %   run these tests with two command-line commands:
    %   >> testCase = testSineFit();
    %   >> res = run(testCase);  
    %
    %   FC must be an existing Freq_Cal class instance passed into the constructor

    properties
        SF              % a SineFit class instance will be instantiated
        Hq = @(x) cos(x)  % custom periodic function handle
        fs = 1.0e6      % sample rate
        t               % time vector
        N = pow2(9)     % size of the test signal
        sigParams = struct()
        noiseParams = struct()
        harmParams = struct()
        y              % test signal
        
        % results
        params
        y_est
        y_resid
        err_rms
        
        % utility properties
        fig = 0;
        
        
    end
    
    %% Constructor
    methods (Access = public)
       
        function obj = testSineFit()
            obj.SF = SineFit();     % instantiate a SineFit object
        end
        
    end
    
%%-------------------------------------------------------------------------
    %% Test Methods
    % These functions will be called on   >> "res = run(testCase);"
    methods (Test)
        function regressionTests (obj)
            defaultSingleChanSignal(obj)
            genSignal(obj)
            %test3PFit(obj);
            test2x3PFit(obj)
        end
        
    end
 
%%------------------------------------------------------------------------- 
    % Public methods to test the fitter
    methods (Access = public)
        function test3PFit(obj)
           fprintf('3-parameter fit test') 
           obj.SF.Y = real(obj.y);
           obj.SF.t = obj.t;
           W = obj.sigParams.f*2*pi; 
           obj.SF.W = W;
           obj.SF.fitter()
           
           % plot the results
           obj.fig = obj.fig+1; figure(obj.fig)
           subplot(2,1,1)
           title('1 Channel 3-P fit')
           hold on
           obj.SF.plot(1,'Y','Y_est')
           hold off
           subplot(2,1,2)
           obj.SF.plot(1,'Y_resid')
           
        end
        
        function test2x3PFit(obj)
           fprintf('3-parameter, 2-channel, fit test\n') 
           %obj.SF.Y(:,1) = real(obj.y);
           %obj.SF.Y(:,2) = imag(obj.y)];
           obj.SF.Y =horzcat(real(obj.y),imag(obj.y));
           obj.SF.t = obj.t;
           W = obj.sigParams.f*2*pi; 
           obj.SF.W = W;
           obj.SF.fitter(1)
           fprintf('iter: %d\n',obj.SF.iter)
           
           % plot the results
           obj.fig = obj.fig+1; figure(obj.fig)
           subplot(2,1,1)
           hold on
           title('2 Channel 3-P fit')
           obj.SF.plot(1,'Y','Y_est')
           hold off
           subplot(2,1,2)
           hold on
           obj.SF.plot(2,'Y','Y_est')
           hold off
            
            
        end
        
    end
    
%%--------------------------------------------------------------------------
    % Private methods for signal generation
    methods (Access = private)
       
        function defaultSingleChanSignal(obj)
            obj.fs = 1.0e6;
            ts = 1/obj.fs;
            obj.t = (0:(obj.N-1)).'*ts;
            
            % single tone quadrature signal
            obj.sigParams.f = 6667;
            obj.sigParams.Ai = 1.2;
            obj.sigParams.Aq = 0.8;
            obj.sigParams.phi = pi/3;
            obj.sigParams.phq = pi/4;
            obj.sigParams.offi = 0.8;
            obj.sigParams.offq = 1.2;
            
            % add some noise and harmonics
            obj.noiseParams.stdJit = ts*1.0e-1;     % standard deviation of jitter noise
            obj.noiseParams.stdNoi = 1.0e-3;        % standard deviation of noise added to signal
            obj.noiseParams.stdPha  = pi/36;       % standard deviation of phase noise
            obj.harmParams.A2I = 5.0e-2;           % 2. harmonic proportional to Ai
            obj.harmParams.A3I = 1.0e-2;           % 3. harmonic proportional to Ai
            obj.harmParams.A2Q = 4.0e-2;           % 2. harmonic proportional to Aq
            obj.harmParams.A3Q = 2.0e-2;         	% 3. harmonic proportional to Aq
            
        end
        
        function genSignal(obj)
            %-------------------
            % generate Noise
            rng(693653) % seed number for repeatability
            n_i     = randn(obj.N,1);	n_i     = obj.noiseParams.stdNoi*n_i/std(n_i);
            n_q     = randn(obj.N,1);   n_q     = obj.noiseParams.stdNoi*n_q/std(n_q);
            n_phi   = randn(obj.N,1);   n_phi   = obj.noiseParams.stdPha*n_phi/std(n_phi);
            n_phq   = randn(obj.N,1);   n_phq   = obj.noiseParams.stdPha*n_phq/std(n_phq);
            n_jiti  = randn(obj.N,1);   n_jiti  = obj.noiseParams.stdJit*n_jiti/std(n_jiti);
            n_jitq  = randn(obj.N,1);   n_jitq  = obj.noiseParams.stdJit*n_jitq/std(n_jitq);
            %--------------------
            
            % Generate noisy single-tones and a noisy complex sinusoid
            % clock jitter
            ti = obj.t + n_jiti;
            tq = obj.t + n_jitq;
            % sinusoid with phase noise and clock jitter
            w = 2*pi*obj.sigParams.f;
            yi = obj.sigParams.Ai*cos(w*ti+obj.sigParams.phi+n_phi);
            yq = obj.sigParams.Aq*sin(w*tq+obj.sigParams.phq+n_phq);     % Hq is the sin(...) function
            % add offset, noise & harmonics
            yi = obj.sigParams.offi + yi + obj.harmParams.A2I*yi.*yi + obj.harmParams.A3I*yi.*yi.*yi + n_i;
            yq = obj.sigParams.offq + yq + obj.harmParams.A2Q*yq.*yq + obj.harmParams.A3Q*yq.*yq.*yq + n_q;
            obj.y = yi+1i*yq;

        end
        
    end
    
end
    