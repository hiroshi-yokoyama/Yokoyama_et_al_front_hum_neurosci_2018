% Kuramoto_sim_v4() - Solve the differential equation of Kuramoto-like N-th phase coupling oscillator  
% [input]
%  - Nosc       : # of coupling oscillotor
%  - F          : frequency parameter for eigen frequecny (Hz)
%  - tspan      : row or column vector of time periods
%  - time_index : [1 x 2] input vector of event time-epoch (sec)
%                 e.g.) If time_index = [0.1, 0.3], 
%                       the oscillator is synchronized in the time period between 0.1 to 0.3 sec. 
%  - rnd_seed   : the settings of the random number generator
%                 The data structure is in accordance with the "rng" functions  (See the details in descriptions of rng() function)
%                 If you ignor this parameter, the oscillator model is calculated with current settings of random number generator.
% [output]
% - Tout        : [sample x 1] vector of time periods (sec)
% - theta       : phase angle in each oscillator  [Nosc x time sample]
% - signal      : time-domain signals in each oscillator [Nosc x time sample]
% - rnd_seed    : the settings of the random number generator
%                 The data structure is in accordance with the "rng" functions  (See the details in descriptions of rng() function)
function [Tout,theta, signal, rnd_seed]=Kuramoto_sim_v4(Nosc, F, tspan, time_indx, rnd_seed)
%     Nosc = 4;

    T = tspan;
    theta_st = pi/2.*ones(Nosc, 1);
    
    if nargin<=4
        rnd_seed.omega = rng;
        omega = 2*pi*unifrnd(0, F, Nosc, length(T));
        rnd_seed.zeta = rng;
        zeta   = normrnd(0,30, [Nosc,length(T)]);
    else
        rng(rnd_seed.omega);
        omega = 2*pi*unifrnd(0, F, Nosc, length(T));
        rng(rnd_seed.zeta);
        zeta   = normrnd(0,30, [Nosc,length(T)]);
    end
    
    [Tout,theta] = ode113(@Nkuramoto,T,theta_st);
    theta = mod(theta, 2*pi)-pi;
    
    signal = real(ifft(exp(1i.*theta')));

    function thetaDot = Nkuramoto(t, theta)
        thetaDot = zeros(size(theta));
        PRC = zeros(size(theta)); 
        K = 0.03;
        
        t_indx = find(t<=tspan,1);
        % Instrinsic part
        for i=1:Nosc
            if t >= time_indx(1) && t <= time_indx(2)
                Noise  = 0;
                OMEGA = 2*pi*F;
                for j=1:Nosc
                    if i~=j
                        PRC(i) = K*sin(theta(j)-theta(i));
                    end
                end
            else
                Noise = zeta(i, t_indx);
                OMEGA = omega(i, t_indx);
                for j=1:Nosc
                    if i~=j
                        
                        PRC(i) = K*sin(theta(j)-theta(i));
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetaDot(i) = OMEGA + PRC(i) + Noise;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end

end
