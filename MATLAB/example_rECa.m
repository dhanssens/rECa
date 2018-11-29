%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%                              EXAMPLE FILE                               %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Use:
%  Calculates the rECa (S/m) of a 1D/2D QP (and IP) dataset (ppm). The 
%  calculation is based on the assumption of a 1D stratified half-space. 
%  Typical sensor characteristics are stored in the sensor structure (S). 
%  The input structure (I) contains both IP (I.IPdata)) and QP data 
%  (I.QPdata; ppm). An initial ECa estimation (S/m) could be included in 
%  order to guide the rECa estimation and ensure an increased process 
%  speed. If no initial ECa estimation is included, an Amplitude based 
%  approach is used. The forward model used, can be found at 
%  github.com/dhanssens.
%
%  Created by Daan Hanssens
%  Ghent University, Belgium
%
%  Cite:
%  Hanssens, D., Delefortrie, S., Bobe, C., Hermans, T., De Smedt, P., 2019. 
%  Improving the reliability of soil EC-mapping: Robust apparent electrical 
%  conductivity (rECa) estimation in ground-based frequency domain electromagnetics. 
%  Geoderma 337, 1155-1163.
%

   clc; clear all; close all;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- User-input -------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % Sensor characteristics (S structure)
    %

        S.x=        4;                                                     % x-coordinate receiver (m)
        S.y=        0;                                                     % y-coordinate receiver (m)
        S.z=        -0.16;                                                 % z-coordinate receiver (m) - positive z-axis pointed down
        S.height=   0.16;                                                  % Height of transmitter (m)
        S.freq=     9000;                                                  % Frequency (Hz)
        S.mom=      1;                                                     % Transmitter moment (A.m^2)
        S.ori=      'ZZ';                                                  % Coil orientation (Transmitter(X,Y,Z),Receiver(X,Y,Z))


    %
    % QP (and IP) input (I structure)
    %

        I.QPdata= QP;                                                      % QP input (ppm)
        I.IPdata= IP;                                                      % IP input (ppm), optional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- Output ------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % Calculate rECa (Amplitude based initial estimation)
    %
    
        rECa_ampl_based = calc_rECa(S,I);
     
        
    %
    % Calculate rECa (User-supplied initial estimation)
    %
    
        init = .4;                                                         % Init EC estimation (S/m)
        rECa_init_based = calc_rECa(S,I,init);







