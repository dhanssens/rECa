%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%        CALCULATE ROBUST APPARENT ELECTRICAL CONDUCTIVITY (rECa)         %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  rECa = calc_rECa(S,I,init)
%
%  Description:
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
%  Input:
%  I.QPdata             QP data (ppm)
%  I.IPdata             IP data (ppm)
%  S                    Sensor structure (see github.com/dhanssens)
%  init                 Initial ECa (S/m) estimation (optional)
%
%  Output:
%  rECa                 rECa (S/m) estimation
%
%  Cite:
%  Hanssens, D., Delefortrie, S., Bobe, C., Hermans, T., and P., De Smedt, 
%  2018. Improving the reliability of soil EC-mapping: robust apparent 
%  electrical conductivity (rECa) estimation in ground-based frequency 
%  domain electromagnetics. Geoderma.
%
%  Additional reference:
%  Huang, H., Won, I.J., 2000. Conductivity and Susceptibility Mapping 
%  Using Broadband Electromagnetic Sensors. Journal of Environmental & 
%  Engineering Geophysics 5(4), 31-41.
%
%  Author:
%  Created by Daan Hanssens.
%  Ghent University, Belgium, 2017.
%
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function rECa = calc_rECa(S,I,init)
    
    %
    % Model characteristics (M structure)
    %

        M.con= [];                                                         % Conductivity of layer(s) (S/m)
        M.sus= 0;                                                          % Susceptibility of layer(s) (SI unit)
        M.perm= 1/(35950207149.4727056*pi);                                % Permittivity of layer(s) (F/m)
        M.thick= 0;                                                        % Layer(s) thickness (m)


    %
    % Get extra information about input structure
    %

        dim= size(I.QPdata);                                               % Dimensions
        blank= I.QPdata .* 0 + 1;                                          % Boolean blank
        rQPdata= reshape(I.QPdata,[],1);                                   % Reshaped QP data

        
    %
    % Generate EC (S/m) - QP (ppm) curve 
    %

        % Set initial EC range (S/m)
        len_r= 5000;
        crv_EC_r= linspace(0.0001, 5, len_r);

        % Get responses (ppm)
        for ii= 1:length(crv_EC_r)

            % Update conductivity of halfspace (S/m)
            M.con= crv_EC_r(ii);

            % Calculate forward response (ppm)
            [FWD_IP(ii) FWD_QP(ii)]= FDEM1DFWD_RC(S,M);                    % Forward model (github.com/dhanssens)

        end

        
    %
    % Get rECa (S/m)
    %
          
        % Check monotony of EC (S/m) - QP (ppm) curve
        if all(diff(FWD_QP) > 0) == 1
            
            % Interpolation (pchip)
            iECa= interp1(FWD_QP,crv_EC_r,rQPdata,'pchip');
                        
        else
            
                   
    %
    % Initial ECa (S/m) estimation
    %

            % Get arguments
            narginchk(2,3); 

            % No initial ECa (S/m) estimation included
            if nargin < 3

                for r=1:dim(1)

                    for c=1:dim(2)

                        if isnan(I.QPdata(r,c)) == 0

                            % Calculate inital ECa (S/m) estimation
                            % Amplitude based
                            E= 1/2 .* sqrt((FWD_QP-I.QPdata(r,c)).^2 + (FWD_IP-I.IPdata(r,c)).^2);
                            init(r,c)= crv_EC_r(E == min(E(:)));

                        else

                            % Set NaN
                            init(r,c)= NaN;

                        end
                    end
                end

                % Reshape
                init= reshape(init,1,[]);

            % Initial ECa (S/m) estimation included
            else 

                init= repmat(init,1,numel(I.QPdata));

            end
                                    
            for ii= 1:numel(I.QPdata)
               
                % Check for NaN values
                if isnan(rQPdata(ii)) == 0
                
                    % Calculate cross-sections
                    [mItx]= polyxpoly(crv_EC_r,FWD_QP,[crv_EC_r(1) crv_EC_r(end)],[rQPdata(ii) rQPdata(ii)]);

                    % Case no cross-sections
                    if isempty(mItx) == 1; mItx= NaN; end
                    
                    % Account for initial ECa estimation
                    [M, ind]= min(abs(mItx-init(ii)),[],1);

                    % Get iECa (S/m)
                    iECa(ii)= mItx(ind);

                    % Clear var
                    clear mItx;
                
                else
                    
                    % Get ECa (S/m)
                    iECa(ii)= NaN;
                    
                end
            end
        end 
        
        % Original shape
        rECa= reshape(iECa,dim(1),dim(2)) .* blank;
        
end

