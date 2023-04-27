function [fitobj, goodness, output] = zesespolfit(polcurvestruct, varargin)
%ZESESPOLFIT Fits empirical polarization curve parameters to ZESES
%   fuel cell data. polcurvestruct must be a polarization curve imported with
%   fcpolread_v2.m or fcpolread_G60II.m and further tweaked with fitAllPatData.m 
%   or similar so it has an extractedConstants feature.
%   [FO, G, O] = zesespolfit(polcurvestruct) creates a fit object, FO, that
%   encapsulates the result of fitting the fuel cell polarization
%   curve given in the structure "polcurvestruct" obtained from calling
%   fcpolread_v2.m on output from the SINTEF fuel cell test stands.
%   Mike Gerhardt originally got this code from Simon Clark & tinkered
%   heavily with it. It came from the AlkalineElectrolyzerDataFitting
%   repository and was originally called wepolfit.m.
%   [FO, G, O] = zesespolfit(I, uI, U, uU, 'param1', value1, 'param2', value2,
%   ...) creates a fit object that considers the user-defined
%   parameter-value pairs defined by e.g. ...,'param1', value1, ... The
%   user can choose from the following parameter value pairs:
%
%   'temp'      Operational temperature, [K]. Default = 323 K. 
%   'Er'        Reversible cell voltage, [V]. Default = 1.229 V.
%   'Rohm'      Ohmic cell resistance, [Ohm cm2]. Default = 0.002 Ohm cm2. 
%   'OERi0'     OER exchange current density, [A/cm2]. Default = 1e-6.
%   'HERi0'     HER exchange current dentisy, [A/cm2]. Default = 1e-6.
%   'i0'        Lumped exchange current density, [A/cm2]. Default = 1e-6.
%   'OERalpha'  OER charge transfer coefficient, [-]. Default = 0.41.
%   'HERalpha'  HER charge transfer coefficient, [-]. Default = 0.37.
%   'n'         Electrons transferred, [-]. Default  = 2.
%   'A'         Tafel slope, [V]. Defualt = 0.040 V.
%   'kinetics'  Kinetic model |'Tafel'|'Butler-Volmer'|. Default = 'Tafel'.
%   'display'   Display output |'true'|'false'|. Default = 'true'.

%% Input checks
I = polcurvestruct.If./1000;
U = polcurvestruct.Uf;
knowns = polcurvestruct.extractedConstants;

% Check that I and U have the same shape
if ~isequal(size(I),size(U))
    error('check_iu: I & U must have same dimensions');
end


%% Initialization
% Initialize default fields array
fields = {  'temp',     353, ...        % Temperature (K)
            'Er',       1.20336, ...      % Voltage (V)
            'Rohm',     knowns.R__M__cm2_./1000, ...      % Ohmic Resistance (Ohm cm2)
            'i0',       knowns.i0_A_cm2_, ...       % Lumped exchange current density (A/cm2)
            'n',        2, ...          % Number of electrons
            'b',        knowns.b_mV_dec_./1000, ...      % Tafel slope (V/dec)
            'in',       knowns.in_mA_cm2_./1000, ...       % H2 cross-over current density (A/cm2)
            'xi',       0.1,...         % Mass transfer/limiting current prefactor (V)
            'method', 'consti0',...    % fit method
            'display',  'true' };      % Display output

% Initialize fitobjput structure
fitobj = struct(fields{:});
    
% Switch user-defined fields
for i = 1:2:length(varargin)
    varn    = varargin{i};
    varval  = varargin{i+1};
    switch varn
        case 'temp'
            fitobj.temp        = varval;
        case 'Er'
            fitobj.Er          = varval;
        case 'Rohm'
            fitobj.Rohm        = varval;
        case 'i0'
            fitobj.i0          = varval;
        case 'n'
            fitobj.n           = varval;
        case 'b'
            fitobj.b           = varval;
        case 'in'
            fitobj.in          = varval;
        case 'xi'
            fitobj.xi          = varval;
        case 'method'
            fitobj.method      = varval;
        case 'kinetics'
            fitobj.kinetics    = varval;
        case 'display'
            fitobj.display     = varval;
    end
end
    
%% Define Physical Constants
R = 8.31446261815324;   % J mol^-1 K^-1
F = 96485.3329;         % A s mol^-1

%% Define Shorter Variable Names
% Note all variables below are treated as "known" (not fitting parameters)
T           = fitobj.temp;
Er          = fitobj.Er; 
% OERalpha    = fitobj.OERalpha;
% OERi0       = fitobj.OERi0;
% HERalpha    = fitobj.HERalpha;
% HERi0       = fitobj.HERi0;
i0          = fitobj.i0; 
Rohm        = fitobj.Rohm; 
in          = fitobj.in;
n           = fitobj.n;
b           = fitobj.b; 

%% Fit Experimental Data
x = I;
y = U;



% Define the modelling method for the fit


% Define the fit type
switch fitobj.method
    case 'consti0'
        FT = fittype(@(xi,iL,x) Er - b.a*log10((x+in)./i0) - Rohm*x + xi.*log10(1-x./iL),...
            'independent',{'x'}, ...
            'coefficients',{'xi','iL'});
        
        % FT = fittype('Er + b.*log10((x+in)./i0) + Rohm*x + xi.*log10(1-x/iL)', ...
        %     'dependent',{'y'},'independent',{'x'},...
        %     'coefficients',{'Er', 'i0', 'b','Rohm', 'xi', 'iL'}, ...
        %     'problem', {'i0'});
        
        % Define fit options, start point, bounds, and tolerances
        options             = fitoptions(FT);
        options.StartPoint  = [0.1,10];
        options.Lower       = [0, 0];
        %options.Upper       = [100, 30];
        options.TolX        = 1e-6;
        options.TolFun      = 1e-6;
        
        % Calculate the fit
        [FO, G, O] = fit(x,y,FT,options);
        
        % Store the fitted values in the fit object
        fitobj.xi      = FO.xi;
        fitobj.iL      = FO.iL;
        
        
        % Store the fitted values in the fit variable object for display
        fitvar.xi      = FO.xi;
        fitvar.iL      = FO.iL;
        
  
    case 'vari0'
        FT = fittype(@(xi,iL,i0v,x) Er - b.*log10((x+in)./i0v) - Rohm*x + xi.*log10(1-x./iL),...
            'independent',{'x'}, ...
            'coefficients',{'xi','iL','i0v'});
        
        % Define fit options, start point, bounds, and tolerances
        options             = fitoptions(FT);
        options.StartPoint  = [0.1,10, fitobj.i0];
        options.Lower       = [0, 0, 0];
        %options.Upper       = [100, 30];
        options.TolX        = 1e-6;
        options.TolFun      = 1e-6;
        
        % Calculate the fit
        [FO, G, O] = fit(x,y,FT,options);
        
        % Store the fitted values in the fit object
        fitobj.xi      = FO.xi;
        fitobj.iL      = FO.iL;
        fitobj.i0v     = FO.i0v;
        
        % Store the fitted values in the fit variable object for display
        fitvar.xi      = FO.xi;
        fitvar.iL      = FO.iL;
        fitvar.i0v     = FO.i0v;
  
    case 'consti0-constxi'
        
        xi = fitobj.xi;
        FT = fittype(@(iL,x) Er - b.*log10((x+in)./i0) - Rohm*x + xi.*log10(1-x./iL),...
            'independent',{'x'}, ...
            'coefficients',{'iL'});
        
        
        % Define fit options, start point, bounds, and tolerances
        options             = fitoptions(FT);
        options.StartPoint  = [2];
        options.Lower       = [1.5];
        options.Upper       = [10];
        options.TolX        = 1e-6;
        options.TolFun      = 1e-6;
        
        % Calculate the fit
        [FO, G, O] = fit(x,y,FT,options);
        
        % Store the fitted values in the fit object
        fitobj.iL      = FO.iL;
        
        % Store the fitted values in the fit variable object for display
        fitvar.iL      = FO.iL;
            
end

% Store the fit goodness
goodness = G;

% Store the fit output
output = O;

% Calculate the interpolated values for the output comparison
fitobj.I = linspace(0,max(I),1000)';
fitobj.U = fitobj.Er - fitobj.b.*log10((fitobj.I+fitobj.in)./fitobj.i0) - fitobj.Rohm*fitobj.I + fitobj.xi.*log10(1-fitobj.I./fitobj.iL);

%% Display results if fitobjput is selected (default = true)
if strcmpi(fitobj.display,'true') == 1
    %make title string
    titlestr = string(polcurvestruct.extractedConstants.Time_hr_) + ' hrs, ' + ...
        string(polcurvestruct.extractedConstants.applied_current) + ' mA/cm^2';
    
    figure
    plot(fitobj.I, fitobj.U, 'LineWidth', 3), hold on
    plot(I, U, 'o', 'LineWidth', 2, 'MarkerSize', 7), hold off
    
    % Format the axes and plot area
    xlabel('Current Density  /  A \cdot cm^{-2}')
    xlim([min(fitobj.I), max(fitobj.I)])
    ylabel('Voltage  /  V')
    legend('Fit', 'Data', 'location', 'northwest')
    title(titlestr)
end

disp(' ')
disp('The following variables were fit: ')
disp(' ')
disp(fitvar)

end

