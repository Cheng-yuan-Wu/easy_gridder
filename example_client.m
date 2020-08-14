% This script provides an example of finding optimal grid structures across
% a range of design constraints. We define a range of cell sizes and
% transparent conductor sheet resistances, then find optimal 2-tier
% straight-line grid geometries for each combination of size and
% resistance.

clear;
%% Results cache for batch runs. Examine this variable later to extract results:
global gridResultsCells
gridResultsCells = cell(0);

%% System Material Constants
% Set based on active layer, metal choice for grid lines, and grid contact
% measurements

mat = 'GaSb';
switch mat
    case 'aSi'
        Jop = 0.012; %A/cm2
        Vop = 0.8; %V
    case 'PVT'
        Jop = 0.020;
        Vop = 0.8;
    case 'OPV'
        Jop = 0.010;
        Vop = 0.5;
    case 'GaSb'
        Jop = 0.016;
        Vop = 0.6;
end

Pwire = 10^-5; %Ohm-cm
Pcontact = 5e-3; %Ohm-cm^2

%% Design Constants
AspectRatio = 0.5; %aspect ratio of grid. 0.5 means h = w/2

Ls = [10 100]; %cm widths of solar cell
    
wireLimits = [0.0001 10]; %cm, Allowed dimensions for the wires
% Psheets = [30 60 100 1000 10000]; %Ohms/sq of the underlying material
Psheets = logspace(log10(30.1), log10(10000), 10);
varySheet = true ; %Calculate appropriate extinction coeff of the sheet?
geometry = 'bars'; %Choose bars, hexes, squares, or triangles

%% Create Transfer Functions
% Sheet resistance:
if varySheet
    % Interpolate sheet extinction values from measured points.
%     x = [40 0.975; 20 0.96; 10 0.92; 5 0.855; 2 0.76]; % ITO Sheet
    x = [10000 1; 300 0.85; 168 0.83; 68 0.78; 55 0.72; 30 0.5]; % FTO Sheet
    y = x(:,2);
    x = x(:,1);
    sheetFun = @(n) interp1(x, y, n);
else
    sheetFun = @(n) 1; %No Extinction
end

%% Optimizer setup - take a best guess at grid geometry
numTiers = 2;  % tiers of grid lines - 2 is usually fine

% range of randomized starting measurements in cm:
smallEnd = 4*min(wireLimits(:));
bigEnd = 0.6;

%% Solve
iters = 200;  % number of random starts per experiment
for L = Ls
    for Psheet = Psheets
        gridSolver = makeMultiTierSolver (Vop, Jop, Pwire, Pcontact, AspectRatio,...
                                          L, Psheet, geometry, sheetFun, wireLimits);
        optimum_v = solverSampler(gridSolver, iters, numTiers, smallEnd, bigEnd);
        
        %% Parse Result
        display = true;
        saveResults = true;
        gridSolver(optimum_v, display, saveResults);
    end
end

if saveResults
    xlswrite([datestr(now(), 'HHMMSS') '.xlsx'], gridResultsCells);
end
