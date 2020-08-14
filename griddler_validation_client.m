clear;
%% Result cache for batch runs:
global gridResultsCells
gridResultsCells = cell(0);
%% System Material Constants
% Set based on active layer, metal choice for grid lines, and grid contact
% measurements

FILE = './recipes/griddler_val_5cm.csv';
recipe = readtable(FILE);

Vop = recipe.value(1);
Jop = recipe.value(2);
Pwire = recipe.value(3);
Psheet = recipe.value(4);
w_min = recipe.value(5);
AspectRatio = recipe.value(7);
L = recipe.value(8);

Pcontact = recipe.value(9); %Ohm-cm^2

%% Set Design Constraints

wireLimits = [w_min 10]; %cm, Allowed dimensions for the wires
geometry = 'bars'; %From bars, hexes, squares, or triangles

sheetFun = @(n) 1; %No sheet model: 100% transmission assumed

%% Optimizer setup - take a best guess at grid geometry
numTiers = 2;

smallEnd = 2*min(wireLimits(:));  %cm
bigEnd = 0.1;  %cm


%% Solve
iters = 500;
gridSolver = makeMultiTierSolver (Vop, Jop, Pwire, Pcontact, AspectRatio,...
                                  L, Psheet, geometry, sheetFun, wireLimits);
optimum_v = solverSampler(gridSolver, iters, numTiers, smallEnd, bigEnd);

disp('Best design identified:')
disp('Width[cm]  Pitch[cm]')
disp(optimum_v)

%% Parse Result
display = true;
saveResults = true;
power_loss = gridSolver(optimum_v, display, saveResults);
power = Jop * Vop * (1 - power_loss);

disp('Power output [mW/cm2]:  ' + string(1e3 * power))
if saveResults
    xlswrite([datestr(now(), 'HHMMSS') '.xlsx'], gridResultsCells);
end

%% Calculate the corresponding Griddler simulation values
% Use single bus bar, single probe point, square wafer, dual print, and no
% finger end joining or gap
% 
% Diode characteristics:
% 1-sun current of 25 mA/cm2
% J01 450 fA/cm2
% J02 10 nA/cm2

params = {'wafer length', L;  %cm
    'wafer width', optimum_v(1, 2);  %cm
    'busbar width', optimum_v(1, 1) * 10;  %mm
    'no fingers', wafer_width/optimum_v(2, 2);
    'finger width', optimum_v(2, 1) * 1e4;  %um
    'finger mohm/sq', 1e3 * Pwire / (optimum_v(2, 1) * AspectRatio); %mohm/sq
    'bus mohm/sq', 1e3 * Pwire / (optimum_v(1, 1) * AspectRatio); %mohm/sq
    'finger contact res', 1e3 * Pcontact;  %mohm-cm2
    'layer sheet res', Psheet}  %#ok<*NOPTS> %ohm/sq

%% Manually-determined geometry and operating conditions
Vmp = 0.560;
Jmp = 0.0196;
target_power = 10.192e-3;
v = [0.0832 0.6669;
     0.0053 0.6669/8];

% Equivalent sim parameters
params = {'wafer length', L;  %cm
    'wafer width', v(1, 2);  %cm
    'busbar width', v(1, 1) * 10;  %mm
    'no fingers', wafer_width/v(2, 2);
    'finger width', v(2, 1) * 1e4;  %um
    'finger mohm/sq', 1e3 * Pwire / (v(2, 1) * AspectRatio); %mohm/sq
    'bus mohm/sq', 1e3 * Pwire / (v(1, 1) * AspectRatio); %mohm/sq
    'finger contact res', 1e3 * Pcontact;  %mohm-cm2
    'layer sheet res', Psheet}  %ohm/sq

power_loss = gridSolver(v, true, false);
power = Jmp * Vmp * (1 - power_loss);
factor_adjustment = target_power/power;
disp('Power output [mW/cm2]:  ' + string(1e3 * power))
disp('Adjusted output [mW/cm2]:  ' + string(1e3 * power * factor_adjustment))

%% Sweep