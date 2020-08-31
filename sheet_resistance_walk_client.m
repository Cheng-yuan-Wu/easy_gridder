%% Initialize and clear the results cache for batch runs:
clear;
global gridResultsCells
gridResultsCells = cell(0);
%% Initialize sheet resistance experiment
% Holding material properties constant, increase sheet resistance by 
% factors of 10. Observe the corresponding reduction in optimal wire width
% and overall power reduction.

FILE = './recipes/griddler_val_5cm.csv';
recipe = readtable(FILE);

Vop = recipe.value(1);
Jop = recipe.value(2);
Pwire = recipe.value(3);
Psheet = recipe.value(4);  % We'll be changing this...
w_min = 1e-7;  % Also override this one.
AspectRatio = recipe.value(7);
L = recipe.value(8);
Pcontact = recipe.value(9); %Ohm-cm^2

% Set Design Constraints
wireLimits = [w_min 10]; %cm, Allowed dimensions for the wires
geometry = 'bars'; %From bars, hexes, squares, or triangles

sheetFun = @(n) 1; %No sheet model: 100% transmission assumed

% Optimizer setup - take a best guess at grid geometry
numTiers = 2;

smallEnd = 2e-4;  %cm. Custom value since we overwrote the wire minimum
bigEnd = 0.1;  %cm

%% Solve
iters = 2000;
sheet_resistances = 10.^linspace(0, 6, 20);
sheet_resistances = 10;

for Psheet = sheet_resistances
    gridSolver = makeMultiTierSolver (Vop, Jop, Pwire, Pcontact,...
                                      AspectRatio, L, Psheet, geometry, ...
                                      sheetFun, wireLimits);
    V = solverSampler(gridSolver, iters, numTiers, smallEnd, bigEnd);
    
    % get power value:
    display = true;
    saveResults = true;
    power_loss = gridSolver(V, display, saveResults);
    power = Jop * Vop * (1 - power_loss);
    
    disp('R: ' + string(Psheet) + ' w: ' + string(V(end, 1))...
        + ' P:' + string(1e3 * power))
end