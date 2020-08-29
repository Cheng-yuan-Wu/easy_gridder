%% Initialize and clear the results cache for batch runs:
clear;
global gridResultsCells
gridResultsCells = cell(0);
%% System Material Constants
% Set based on active layer, metal choice for grid lines, and grid contact
% measurements. For Griddler validation, we take J_max_power from
% simulation and Vop is the highest sheet voltage in the FEA.

FILE = './recipes/griddler_val_20cm.csv';
recipe = readtable(FILE);

Vop = recipe.value(1);
Jop = recipe.value(2);
Pwire = recipe.value(3);
Psheet = recipe.value(4);
w_min = recipe.value(5);
AspectRatio = recipe.value(7);
L = recipe.value(8);

Pcontact = recipe.value(9); %Ohm-cm^2

% Set Design Constraints
wireLimits = [w_min 10]; %cm, Allowed dimensions for the wires
geometry = 'bars'; %From bars, hexes, squares, or triangles

sheetFun = @(n) 1; %No sheet model: 100% transmission assumed

% Optimizer setup - take a best guess at grid geometry
numTiers = 2;

smallEnd = 2*min(wireLimits(:));  %cm
bigEnd = 0.1;  %cm


%% Solve
iters = 500;
gridSolver = makeMultiTierSolver (Vop, Jop, Pwire, Pcontact, AspectRatio,...
                                  L, Psheet, geometry, sheetFun, wireLimits);
% Carry this measurement vector through the rest of the problem
V = solverSampler(gridSolver, iters, numTiers, smallEnd, bigEnd);

disp('Best design identified:')
disp('Width[cm]  Pitch[cm]')
disp(V)

% Parse Result
display = true;
saveResults = true;
power_loss = gridSolver(V, display, saveResults);
power = Jop * Vop * (1 - power_loss);

disp('Power loss fraction: ' + string(power_loss))
disp('Power output [mW/cm2]:  ' + string(1e3 * power))
% if saveResults
%     xlswrite([datestr(now(), 'HHMMSS') '.xlsx'], gridResultsCells);
% end

% Calculate the Griddler simulation values to recreate this design.
% Use single bus bar, single probe point, square wafer, dual print, and no
% finger end joining or gap
% 
% Typical mono-Si cell characteristics:
% 1-sun current of 25 mA/cm2
% J01 450 fA/cm2
% J02 10 nA/cm2

params = {'wafer length', L;  %cm
    'wafer width', V(1, 2);  %cm
    'busbar width', V(1, 1) * 10;  %mm
    'no fingers', V(1, 2)/V(2, 2);
    'finger width', V(2, 1) * 1e4;  %um
    'finger mohm/sq', 1e3 * Pwire / (V(2, 1) * AspectRatio); %mohm/sq
    'bus mohm/sq', 1e3 * Pwire / (V(1, 1) * AspectRatio); %mohm/sq
    'finger contact res', 1e3 * Pcontact;  %mohm-cm2
    'layer sheet res', Psheet}  %#ok<*NOPTS> %ohm/sq

%% Manually-determined geometry and operating conditions
num_fingers = 19   % Take a rounded value from the optimization run
Vmp = 0.499;   % V max power from griddler
Jmp = 0.01683;  % J max power from griddler
target_power = Vmp * Jmp;  %mohm/cm2
V(2, 2) = V(1, 2) / num_fingers
V0 = V;   % Save this optimal geometry - return at the start of each

% Equivalent sim parameters
params = {'wafer length', L;  %cm
    'wafer width', V(1, 2);  %cm
    'busbar width', V(1, 1) * 10;  %mm
    'no fingers', V(1, 2)/V(2, 2);
    'finger width', V(2, 1) * 1e4;  %um
    'finger mohm/sq', 1e3 * Pwire / (V(2, 1) * AspectRatio); %mohm/sq
    'bus mohm/sq', 1e3 * Pwire / (V(1, 1) * AspectRatio); %mohm/sq
    'finger contact res', 1e3 * Pcontact;  %mohm-cm2
    'layer sheet res', Psheet}  %ohm/sq

power_loss = gridSolver(V, true, false);
power = Jmp * Vmp * (1 - power_loss);
factor_adjustment = target_power/power;  % Carry this factor to the next section
disp('Power output [mW/cm2]:  ' + string(1e3 * power))
disp('Tuned output [mW/cm2]:  ' + string(1e3 * power * factor_adjustment))

%% Sweep line pitches
V = V0;
n = [10 12 14 16 18 19 20 25 30 35 40];
pitches = V(1,2)./n

for p = pitches
    disp('>>>>>>>>>>>>')
    disp('Pitch: ' + string(p))
    V(2,2) = p;
    power_loss = gridSolver(V, false, false);
    power = Jmp * Vmp * (1 - power_loss) * factor_adjustment;
    disp('power:  ' + string(power * 1e3))
    
    % Equivalent sim parameters
    params = {'wafer length', L;  %cm
        'wafer width', V(1, 2);  %cm
        'busbar width', V(1, 1) * 10;  %mm
        'no fingers', V(1, 2)/V(2, 2);
        'finger width', V(2, 1) * 1e4;  %um
        'finger mohm/sq', 1e3 * Pwire / (V(2, 1) * AspectRatio); %mohm/sq
        'bus mohm/sq', 1e3 * Pwire / (V(1, 1) * AspectRatio); %mohm/sq
        'finger contact res', 1e3 * Pcontact;  %mohm-cm2
        'layer sheet res', Psheet}  %ohm/sq
    disp('>>>>>>>>>>>>')
end

%% High-resolution pitch sweep
V = V0;
n = linspace(9.5, 41);
pitches = V(1,2)./n;
result = zeros(100, 2);
result(:, 1) = pitches;
i = 1;

for p = pitches
    V(2,2) = p;
    power_loss = gridSolver(V, false, false);
    power = Jmp * Vmp * (1 - power_loss) * factor_adjustment;
    result(i, 2) = power;
    i = i + 1;
end

%% Sweep line widths
V = V0;
widths = [70 90 110 130 145 147 150 170 200 250 300] .* 1e-4;

for w = widths
    disp('>>>>>>>>>>>>')
    disp('Width: ' + string(w) + ' cm')
    V(2, 1) = w;
    power_loss = gridSolver(V, false, false);
    power = Jmp * Vmp * (1 - power_loss) * factor_adjustment;
    disp('power:  ' + string(power * 1e3))
    
    % Equivalent sim parameters
    params = {
        'finger width', V(2, 1) * 1e4;  %um
        'finger mohm/sq', 1e3 * Pwire / (V(2, 1) * AspectRatio); %mohm/sq
        'wafer length', L;  %cm
        'wafer width', V(1, 2);  %cm
        'busbar width', V(1, 1) * 10;  %mm
        'no fingers', V(1, 2)/V(2, 2);
        'bus mohm/sq', 1e3 * Pwire / (V(1, 1) * AspectRatio); %mohm/sq
        'finger contact res', 1e3 * Pcontact;  %mohm-cm2
        'layer sheet res', Psheet}  %ohm/sq
    disp('>>>>>>>>>>>>')
end

%% High-resolution width sweep
V = V0;
widths = linspace(65, 310) .* 1e-4;
result = zeros(100, 2);
result(:, 1) = widths;
i = 1;

for w = widths
    V(2, 1) = w;
    power_loss = gridSolver(V, false, false);
    power = Jmp * Vmp * (1 - power_loss) * factor_adjustment;
    result(i, 2) = power;
    i = i + 1;
end