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
Psheet = recipe.value(4);
w_min = recipe.value(5);
AspectRatio = recipe.value(7);
L = recipe.value(8);

Pcontact = recipe.value(9); %Ohm-cm^2