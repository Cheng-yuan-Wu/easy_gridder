function result = solverSampler(solver, numRuns, numTiers, vMin, vMax)
% Run solver numRuns times with numTiers measurements set randomly between
% vMin and vMax. Return the single vector v that yielded the smallest loss
% among all runs.

maxLoss = 1;  % 100% power loss is the starting baseline
result = NaN;
for i = 1:numRuns
    % Initialize v with a random monotonic value spread in [vMin, vmax]
    v = sort((exp(rand(2*numTiers,1))-1)/(exp(1)-1)*vMax, 'descend');
    v = max(v, vMin);
    v = flip(reshape(v,[],2), 2);
    
    % Attempt to optimize v. Some attempts will fail; naivly ignore errors.
    try
        opt_v = fminsearch(@(n) solver(n, false, false), v,...
                        struct('TolFun', 1e-3,...
                               'MaxFunEvals', 300*6,...
                               'Display', 'off'));
%         opt_v = fmincon(@(n) solver(n, false, false), v, -ones(1, numel(v)), -vMin);

        loss = solver(opt_v, false, false);
        if abs(loss) < maxLoss
            result = opt_v;
        end
    catch
    end
end
end