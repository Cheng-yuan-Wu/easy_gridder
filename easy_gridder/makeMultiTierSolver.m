function mySolver = makeMultiTierSolver(Vop, Jop, Pwire, Pcontact, AspectRatio,...
                                        L, Psheet, geometry, sheetFun, wireLimits)
    function powerLoss = gridFunction(v, display, saveResults)
        % v is a vector containing wire width and spacing values for some
        % number of current collection lattices:
        % [ Width_1 Pitch_1; ... %<we would call this tier the "busbars"
        %   Width_2 Pitch_2; ...
        %      :       :
        %   Width_n Pitch_n]  %<this tier is the only one that interacts
        %                       with the underlying sheet material
                                        
        % If any wire widths are out of bounds, confine them to preset size
        % limits:
        v(:,1) = confine(v(:,1), wireLimits);
        
%         v(1,2) = L;  % Lock in dimension values here
        v = abs(vertcat([inf L], v));   % Slap that top level L right on thar
                                        % disallow any negative numbers 
                                        % (they're non-physical, but may result from optimizer)

        % Hex and Square calculations are equivalent:
        if strcmp(geometry, 'squares')
            geometry = 'hexes';
        end
        
        % Assign loss calculation functions based on selected geometry
        switch geometry
            case 'bars'
                shadowFun =  @barShadow;
                contactFun = @barContact;
                lineFun =    @barLine;
            case 'hexes'
                shadowFun =  @hexShadow;
                contactFun = @hexContact;
                lineFun =    @hexLine;
            case 'triangles'
                shadowFun =  @triShadow;
                contactFun = @triContact;
                lineFun =    @triLine;
            otherwise
                error('Not a valid grid geometry');
        end
        
        % Compute individual loss contributions:
        gridLosses = zeros(size(v));
        for i = 2:size(v,1)
            gridLosses(i,1) = shadowFun(v(i,1), v(i,2)); %Shadow loss f(w, l)
            gridLosses(i,2) = lineFun(v(i,1), v(i,2), v(i-1,2)); %Linegrid loss f(w, l, L)
        end
        gridLosses = gridLosses(2:end,:); %Truncate empty top row
        
        % Bottom Layer Losses
        sheetExtinction = 1 - sheetFun(Psheet);
        sheetDrop = Jop/Vop * 1/12 * Psheet * v(end,end)^2;
        contact = contactFun(v(end,1), v(end,end));
        
        % Sum total power loss
        powerLoss = 1 - prod(1-gridLosses(:))*...
                        (1-sheetExtinction)*...
                        (1-sheetDrop)*...
                        (1-contact);
        
        if display || saveResults % Prep outputs
            gridStats = cell(1,5*(size(v,1)-1));
            statLabels = cell(1,5*(size(v,1)-1));
            for i = 1:size(v,1)-1
                gridStats{5*(i-1)+1} = v(i+1,1);% Width
                statLabels{5*(i-1)+1} = ['Width ' num2str(i) ' [cm]'];% Width
                
                gridStats{5*(i-1)+2} = v(i+1,2);% Spacing
                statLabels{5*(i-1)+2} = ['Space ' num2str(i) ' [cm]'];% Spacing
                
                gridStats{5*(i-1)+3} = 100*gridLosses(i,1);% Shadow
                statLabels{5*(i-1)+3} = ['Shadow ' num2str(i) ' [%]'];% Shadow
                
                gridStats{5*(i-1)+4} = 100*gridLosses(i,2);% Line
                statLabels{5*(i-1)+4} = ['Line ' num2str(i) ' [%]'];% Line
                
                gridStats{5*(i-1)+5} = Pwire*v(i+1,2)/(v(i+1,1)^2*AspectRatio);% Ohm/sq
                statLabels{5*(i-1)+5} = ['Ohm/sq ' num2str(i)];% Ohm/sq
            end
            gridStats = horzcat({Psheet, AspectRatio, L, wireLimits(1), wireLimits(2), 100*powerLoss},...
                gridStats,...
                {100*sheetDrop, 100*contact, 100*sheetExtinction});
            statLabels = horzcat({'Psheet  ', 'Aspect Ratio', 'L [cm]  ','wire min [cm]', 'wire max [cm]', 'Power Loss [%]'},...
                statLabels,...
                {'Sheet Drop %', 'Contact Loss %', 'Sheet Extinction %'}); 
        end
        
        if display
        % Display results in the terminal
            fprintf('\n');
            fprintf('with %d tiers\n', size(v,1)-1);
            for i = 1:size(statLabels,2)
%             for i = 6 %To only display power loss
                fprintf('%s\t% f\n',statLabels{i},gridStats{i});
            end
            fprintf('\n');
        end
        
        if saveResults
        % Append results to a global results variable
            global gridResultsCells %#ok<TLEV>
            [tableSize,~] = size(gridResultsCells);
            if tableSize < 1
                gridResultsCells = statLabels;
            end
            gridResultsCells = vertcat(gridResultsCells, gridStats);
        end
    end
    mySolver = @gridFunction;
    
    function [loss] = barShadow(w, l)
        loss = w/l;
    end
    function [loss] = barContact(w, l)
        loss = Jop/Vop * Pcontact * l/w;
    end
    function [loss] = barLine(w, l, L)
        loss = Jop/Vop * 1/12 * (Pwire*l)/(w^2*AspectRatio) * L^2;
    end
% TODO implement these:
%     case 'hexes'
%         shadowFun =  @hexShadow;
%         contactFun = @hexContact;
%         lineFun =    @hexLine;
%     case 'triangles'
%         shadowFun =  @triShadow;
%         contactFun = @triContact;
%         lineFun =    @triLine;
end

function x = confine(x, limits)
    for i = 1:size(x,1)
        if x(i,:) < limits(1)
            x(i,:) = limits(1);
        elseif x(i,:) > limits(2)
            x(i,:) = limits(2);
        end
    end
end