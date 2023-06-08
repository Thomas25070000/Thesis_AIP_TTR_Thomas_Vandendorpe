function runningtimes = generateRunningTimes(blocksections, settings)

vyellow = settings.infrastructure.vyellow / 3.6;

buffer_factor = 1 + settings.TT.buffer;

%% Start with the disruption for L trains.
vmax = settings.trains.speed.R;
acc = settings.trains.acc.R;
decc = -settings.trains.decc.R;

direction = 12;
blockSpeed = min(blocksections.vmax_disruption_dir1, vmax) / 3.6;     % Include conversion to m/s
runningtimes.L12.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L12.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L12.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 13;
blockSpeed = min(blocksections.vmax_disruption_dir1, vmax) / 3.6;     % Include conversion to m/s
runningtimes.L13.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L13.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L13.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
blockSpeed = min(blocksections.vmax_disruption_dir0, vmax) / 3.6;
runningtimes.L02.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L02.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L02.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 3;
blockSpeed = min(blocksections.vmax_disruption_dir0, vmax) / 3.6;
runningtimes.L03.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L03.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L03.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);



%% L in normal case
blockSpeed = min(blocksections.vmax, vmax) / 3.6;

direction = 12;
runningtimes.L12.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L12.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L12.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 13;
runningtimes.L13.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L13.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L13.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
runningtimes.L02.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L02.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L02.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 3;
runningtimes.L03.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.L03.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.L03.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);




%% IC in disruption
vmax = settings.trains.speed.IC;
acc = settings.trains.acc.IC;
decc = -settings.trains.decc.IC;


direction = 12;
blockSpeed = min(blocksections.vmax_disruption_dir1, vmax) / 3.6;     % Include conversion to m/s
runningtimes.IC12.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC12.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC12.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 13;
blockSpeed = min(blocksections.vmax_disruption_dir1, vmax) / 3.6;     % Include conversion to m/s
runningtimes.IC13.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC13.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC13.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
blockSpeed = min(blocksections.vmax_disruption_dir0, vmax) / 3.6;     % Include conversion to m/s
runningtimes.IC02.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC02.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC02.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 3;
blockSpeed = min(blocksections.vmax_disruption_dir0, vmax) / 3.6;     % Include conversion to m/s
runningtimes.IC03.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC03.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC03.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

%% IC in normal case
blockSpeed = min(blocksections.vmax, vmax) / 3.6;

direction = 12;
runningtimes.IC12.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC12.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC12.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 13;
runningtimes.IC13.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC13.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC13.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
runningtimes.IC02.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC02.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC02.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 3;
runningtimes.IC03.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
runningtimes.IC03.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
runningtimes.IC03.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);


end

%% Regular
function RT = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction)

if direction == 2 | direction == 3
    blockSpeed(1:end) = blockSpeed(end:-1:1);
    blocksections = sortrows(blocksections,'id','descend');
end

% Generate the basic one, running at full speeds all the time.
Nblocks = size(blocksections,1);
RT = zeros(Nblocks,1);

for bb = 1:Nblocks-1
    if blockSpeed(bb) == blockSpeed(bb+1) || blockSpeed(bb) < blockSpeed(bb+1)
        % We can continue at the same speed.
        try
            if blockSpeed(bb) > blockSpeed(bb-1)
                % We could accelerate from the previous one!
                ta = (blockSpeed(bb) - blockSpeed(bb-1)) / acc;
                xa = (blockSpeed(bb) - blockSpeed(bb-1))^2 / 2 / acc;
                RT(bb) = (blocksections.length(bb) - xa) / blockSpeed(bb) + ta;
            else
                RT(bb) = blocksections.length(bb) / blockSpeed(bb);
            end
        catch
            RT(bb) = blocksections.length(bb) / blockSpeed(bb);
        end
    elseif blockSpeed(bb) > blockSpeed(bb+1)
        % We have to deccelerate to the indicated speed from the start of
        % this section on!
        td = (blockSpeed(bb) - blockSpeed(bb+1)) / decc;
        xd = (blockSpeed(bb) - blockSpeed(bb+1))^2 / 2 / decc;
        RT(bb) = (blocksections.length(bb) - xd) / blockSpeed(bb+1) + td;
    else
        try
            if blockSpeed(bb) > blockSpeed(bb-1)
                % We could accelerate!
                ta = (blockSpeed(bb) - blockSpeed(bb-1)) / acc;
                xa = (blockSpeed(bb) - blockSpeed(bb-1))^2 / 2 / acc;
                RT(bb) = (blocksections.length(bb) - xa) / blockSpeed(bb) + ta;
            elseif blockSpeed(bb) == blockSpeed(bb-1)
                % Run at the max speed.
                RT(bb) = blocksections.length(bb) / blockSpeed(bb);
            end
        end
    end
end



% BUT!: if D contains only one block section, it's max. speed will depend
% on the characteristics at the switches, i.e. the blocks before and after
% it.
indexD = find(strcmp(blocksections.type,'D'));
nrBlocksD = length(indexD);
if nrBlocksD == 1
    bb = indexD;
    vS1 = blockSpeed(indexD-1);
    vD = blockSpeed(indexD);
    vS2 = blockSpeed(indexD+1);
	% Compare with the switch before
    if (vS1 < vD && vS2 < vD) || (vS1 < vD && vS2 == vD)
        % Change speed from vS1 to vS2
        if vS1 <= vS2
            % Accelerate!
            ta = (vS2 - vS1) / acc;
            xa = (vS2 - vS1)^2 / 2 / acc;
            RT(bb) = (blocksections.length(bb) - xa) / vS2 + ta;
        else
            % Deccelerate!
            td = (vS1 - vS2) / decc;
            xd = (vS1 - vS2)^2 / 2 / decc;
            RT(bb) = (blocksections.length(bb) - xd) / vS2 + td;
        end
    elseif (vS1 < vD && vS2 > vD)
        % Accelerate!
        ta = (vD - vS1) / acc;
        xa = (vD - vS1)^2 / 2 / acc;
        RT(bb) = (blocksections.length(bb) - xa) / vD + ta;
    elseif (vS1 >= vD && vS2 > vD) || (vS1 == vD && vS2 == vD)
        % Stay at same speed!
        RT(bb) = blocksections.length(bb) / blockSpeed(bb);
    elseif (vS1 >= vD && vS2 < vD)
        % Deccelerate!
        td = (vD - vS2) / decc;
        xd = (vD - vS2)^2 / 2 / decc;
        RT(bb) = (blocksections.length(bb) - xd) / vS2 + td;
    end
end

% For the last section.
if blockSpeed(end) == blockSpeed(end-1) || blockSpeed(end) < blockSpeed(end-1)
    RT(end) = blocksections.length(end) / blockSpeed(end);
elseif blockSpeed(end) > blockSpeed(end-1)
    ta = (blockSpeed(end) - blockSpeed(end-1)) / acc;
    xa = (blockSpeed(end) - blockSpeed(end-1))^2 / 2 / acc;
    RT(end) = (blocksections.length(end) - xa) / blockSpeed(end) + ta;
end

if direction == 2 | direction == 3
    RT(1:end) = RT(end:-1:1);
end

RT = round((RT * buffer_factor)');

if ~all(RT)
    disp('Problem');
    disp(RT);
    pause;
end

end

%% Braking
function RT = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction)

if direction == 2 | direction == 3
    blockSpeed(1:end) = blockSpeed(end:-1:1);
    blocksections = sortrows(blocksections,'id','descend');
end

% Generate the basic one, running at full speeds all the time.
Nblocks = size(blocksections,1);
RT = zeros(Nblocks,1);

% For the first block, assuming incoming at full speed.
RT(1) = blockSpeed(1) / abs(decc) + blocksections.length(1) / vyellow;
RT(2:end) = blockSpeed(1:end-1) / abs(decc) + blocksections.length(2:end) / vyellow;

if direction == 0
    RT(1:end) = RT(end:-1:1);
end

RT = round((RT * buffer_factor)');

if ~all(RT)
    disp('Problem');
    disp(RT);
    pause;
end

end

%% Accelerating
function RT = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor, direction)

if direction == 2 | direction == 3
    blockSpeed(1:end) = blockSpeed(end:-1:1);
    blocksections = sortrows(blocksections,'id','descend');
end

% Generate the basic one, running at full speeds all the time.
Nblocks = size(blocksections,1);
RT = zeros(Nblocks,1);

% For the first block, assuming incoming at full speed.
targetSpeed = min(blockSpeed(2:end), blockSpeed(1:end-1));
targetSpeed(Nblocks) = blockSpeed(end);
RT = targetSpeed / acc + blocksections.length ./ targetSpeed;

if direction == 2 | direction == 3
    RT(1:end) = RT(end:-1:1);
end

RT = round((RT * buffer_factor)');

if ~all(RT)
    disp('Problem');
    disp(RT);
    pause;
end

end
