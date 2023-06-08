function [runningtimes, clearingtimes] = generateRunningAndClearingTimes(blocksections, settings)

vyellow = settings.infrastructure.vyellow / 3.6;

buffer_factor = 1 + settings.TT.buffer;

%% Start with the disruption for L trains.
vmax = settings.trains.speed.R;
acc = settings.trains.acc.R;
decc = -settings.trains.decc.R;
length = settings.trains.length.R;

direction = 12;
blockSpeed = min(blocksections.vmax_disruption_dir12, vmax) / 3.6;     % Include conversion to m/s
[runningtimes.L12.disrupted, clearingtimes.L12.disrupted] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L12.disrupted = clearingtimes.L12.disrupted + settings.TT.blocktimes.releaseT;

direction = 13;
blockSpeed = min(blocksections.vmax_disruption_dir13, vmax) / 3.6;     % Include conversion to m/s
[runningtimes.L13.disrupted, clearingtimes.L13.disrupted] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L13.disrupted = clearingtimes.L13.disrupted + settings.TT.blocktimes.releaseT;

% runningtimes.L1.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.L1.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
blockSpeed = min(blocksections.vmax_disruption_dir02, vmax) / 3.6;
[runningtimes.L02.disrupted, clearingtimes.L02.disrupted] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L02.disrupted = clearingtimes.L02.disrupted + settings.TT.blocktimes.releaseT;

direction = 3;
blockSpeed = min(blocksections.vmax_disruption_dir03, vmax) / 3.6;
[runningtimes.L03.disrupted, clearingtimes.L03.disrupted] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L03.disrupted = clearingtimes.L03.disrupted + settings.TT.blocktimes.releaseT;

% runningtimes.L0.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.L0.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);



%% L in normal case
blockSpeed = min(blocksections.vmax, vmax) / 3.6;

direction = 12;
[runningtimes.L12.regular, clearingtimes.L12.regular] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L12.regular = clearingtimes.L12.regular + settings.TT.blocktimes.releaseT; 
% runningtimes.L1.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
% runningtimes.L1.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.L1.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 13;
[runningtimes.L13.regular, clearingtimes.L13.regular] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L13.regular = clearingtimes.L13.regular + settings.TT.blocktimes.releaseT; 
% runningtimes.L1.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
% runningtimes.L1.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.L1.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
[runningtimes.L02.regular, clearingtimes.L02.regular] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L02.regular = clearingtimes.L02.regular + settings.TT.blocktimes.releaseT; 

direction = 3;
[runningtimes.L03.regular, clearingtimes.L03.regular] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.L03.regular = clearingtimes.L03.regular + settings.TT.blocktimes.releaseT;
% runningtimes.L0.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
% runningtimes.L0.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.L0.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);




%% IC in disruption
vmax = settings.trains.speed.IC;
acc = settings.trains.acc.IC;
decc = -settings.trains.decc.IC;


direction = 12;
blockSpeed = min(blocksections.vmax_disruption_dir12, vmax) / 3.6;     % Include conversion to m/s
[runningtimes.IC12.disrupted , clearingtimes.IC12.disrupted ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC12.disrupted = clearingtimes.IC12.disrupted + settings.TT.blocktimes.releaseT; 

direction = 13;
blockSpeed = min(blocksections.vmax_disruption_dir13, vmax) / 3.6;     % Include conversion to m/s
[runningtimes.IC13.disrupted , clearingtimes.IC13.disrupted ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC13.disrupted = clearingtimes.IC13.disrupted + settings.TT.blocktimes.releaseT; 
% runningtimes.IC1.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
% runningtimes.IC1.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.IC1.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
blockSpeed = min(blocksections.vmax_disruption_dir02, vmax) / 3.6;     % Include conversion to m/s
[runningtimes.IC02.disrupted , clearingtimes.IC02.disrupted ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC02.disrupted = clearingtimes.IC02.disrupted + settings.TT.blocktimes.releaseT; 

direction = 3;
blockSpeed = min(blocksections.vmax_disruption_dir03, vmax) / 3.6;     % Include conversion to m/s
[runningtimes.IC03.disrupted , clearingtimes.IC03.disrupted ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC03.disrupted = clearingtimes.IC03.disrupted + settings.TT.blocktimes.releaseT; 
% runningtimes.IC0.disrupted = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
% runningtimes.IC0.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.IC0.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);


%% IC in normal case
blockSpeed = min(blocksections.vmax, vmax) / 3.6;

direction = 12;
[runningtimes.IC12.regular , clearingtimes.IC12.regular ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC12.regular = clearingtimes.IC12.regular + settings.TT.blocktimes.releaseT; 

direction = 13;
[runningtimes.IC13.regular , clearingtimes.IC13.regular ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC13.regular = clearingtimes.IC13.regular + settings.TT.blocktimes.releaseT; 

% runningtimes.IC1.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
% runningtimes.IC1.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.IC1.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

direction = 2;
[runningtimes.IC02.regular , clearingtimes.IC02.regular ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC02.regular = clearingtimes.IC02.regular + settings.TT.blocktimes.releaseT; 

direction = 3;
[runningtimes.IC03.regular , clearingtimes.IC03.regular ] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction);
clearingtimes.IC03.regular = clearingtimes.IC03.regular + settings.TT.blocktimes.releaseT; 

% runningtimes.IC0.regular = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
% runningtimes.IC0.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
% runningtimes.IC0.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);


end

%% Regular
function [RT, CT] = determineRunningAndClearingTime(blocksections, blockSpeed, acc, decc, length, buffer_factor, direction)

if direction == 2 | direction == 3
    blockSpeed(1:end) = blockSpeed(end:-1:1);
    blocksections = sortrows(blocksections,'id','descend');
end

% Generate the basic one, running at full speeds all the time.
Nblocks = size(blocksections,1);
RT = zeros(Nblocks,1);
CT = zeros(Nblocks,1);

for bb = 1:Nblocks-1
    if blockSpeed(bb) == blockSpeed(bb+1) || blockSpeed(bb) < blockSpeed(bb+1)
        % We can continue at the same speed.
        try
            if blockSpeed(bb) > blockSpeed(bb-1)
                % We could accelerate from the previous one!
                ta = (blockSpeed(bb) - blockSpeed(bb-1)) / acc;
                xa = (blockSpeed(bb) - blockSpeed(bb-1))^2 / 2 / acc;
                RT(bb) = (blocksections.length(bb) - xa) / blockSpeed(bb) + ta;
                
                if xa > length
                    CT(bb) = (-blockSpeed(bb-1) + sqrt(blockSpeed(bb-1)^2 + 2*acc*length))/acc;
                else
                    CT(bb) = (length - xa) / blockSpeed(bb-1) + ta;
                end
                
            else
                RT(bb) = blocksections.length(bb) / blockSpeed(bb);
                CT(bb) = length / blockSpeed(bb);
            end
        catch
            RT(bb) = blocksections.length(bb) / blockSpeed(bb);
            CT(bb) = length / blockSpeed(bb);
        end
    elseif blockSpeed(bb) > blockSpeed(bb+1)
        % We have to deccelerate to the indicated speed from the start of
        % this section on!
        td = (blockSpeed(bb) - blockSpeed(bb+1)) / -decc;
        xd = (blockSpeed(bb) - blockSpeed(bb+1))^2 / 2 / -decc;
        RT(bb) = (blocksections.length(bb) - xd) / blockSpeed(bb+1) + td;
        
        if xd > length
            CT(bb) = (-blockSpeed(bb) + sqrt(blockSpeed(bb)^2 + 2*decc*length))/decc;
        else
            CT(bb) = (length - xd) / blockSpeed(bb) + td;
        end
        
    else
        try
            if blockSpeed(bb) > blockSpeed(bb-1)
                % We could accelerate!
                ta = (blockSpeed(bb) - blockSpeed(bb-1)) / acc;
                xa = (blockSpeed(bb) - blockSpeed(bb-1))^2 / 2 / acc;
                RT(bb) = (blocksections.length(bb) - xa) / blockSpeed(bb) + ta;
                
                if xa > length
                    CT(bb) = (-blockSpeed(bb-1) + sqrt(blockSpeed(bb-1)^2 + 2*acc*length))/acc;
                else
                    CT(bb) = (length - xa) / blockSpeed(bb-1) + ta;
                end
                
                
            elseif blockSpeed(bb) == blockSpeed(bb-1)
                % Run at the max speed.
                RT(bb) = blocksections.length(bb) / blockSpeed(bb);
                CT(bb) = length / blockSpeed(bb);
            end
        end
    end
end




% BUT!: if D contains only one block section, it's max. speed will depend
% on the characteristics at the switches, i.e. the blocks before and after
% it.
indexD = find(strcmp(blocksections.type,'D'));
nrBlocksD = size(indexD,1);
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
    CT(end) = length / blockSpeed(end);
elseif blockSpeed(end) > blockSpeed(end-1)
    ta = (blockSpeed(end) - blockSpeed(end-1)) / acc;
    xa = (blockSpeed(end) - blockSpeed(end-1))^2 / 2 / acc;
    RT(end) = (blocksections.length(end) - xa) / blockSpeed(end) + ta;
    
    if xa > length
        CT(bb) = (-blockSpeed(bb-1) + sqrt(blockSpeed(bb-1)^2 + 2*acc*length))/acc;
    else
        CT(bb) = (length - xa) / blockSpeed(bb-1) + ta;
    end
end

if direction == 2 | direction == 3
    RT(1:end) = RT(end:-1:1);
end

% In case of diverging tracks

if direction == 12 | direction == 13
    
   S2_index = find(strcmp(blocksections.type,'S2'));
    
   CT(S2_index) = 0;
   CT(S2_index-1) = 0;
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
