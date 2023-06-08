function runningtimes = generateRunningTimes_case3(blocksections, settings)

vyellow = settings.infrastructure.vyellow / 3.6;

buffer_factor = 1 + settings.TT.buffer;

Ntracks = size(settings.tracks,1);
Nblocks = max(blocksections.id);

for tt = 1:Ntracks
    trbs = find(blocksections.track == settings.tracks.nr(tt));
    
    if ~settings.tracks.closed(tt)
        %% Start with the disruption for L trains.
        vmax = settings.trains.speed.R;
        acc = settings.trains.acc.R;
        decc = -settings.trains.decc.R;

        direction = 12;
        blockSpeed = min(blocksections.vdis(trbs), vmax) / 3.6;     % Include conversion to m/s
        runningtimes.L12.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
    %     runningtimes.L1.disrupted_acc(tt,:) = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
    %     runningtimes.L1.disrupted_stop(tt,:) = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
        direction = 13;
        blockSpeed = min(blocksections.vdis(trbs), vmax) / 3.6;     % Include conversion to m/s
        runningtimes.L13.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
    %     runningtimes.L1.disrupted_acc(tt,:) = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
    %     runningtimes.L1.disrupted_stop(tt,:) = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

        direction = 02;
        runningtimes.L02.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
    %     runningtimes.L0.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
    %     runningtimes.L0.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
        direction = 03;
        runningtimes.L03.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
    %     runningtimes.L0.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
    %     runningtimes.L0.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
    else
        runningtimes.L12.disrupted(tt,:) = zeros(1,Nblocks);
        runningtimes.L02.disrupted(tt,:) = zeros(1,Nblocks);
        runningtimes.L13.disrupted(tt,:) = zeros(1,Nblocks);
        runningtimes.L03.disrupted(tt,:) = zeros(1,Nblocks);
    end


    %% L in normal case
    vmax = settings.trains.speed.R;
    acc = settings.trains.acc.R;
    decc = -settings.trains.decc.R;
    blockSpeed = min(blocksections.vreg(trbs), vmax) / 3.6;

    direction = 12;
    runningtimes.L12.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.L1.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.L1.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
    direction = 13;
    runningtimes.L13.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.L1.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.L1.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

    direction = 02;
    runningtimes.L02.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.L0.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.L0.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
    direction = 03;
    runningtimes.L03.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.L0.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.L0.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);



    %% IC in disruption
    vmax = settings.trains.speed.IC;
    acc = settings.trains.acc.IC;
    decc = -settings.trains.decc.IC;

    if ~settings.tracks.closed(tt)
        direction = 12;
        blockSpeed = min(blocksections.vdis(trbs), vmax) / 3.6;     % Include conversion to m/s
        runningtimes.IC12.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%         runningtimes.IC1.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%         runningtimes.IC1.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
        direction = 13;
        blockSpeed = min(blocksections.vdis(trbs), vmax) / 3.6;     % Include conversion to m/s
        runningtimes.IC13.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%         runningtimes.IC1.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%         runningtimes.IC1.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

        direction = 02;
        runningtimes.IC02.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%         runningtimes.IC0.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%         runningtimes.IC0.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
        direction = 03;
        runningtimes.IC03.disrupted(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%         runningtimes.IC0.disrupted_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%         runningtimes.IC0.disrupted_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
    else
        runningtimes.IC12.disrupted(tt,:) = zeros(1,Nblocks);
        runningtimes.IC02.disrupted(tt,:) = zeros(1,Nblocks);
        runningtimes.IC13.disrupted(tt,:) = zeros(1,Nblocks);
        runningtimes.IC03.disrupted(tt,:) = zeros(1,Nblocks);
    end

    %% IC in normal case
    vmax = settings.trains.speed.IC;
    acc = settings.trains.acc.IC;
    decc = -settings.trains.decc.IC;
    
    blockSpeed = min(blocksections.vreg(trbs), vmax) / 3.6;

    direction = 12;
    runningtimes.IC12.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.IC1.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.IC1.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
    direction = 13;
    runningtimes.IC13.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.IC1.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.IC1.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);

    direction = 02;
    runningtimes.IC02.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.IC0.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.IC0.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
   
    direction = 03;
    runningtimes.IC03.regular(tt,:) = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction);
%     runningtimes.IC0.regular_acc = determineRunningTimeAcc(blocksections, blockSpeed, acc, buffer_factor,direction);
%     runningtimes.IC0.regular_stop = determineRunningTimeBrake(blocksections, blockSpeed, decc, vyellow, buffer_factor, direction);
   

end

end

%% Regular
function RT = determineRunningTime(blocksections, blockSpeed, acc, decc, buffer_factor, direction)

if direction == 02 || direction == 03
    blockSpeed(1:end) = blockSpeed(end:-1:1);
    blocksections = sortrows(blocksections,'id','descend');
end

% Generate the basic one, running at full speeds all the time.
Nblocks = max(blocksections.id);
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

if direction == 02 || direction == 03
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

if direction == 02 || direction == 03
    blockSpeed(1:end) = blockSpeed(end:-1:1);
    blocksections = sortrows(blocksections,'id','descend');
end

% Generate the basic one, running at full speeds all the time.
Nblocks = max(blocksections.id);
RT = zeros(Nblocks,1);

% For the first block, assuming incoming at full speed.
RT(1) = blockSpeed(1) / abs(decc) + blocksections.length(1) / vyellow;
RT(2:end) = blockSpeed(1:end-1) / abs(decc) + blocksections.length(2:end) / vyellow;

if direction == 02 || direction == 03
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

if direction == 02 || direction == 03
    blockSpeed(1:end) = blockSpeed(end:-1:1);
    blocksections = sortrows(blocksections,'id','descend');
end

% Generate the basic one, running at full speeds all the time.
Nblocks = max(blocksections.id);
RT = zeros(Nblocks,1);

% For the first block, assuming incoming at full speed.
targetSpeed = min(blockSpeed(2:end), blockSpeed(1:end-1));
targetSpeed(Nblocks) = blockSpeed(end);
RT = targetSpeed / acc + blocksections.length ./ targetSpeed;

if direction == 02 || direction == 03
    RT(1:end) = RT(end:-1:1);
end

RT = round((RT * buffer_factor)');

if ~all(RT)
    disp('Problem');
    disp(RT);
    pause;
end

end
