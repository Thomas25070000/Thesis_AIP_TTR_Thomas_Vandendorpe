function [time, regtime] = returnFirstApproachTime_case3(blocksections, settings, traintype, direction, track)

blocksections(blocksections.track ~= track,:) = [];

lengthB = max(blocksections.length);

if ~ischar(traintype)
    traintype = traintype{:};
end

switch traintype
    case {'IC','THA'}
        acc = settings.trains.acc.IC;
    case 'L'
        acc = settings.trains.acc.R;
end

if direction == 12 || direction == 13
    % Entering from the left.
    if blocksections.vdis(1) >= blocksections.vdis(2)
        speed = blocksections.vdis(2) / 3.6;
    else
        speed = blocksections.vdis(1) / 3.6;
    end
    if blocksections.vreg(1) >= blocksections.vreg(2)
        regspeed = blocksections.vreg(2) / 3.6;
    else
        regspeed = blocksections.vreg(1) / 3.6;
    end
else
    % From right
    if blocksections.vdis(end) >= blocksections.vdis(end-1)
        speed = blocksections.vdis(end-1) / 3.6;
    else
        speed = blocksections.vdis(end) / 3.6;
    end
    if blocksections.vreg(end) >= blocksections.vreg(end-1)
        regspeed = blocksections.vreg(end-1) / 3.6;
    else
        regspeed = blocksections.vreg(end) / 3.6;
    end
end

% For the disrupted situation.
ta = speed / acc;
xa = speed^2 / 2 / acc;
time = (lengthB - xa) / speed + ta;
time = round(time * (1 + settings.TT.buffer));
    
% For the regular situation.
ta = regspeed / acc;
xa = regspeed^2 / 2 / acc;
regtime = (lengthB - xa) / regspeed + ta;
regtime = round(regtime * (1 + settings.TT.buffer));

if settings.TT.noFirstApproach
    time = 0;
    regtime = 0;
end


end