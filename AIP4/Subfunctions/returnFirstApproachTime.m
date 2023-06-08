function [time, regtime] = returnFirstApproachTime(blocksections, settings, traintype, direction)

lengthB = max(blocksections.length);

if ~ischar(traintype)
    traintype = traintype{:};
end

switch traintype
    case {'IC','THA'}
        acc = settings.trains.acc.IC;
    case 'R'
        acc = settings.trains.acc.R;
end

if direction
    % Entering from the left.
    if blocksections.vmax_disruption_dir1(1) >= blocksections.vmax_disruption_dir1(2)
        speed = blocksections.vmax_disruption_dir1(2) / 3.6;
    else
        speed = blocksections.vmax_disruption_dir1(1) / 3.6;
    end
    if blocksections.vmax(1) >= blocksections.vmax(2)
        regspeed = blocksections.vmax(2) / 3.6;
    else
        regspeed = blocksections.vmax(1) / 3.6;
    end
else
    % From right
    if blocksections.vmax_disruption_dir0(end) >= blocksections.vmax_disruption_dir0(end-1)
        speed = blocksections.vmax_disruption_dir0(end-1) / 3.6;
    else
        speed = blocksections.vmax_disruption_dir0(end) / 3.6;
    end
    if blocksections.vmax(end) >= blocksections.vmax(end-1)
        regspeed = blocksections.vmax(end-1) / 3.6;
    else
        regspeed = blocksections.vmax(end) / 3.6;
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