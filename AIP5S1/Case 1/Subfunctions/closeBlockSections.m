function blocksections = closeBlockSections(blocksections, settings)

% Identify which block sections will have to be shared.
blocksections.closed(1) = 0;
closed_dir = settings.disruption.direction;

first_S1 = 1;
first_S2 = 1;
for bb = 1:size(blocksections,1)
    switch blocksections.type{bb}
        case 'A'
            blocksections.closed(bb) = 0;
        case 'S1'
%             if first_S1
%                 % First part can always be used.
%                 blocksections.closed(bb) = 0;
%                 first_S1 = 0;
%             elseif blocksections.direction{bb} ~= closed_dir
                % This is the second part, and they cannot be used if there
                % is a mismatch between the CLOSED direction and their own
                % direction.
                blocksections.closed(bb) = 1;
%             else
%                 blocksections.closed(bb) = 0;
%             end
        case 'D'
            blocksections.closed(bb) = 1;
        case 'S2'
%             if ~first_S2
%                 % Last part can always be used.
%                 blocksections.closed(bb) = 0;
%             elseif blocksections.direction{bb} ~= closed_dir
                % This is the second part, and they cannot be used if there
                % is a mismatch between the CLOSED direction and their own
                % direction.
                blocksections.closed(bb) = 1;
%                 first_S2 = 0;
%             else
%                 blocksections.closed(bb) = 0;
%                 first_S2 = 0;
%             end
        case 'B'
            % Can just be used.
    end
end

closeddir = settings.disruption.direction;

blocksections.vmax_disruption_dir12(1) = 0;
blocksections.vmax_disruption_dir13(1) = 0;
blocksections.vmax_disruption_dir02(1) = 0;
blocksections.vmax_disruption_dir03(1) = 0;
for bb = 1:size(blocksections,1)
    if blocksections.closed(bb) == 1
        if ismember(blocksections.type(bb),{'S1'})
            if closeddir
                % If train goes in closed direction, it has to switch to
                % the other track over the switch
                blocksections.vmax_disruption_dir12(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS1);
                blocksections.vmax_disruption_dir13(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS1);
                blocksections.vmax_disruption_dir02(bb) = blocksections.vmax(bb);
                blocksections.vmax_disruption_dir03(bb) = blocksections.vmax(bb);
            else
                blocksections.vmax_disruption_dir12(bb) = blocksections.vmax(bb);
                blocksections.vmax_disruption_dir13(bb) = blocksections.vmax(bb);
                blocksections.vmax_disruption_dir02(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS1);
                blocksections.vmax_disruption_dir03(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS1);
            end
        elseif ismember(blocksections.type(bb),{'S2'})
            if closeddir
                % If train goes in closed direction, it has to switch to
                % the other track over the switch
                blocksections.vmax_disruption_dir12(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS2);
                blocksections.vmax_disruption_dir13(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS2);
                blocksections.vmax_disruption_dir02(bb) = blocksections.vmax(bb);
                blocksections.vmax_disruption_dir03(bb) = blocksections.vmax(bb);
            else
                blocksections.vmax_disruption_dir12(bb) = blocksections.vmax(bb);
                blocksections.vmax_disruption_dir13(bb) = blocksections.vmax(bb);
                blocksections.vmax_disruption_dir02(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS2);
                blocksections.vmax_disruption_dir03(bb) = min(settings.disruption.maxspeed, settings.infrastructure.switch.maxvS2);
            end
        else
            blocksections.vmax_disruption_dir12(bb) = settings.disruption.maxspeed;
            blocksections.vmax_disruption_dir13(bb) = settings.disruption.maxspeed;
            blocksections.vmax_disruption_dir02(bb) = settings.disruption.maxspeed;
            blocksections.vmax_disruption_dir03(bb) = settings.disruption.maxspeed;
        end
    else
        blocksections.vmax_disruption_dir12(bb) = blocksections.vmax(bb);
        blocksections.vmax_disruption_dir13(bb) = blocksections.vmax(bb);
        blocksections.vmax_disruption_dir02(bb) = blocksections.vmax(bb);
        blocksections.vmax_disruption_dir03(bb) = blocksections.vmax(bb);
    end
end

end