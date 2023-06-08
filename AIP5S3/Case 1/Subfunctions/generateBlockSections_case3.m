function blocksections = generateBlockSections_case3(settings)

% Originally, they can just cross each other without any problems
% Assume trains run enter the section at full speed and also exit it like
% that.
Ntracks = size(settings.tracks,1);
blocks = [];
counter = 0;

for tt = 1:Ntracks
    % Additionally, generate the 'block table'
    counter_track = 0;
    % Start with the ones of A.
    aa = 0;
    nrA = settings.infrastructure.blocks.A;
    lengthA = settings.infrastructure.length.A;
    while (aa < settings.infrastructure.blocks.A)
        aa = aa+1;
        counter = counter + 1;
        counter_track = counter_track + 1;
        blocks(counter).track = settings.tracks.nr(tt);
        blocks(counter).id = counter_track;
        blocks(counter).ID = settings.tracks.nr(tt) * 10 + counter_track;
        blocks(counter).length = lengthA/nrA;
        try
            blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
        catch
            % It is the first block!
            blocks(counter).distance = blocks(counter).length;
        end
        blocks(counter).type = 'A';
        blocks(counter).vreg = settings.tracks.vreg(tt);
        blocks(counter).vdis = settings.tracks.vdis(tt);
        blocks(counter).closed = settings.tracks.closed(tt);
    end

%     %Switch S1
%     %Note that a switch consists of two block sections!
%     for ii=1:1
%         counter = counter + 1;
%         blocks(counter).id = counter;
%         blocks(counter).length = settings.infrastructure.length.S1;
%         try
%             blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
%         catch
%             % It is the first block!
%             blocks(counter).distance = blocks(counter).length;
%         end
%         blocks(counter).type = 'S1';
%         %blocks(counter).direction = settings.infrastructure.switch.dirS1;
%         blocks(counter).direction = 1;
%         blocks(counter).vmax = settings.infrastructure.vmax;
%     end

    % Part D
    dd = 0;
    nrD = settings.infrastructure.blocks.D;
    lengthD = settings.infrastructure.length.D;
    while (dd < settings.infrastructure.blocks.D)
        dd = dd+1;
        counter = counter + 1;
        counter_track = counter_track + 1;
        blocks(counter).track = settings.tracks.nr(tt);
        blocks(counter).id = counter_track;
        blocks(counter).ID = settings.tracks.nr(tt) * 10 + counter_track;
        blocks(counter).length = lengthD/nrD;
        try
            if blocks(counter).track == blocks(counter-1).track
                blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
            else
                blocks(counter).distance = blocks(counter).length;
            end
        catch
            % It is the first block!
            blocks(counter).distance = blocks(counter).length;
        end
        blocks(counter).type = 'D';
        blocks(counter).vreg = settings.tracks.vreg(tt);
        blocks(counter).vdis = settings.tracks.vdis(tt);
        blocks(counter).closed = settings.tracks.closed(tt);
    end
    

%     % Switch S2
%     % Note that a switch consists of two block sections!
%     for ii=1:1
%         counter = counter + 1;
%         blocks(counter).id = counter;
%         blocks(counter).length = settings.infrastructure.length.S2;
%         try
%             blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
%         catch
%             % It is the first block!
%             blocks(counter).distance = blocks(counter).length;
%         end
%         blocks(counter).type = 'S2';
%         %blocks(counter).direction = settings.infrastructure.switch.dirS2;
%         blocks(counter).direction = 1;
%         blocks(counter).vmax = settings.infrastructure.vmax;
%     end

    % Part B
    bb = 0;
    counter_before_B = counter;
    nrB = settings.infrastructure.blocks.B;
    lengthB = settings.infrastructure.length.B;
    while (bb < settings.infrastructure.blocks.B)
        bb = bb+1;
        counter = counter + 1;
        counter_track = counter_track + 1;
        blocks(counter).track = settings.tracks.nr(tt);
        blocks(counter).id = counter_track;
        blocks(counter).ID = settings.tracks.nr(tt) * 10 + counter_track;
        blocks(counter).length = lengthB/nrB;
        try
            blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
        catch
            % It is the first block!
            blocks(counter).distance = blocks(counter).length;
        end
        blocks(counter).type = 'B';
        blocks(counter).vreg = settings.tracks.vreg(tt);
        blocks(counter).vdis = settings.tracks.vdis(tt);
        blocks(counter).closed = settings.tracks.closed(tt);
    end

    % Part C
    bb = 0;
    nrC = settings.infrastructure.blocks.C;
    lengthC = settings.infrastructure.length.C;
    while (bb < settings.infrastructure.blocks.C)
        bb = bb+1;
        counter = counter + 1;
        counter_track = counter_track + 1;
        blocks(counter).track = settings.tracks.nr(tt);
        blocks(counter).id = counter_track;
        blocks(counter).ID = settings.tracks.nr(tt) * 10 + counter_track;
        blocks(counter).length = lengthC/nrC;
        if counter == 8 || counter == 19
            blocks(counter).distance = blocks(3).distance+blocks(counter).length;
        else
            blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
        end
        blocks(counter).type = 'C';
        blocks(counter).vreg = settings.tracks.vreg(tt);
        blocks(counter).vdis = settings.tracks.vdis(tt);
        blocks(counter).closed = settings.tracks.closed(tt);
    end
end
for alpha = 1:length(blocks)
    blocks(alpha).distance = blocks(alpha).distance - 400;
end
blocks(1).length = blocks(1).length - 400;
blocks(12).length = blocks(1).length - 400;
blocksections = struct2table(blocks);


end