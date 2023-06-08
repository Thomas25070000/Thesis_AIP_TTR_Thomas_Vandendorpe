function original_timetable = generateOriginalTimetable(settings)

% Originally, they can just cross each other without any problems
% Assume trains run enter the section at full speed and also exit it like
% that.

% Additionally, generate the 'block table'
blocks = [];
counter = 0;
% Start with the ones of A.
aa = 0;
nrA = settings.infrastructure.blocks.A;
lengthA = settings.infrastructure.length.A;
while (aa < settings.infrastructure.blocks.A)
    aa = aa+1;
    counter = counter + 1;
    blocks(counter).id = counter;
    blocks(counter).length = lengthA/nrA;
    try
        blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
    catch
        % It is the first block!
        blocks(counter).distance = blocks(counter).length;
    end
    blocks(counter).type = 'A';
end

% Switch S1
counter = counter + 1;
blocks(counter).id = counter;
blocks(counter).length = settings.infrastructure.length.S1;
try
    blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
catch
    % It is the first block!
    blocks(counter).distance = blocks(counter).length;
end
blocks(counter).type = 'S1';
blocks(counter).direction = settings.infrastructure.switch.dirS1;

% Part D
dd = 0;
nrD = settings.infrastructure.blocks.D;
lengthD = settings.infrastructure.length.D;
while (dd < settings.infrastructure.blocks.D)
    dd = dd+1;
    counter = counter + 1;
    blocks(counter).id = counter;
    blocks(counter).length = lengthD/nrD;
    try
        blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
    catch
        % It is the first block!
        blocks(counter).distance = blocks(counter).length;
    end
    blocks(counter).type = 'D';
end

% Switch S2
counter = counter + 1;
blocks(counter).id = counter;
blocks(counter).length = settings.infrastructure.length.S2;
try
    blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
catch
    % It is the first block!
    blocks(counter).distance = blocks(counter).length;
end
blocks(counter).type = 'S2';
blocks(counter).direction = settings.infrastructure.switch.dirS2;

% Part B
bb = 0;
nrB = settings.infrastructure.blocks.B;
lengthB = settings.infrastructure.length.B;
while (bb < settings.infrastructure.blocks.B)
    bb = bb+1;
    counter = counter + 1;
    blocks(counter).id = counter;
    blocks(counter).length = lengthB/nrB;
    try
        blocks(counter).distance = blocks(counter-1).distance + blocks(counter).length;
    catch
        % It is the first block!
        blocks(counter).distance = blocks(counter).length;
    end
    blocks(counter).type = 'B';
end




end