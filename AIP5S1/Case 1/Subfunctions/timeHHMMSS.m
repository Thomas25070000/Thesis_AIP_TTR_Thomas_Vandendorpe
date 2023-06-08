function text = timeHHMMSS(time)
text = [];
    if time >= 3600
        if mod(time,60) >= 10
            if mod(time,3600) >= 10
                text = [int2str(floor(time/3600)) ':' int2str(floor(mod(time,3600)/60))...
                    ':' int2str(mod(time,60))];
            else
                text = [int2str(floor(time/3600)) ':' int2str(floor(mod(time,3600)/60))...
                    ':0' int2str(mod(time,60))];
            end
        else
            if mod(time,3600) >= 10
                text = [int2str(floor(time/3600)) ':' int2str(floor(mod(time,3600)/60))...
                    ':0' int2str(mod(time,60))];
            else
                text = [int2str(floor(time/3600)) ':0' int2str(floor(mod(time,3600)/60))...
                    ':0' int2str(mod(time,60))];
            end
        end 
    elseif time >= 60 && time < 600
        if mod(time,60) >= 10
            text = ['0:0' int2str(floor(mod(time,3600)/60)) ':' int2str(mod(time,60))];
        else
            text = ['0:0' int2str(floor(mod(time,3600)/60)) ':0' int2str(mod(time,60))];
        end
    elseif time >= 60
        if mod(time,60) >= 10
            text = ['0:' int2str(floor(mod(time,3600)/60)) ':' int2str(mod(time,60))];
        else
            text = ['0:' int2str(floor(mod(time,3600)/60)) ':0' int2str(mod(time,60))];
        end
    elseif mod(time,60) >= 10
        text = ['0:00:' int2str(mod(time,60))];
    elseif mod(time,60) >= 0
        text = ['0:00:0' int2str(mod(time,60))];
    end
end