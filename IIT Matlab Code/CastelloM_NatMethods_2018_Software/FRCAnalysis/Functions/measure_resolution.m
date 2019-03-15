function [cutOff, resolution] = measureResolution(Scale, cutOff, params)
fields = fieldnames(cutOff);

for i=1:numel(fields)
    if cutOff.(fields{i}) == -1
        if ( ...
                params.theta ~= 0 && ~isempty(strfind(fields{i}, '_smallAngles')) ...
                || ...
                params.theta ~= pi/2 && ~isempty(strfind(fields{i}, '_largeAngles')) ...
                )
            warning(strcat('FRC curve is not going under threshold: ', fields{i}));
        end
        cutOff.(fields{i}) = NaN;
        resolution.(fields{i}) = NaN;
    else
        resolution.(fields{i}) = 1000/Scale(cutOff.(fields{i}));
    end
end
end
