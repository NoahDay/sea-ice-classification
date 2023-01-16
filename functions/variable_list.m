function list_out = variable_list(list_in)
% list_in can be one of the following:
    % "static": variables that are not time derivatives
    % "dynamics": dynamic and thermodynamic variables
    % "thermodynamics": thermodynamic terms only
    if list_in == "static"
        list_out = {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'};
    elseif list_in == "dynamics"
        list_out = {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg', 'uvel','vvel','strength','divu','shear','daidtd','daidtt','dagedtd','dagedtt','dafsd_latg','dafsd_latm','dafsd_newi','dafsd_weld','dafsd_wave'};
    end
end

% Static: {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg'}
% Dynamics: {'aice','hi','hs','fsdrad','sice','iage','vlvl','vrdg', 'uvel','vvel','strength','divu','shear','daidtd','daidtt','dagedtd','dagedtt','dafsd_latg','dafsd_latm','dafsd_newi','dafsd_weld','dafsd_wave'};
