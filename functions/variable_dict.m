function var_out = variable_dict(var_in)
    for i = 1:numel(var_in)
        if var_in{i} == "aice"
            var_out{i} = 'SIC';
        elseif var_in{i} == "hi"
            var_out{i} = 'Ice thickness';
        elseif var_in{i} == "iage"
            var_out{i} = 'Ice age';
        elseif var_in{i} == "fsdrad"
            var_out{i} = 'Floe mean radius';
        elseif var_in{i} == "pancake"
            var_out{i} = 'Pancake ice';
        elseif var_in{i} == "large"
            var_out{i} = 'Sheet ice';
        elseif var_in{i} == "wave_sig_ht"
            var_out{i} = 'Sig. wave height';
        elseif var_in{i} == "frazil"
            var_out{i} = 'Frazil';
        elseif var_in{i} == "frzmlt"
            var_out{i} = 'Freeze/melt pot.';
        elseif var_in{i} == "congel"
            var_out{i} = 'Congelation rate';
        elseif var_in{i} == "ice_present"
            var_out{i} = 'Ice present';
        elseif var_in{i} == "peak_period"
            var_out{i} = 'Peak period';
        elseif var_in{i} == "vel"
            var_out{i} = 'Ice speed';
        elseif var_in{i} == "strint"
            var_out{i} = 'Internal stress';
        elseif var_in{i} == "Tair"
            var_out{i} = 'Air temp.';
        elseif var_in{i} == "dafsd_latg"
            var_out{i} = 'Lat growth';
        elseif var_in{i} == "dafsd_latm"
            var_out{i} = 'Lat melt';
        elseif var_in{i} == "dafsd_newi"
            var_out{i} = 'New ice';
        elseif var_in{i} == "dafsd_weld"
            var_out{i} = 'Welding';
        elseif var_in{i} == "dafsd_wave"
            var_out{i} = 'Wave breakup';
        elseif var_in{i} == "dafsd_wave"
            var_out{i} = 'Wave breakup';
        elseif var_in{i} == "sice"
            var_out{i} = 'Ice salinity';
        elseif var_in{i} == "alvl"
            var_out{i} = 'Level ice area frac.';
        elseif var_in{i} == "vlvl"
            var_out{i} = 'Level ice vol. frac.';
        elseif var_in{i} == "krdgn"
            var_out{i} = 'Mean ridge thickness per thickness of ridging ice';
        elseif var_in{i} == "ardg"
            var_out{i} = 'Fractional area of ridged ice';
        elseif var_in{i} == "vrdg"
            var_out{i} = 'Fractional volume of ridged ice';
        elseif var_in{i} == "hs"
            var_out{i} = 'Snow thickness';
        elseif var_in{i} == "FYarea"
            var_out{i} = 'First year ice area';
        elseif var_in{i} == "full_fsd"
            for j = 1:12
                var_out{i+j-1} = sprintf('FSD cat. %g',j);
            end
        elseif var_in{i} == "full_fstd"
            count = 1;
           for nc = 1:5
                for j = 1:12
                    var_out{i+count-1} = sprintf('FSTD cat. (%g %g)',j,nc);
                    count = count + 1;
                end
           end
        elseif var_in{i} == "uvel"
            var_out{i} = 'Ice u-speed';
        elseif var_in{i} == "vvel"
            var_out{i} = 'Ice v-speed';
        elseif var_in{i} == "strength"
            var_out{i} = 'Ice strength';
        elseif var_in{i} == "divu"
            var_out{i} = 'Divergence';
        elseif var_in{i} == "shear"
            var_out{i} = 'Shear';
        elseif var_in{i} == "daidtt"
            var_out{i} = 'Change in SIC (thermo.)';
        elseif var_in{i} == "daidtd"
            var_out{i} = 'Change in SIC (dynamics)';
        elseif var_in{i} == "dagedtt"
            var_out{i} = 'Change in age (thermo.)';
        elseif var_in{i} == "dagedtd"
            var_out{i} = 'Change in age (dynamics)';
        else
            var_out{i} = var_in{i};
        end
    end
end