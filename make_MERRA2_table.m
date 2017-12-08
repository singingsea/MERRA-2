function T = make_MERRA2_table(year)

% grab Dynamic tropopause or WMO thermal one 
Dynamic = 0; % 1 = dynamic, 0 = WMO thermal

if Dynamic == 1
    %var_names = {'Dyn_Tropopauses', 'Temp_at_Dyn_Tropopause', 'Press_at_Dyn_Tropopause', 'Theta_at_Dyn_Tropopause', 'DynSurf'};
    var_names = {'Dyn_Tropopauses'};
    N = size(var_names);
    for i = 1:1:N(2)
        Output = plot_MERRA2_multiyears_v2(year,var_names{1,i})
        data(:,i) = Output(:,2);
        PV = Output(:,3);
        UTC = Output(:,1);
    end
    
    T = array2table(data, 'VariableNames',var_names);
    T.PV = PV;
    T.UTC = UTC;
    T = [T(:,3),T(:,1:2)];
    T = sortrows(T,1);
    
    T = sortrows(T,3);
    TF = T.PV == 1.5;
    T2 = table(T.UTC(TF),'VariableNames',{'UTC'});
    %T2 = table(T.Dyn_Tropopauses(TF),'VariableNames',{'Dyn_Tropopauses_PV1_5'});
    T2.Dyn_Tropopauses_PV1_5 = T.Dyn_Tropopauses(TF);
    TF = T.PV == 2.0;
    T2.Dyn_Tropopauses_PV2 = T.Dyn_Tropopauses(TF);
    TF = T.PV == 3.5;
    T2.Dyn_Tropopauses_PV3_5 = T.Dyn_Tropopauses(TF);
    TF = T.PV == 4.5;
    T2.Dyn_Tropopauses_PV4_5 = T.Dyn_Tropopauses(TF);
    TF = T.PV == 6;
    T2.Dyn_Tropopauses_PV6 = T.Dyn_Tropopauses(TF);
    
    T = T2;
else
    var_names = {'WMO_Tropopauses', 'Temp_at_WMO_Tropopause','Press_at_WMO_Tropopause','Theta_at_WMO_Tropopause','WMOSurf'}
    N = size(var_names);
    
    for i = 1:1:N(2)
        Output = plot_MERRA2_WMO_v2(year,var_names{1,i});
        data(:,i) = Output(:,2);
        UTC = Output(:,1);
    end
    
    T = array2table(data, 'VariableNames',var_names);
    T.UTC = UTC;
    T = [T(:,6),T(:,1:5)];
    T = sortrows(T,1);
    
end