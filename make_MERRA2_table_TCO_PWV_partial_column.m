function MERRA2_column = make_MERRA2_table_TCO_PWV_partial_column(year)
% this function is to make MERRA2 table data, written by Xiaoyi 15.June,2016

plot_columns = 1; % plot figure = 1; 
size_fig = 1/2; save_fig = 1;
Dyn_trop = 0; % 1= dynamic tropopause (PV=3.5), 0 = WMO thermal trop

year = num2str(year);
%data_dir = ['H:\work\MERRA\MERRA2_from_Sophie\' year];
data_dir = ['E:\H\work\MERRA\MERRA2_from_Sophie\' year];
working_dir = ['E:\H\work\MERRA\MERRA2_from_Sophie\' year '\partial_columns\table\'];
mkdir(working_dir);cd(data_dir);


start_time = [year '-01-01 00'];
end_time = [year '-12-31 21'];
if strcmp(year, '2017') % this change is due to we only have 2017 data until end of Oct 2017
    end_time = [year '-10-31 21'];
end
time_start = datenum(start_time,'yyyy-mm-dd HH');
time_end = datenum(end_time,'yyyy-mm-dd HH');
delta_time = 3/24;
time = time_start:delta_time:time_end;


%% read in profiles data
T = read_MERRA2_profiles('Temperature',year); % [K]
H = read_MERRA2_profiles('GPH',year);% [km]
p = read_MERRA2_profiles('Pressure',year); % [hPa]
O3 = read_MERRA2_profiles('O3',year);% this is mass mixing ratio! [kg/kg]
Qv = read_MERRA2_profiles('Qv',year);% specific humidity, this is mass mixing ratio! [kg/kg]

%% read in tropopause data
if Dyn_trop == 1
    combined_tropopause = read_MERRA2_tropopause('Dyn_Tropopauses',year);
else
    combined_tropopause = read_MERRA2_tropopause('WMO_Tropopauses',year);
end
combined_tropopause = combined_tropopause(1,:);

%% height grids
H_grid = 0.1:0.1:76; % Height from 0.1 km to 76 km on 0.1 km grids
N = size(T);

cd(working_dir);

%% main loop
for i =1:1:N(2)
    T_interp(:,i) = interp1(H(:,i),T(:,i),H_grid);
    O3_interp(:,i) = interp1(H(:,i),O3(:,i),H_grid).*(28.8/48);% this convert MMR to VMR
    
    Qv_interp(:,i) = interp1(H(:,i),Qv(:,i),H_grid);% this is MMR of H2O
    %%% int. the specific humidity to get the pwv
    p_interp(:,i) = interp1(H(:,i),p(:,i),H_grid);
    pwv_interp(:,i) = 28.8./8.314./1000.*Qv_interp(:,i).*p_interp(:,i)*100.*100./(T_interp(:,i));% calculate H2O pwv value
    TF = isnan(pwv_interp(:,i));
    pwv_interp_4int = pwv_interp(:,i);
    pwv_interp_4int(TF,:) = [];
    pwv_int(i) = sum(pwv_interp_4int);
    
    %%% int. the ozone VMR to get the TCO
    O3_DU_interp(:,i) = O3_interp(:,i).*p_interp(:,i).*100.*100.*6.02e23./(8.314.*T_interp(:,i));% Ozone (modlec/m^2) = VMR*pressure*Aav*d_h/(R*T)
    TF = isnan(O3_DU_interp(:,i));% calculate O3 mm value
    O3_DU_interp_4int = O3_DU_interp(:,i);
    O3_DU_interp_4int(TF,:) = [];
    O3_DU_int(i) = sum(O3_DU_interp_4int)/2.69e20;% calculate O3 in DU value, 1DU = 2.69e20 molec/m^2
    
    %% int. the ozone VMR to get trop column
    TF_trop = combined_tropopause(i) < H_grid;
    O3_DU_interp_4int = O3_DU_interp(:,i);
    O3_DU_interp_4int(TF_trop,:) = [];
    TF = isnan(O3_DU_interp_4int(:,1));
    O3_DU_interp_4int(TF,:) = [];
    O3_DU_int_trop(i) = sum(O3_DU_interp_4int)/2.69e20;
    
    %% int. the ozone VMR to get strat column
    O3_DU_interp_4int = O3_DU_interp(:,i);
    O3_DU_interp_4int(~TF_trop,:) = [];
    TF = isnan(O3_DU_interp_4int(:,1));
    O3_DU_interp_4int(TF,:) = [];
    O3_DU_int_strat(i) = sum(O3_DU_interp_4int)/2.69e20;
    
    %% int. the ozone VMR to get trop column (0-2 km)
    TF_trop = 2 < H_grid;
    O3_DU_interp_4int = O3_DU_interp(:,i);
    O3_DU_interp_4int(TF_trop,:) = [];
    TF = isnan(O3_DU_interp_4int(:,1));
    O3_DU_interp_4int(TF,:) = [];
    O3_DU_int_trop_2km(i) = sum(O3_DU_interp_4int)/2.69e20;
    
    %% int. the ozone VMR to get trop column (0-3 km)
    TF_trop = 3 < H_grid;
    O3_DU_interp_4int = O3_DU_interp(:,i);
    O3_DU_interp_4int(TF_trop,:) = [];
    TF = isnan(O3_DU_interp_4int(:,1));
    O3_DU_interp_4int(TF,:) = [];
    O3_DU_int_trop_3km(i) = sum(O3_DU_interp_4int)/2.69e20;
end

if plot_columns ~= 0
    %%%%% H2O total pwv time-series %%%%%%%
    figure;
    plot(time,pwv_int);
    datetick('x','mmm-dd','keeplimits');
    ylabel('H_2O pwv [mm]');
    xlabel(year);
    title('MERRA2 H_2O pwv');
    print_setting(size_fig,save_fig,['MERRA2_H2O_pwv_' year]);
    
    %%%%% O3 total DU time-series %%%%%%%
    figure;hold all;
    plot(time,O3_DU_int);
    plot(time,O3_DU_int_trop);
    plot(time,O3_DU_int_strat);
    plot(time,O3_DU_int_trop_2km);
    plot(time,O3_DU_int_trop_3km);
    datetick('x','mmm-dd','keeplimits');
    ylabel('O_3 [DU]');
    xlabel(year);
    legend('TCO','trop','strat','trop 0-2 km','trop 0-3 km');
    if Dyn_trop == 1
        title('MERRA2 O_3 DU (dynamic tropopause)');
        print_setting(size_fig,save_fig,['MERRA2_O3_DU_dyn_trop_and_surface_o3' year]);
    else
        title('MERRA2 O_3 DU (thermal tropopause)');
        print_setting(size_fig,save_fig,['MERRA2_O3_DU_thermal_trop_' year]);
    end
    

end

%%%%% save table data %%%%%%%%%%
MERRA2_column = table;
MERRA2_column.UTC = time';
MERRA2_column.PWV = pwv_int';
MERRA2_column.Ozone = O3_DU_int';
MERRA2_column.Ozone_trop = O3_DU_int_trop';
MERRA2_column.Ozone_strat = O3_DU_int_strat';
MERRA2_column.Ozone_trop_2km = O3_DU_int_trop_2km';
MERRA2_column.Ozone_trop_3km = O3_DU_int_trop_3km';
MERRA2_column.Properties.VariableUnits = {'','mm','DU','DU','DU','DU','DU'};
save('MERRA2_column.mat','MERRA2_column');


function combined_profiles = read_MERRA2_profiles(var,year)
%data_dir = ['H:\work\MERRA\MERRA2_from_Sophie\' year];
data_dir = ['E:\H\work\MERRA\MERRA2_from_Sophie\' year];
cd(data_dir);
%files = ls(['GEOS5MERRA2_DynVarsAtEureka_' year '_*.nc4']);
files = ls(['GEOS5MERRA2_DynVarsAtEureka_' year '*.nc4']);
N = size(files);

for i =1:1:N(1)
    profile{i} = ncread(files(i,:),var);
end
NN = size(profile{1,1});
combined_size = NN(2)*8;
for i =1:1:N(1)
    end_point = combined_size -8 + i;
    combined_profiles(:,i:8:end_point) = profile{i};
end


function combined_tropopause = read_MERRA2_tropopause(var,year)
% var = 'WMO_Tropopauses'
% var = 'Dyn_Tropopauses'
%data_dir = ['H:\work\MERRA\MERRA2_from_Sophie\' year];
data_dir = ['E:\H\work\MERRA\MERRA2_from_Sophie\' year];
cd(data_dir);
files = ls(['GEOS5MERRA2_TropopausesAtEureka_' year '*.nc4']);
N = size(files);

for i =1:1:N(1)
    trop_height{i} = ncread(files(i,:),var);
    if strcmp(var,'Dyn_Tropopauses')
        trop_height{i} = squeeze(trop_height{i}(3,:,:)); % only use PV = 3.5 as Dyn trop
    end
end

NN = size(trop_height{1,1});
combined_size = NN(2)*8;
for i =1:1:N(1)
    end_point = combined_size -8 + i;
    combined_tropopause(:,i:8:end_point) = trop_height{i};
end

