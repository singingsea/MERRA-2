function plot_MERRA_2D_time_series(year)
% this function is to plot MERRA2 data, written by Xiaoyi 30.May,2016

plot_2D_profiles_time_series = 1; % plot 2D profiles? 1 = yes, 0 = no;
plot_ozone_partial_column = 1; % plot partial column? 1 = yes, 0 = no
year = num2str(year);
%load('H:\work\MERRA\MERRA2_from_Sophie\example_profile.mat');
data_dir = ['E:\H\work\MERRA\MERRA2_from_Sophie\' year];
working_dir = ['E:\H\work\MERRA\MERRA2_from_Sophie\' year '\columns\'];
mkdir(working_dir);cd(data_dir);
save_fig = 1;size_fig = 1/2;
top_height = 50; % top of the GPH height on plots [km]


start_time = [year '-01-01 00'];
end_time = [year '-12-31 21'];
time_start = datenum(start_time,'yyyy-mm-dd HH');
time_end = datenum(end_time,'yyyy-mm-dd HH');
delta_time = 3/24;
time = time_start:delta_time:time_end;

%time = 1:1:2920;
T = read_MERRA2_profiles('Temperature',year); % [K]
H = read_MERRA2_profiles('GPH',year);% [km]
PV = read_MERRA2_profiles('PV',year);% unit: [K m^2 (kg s)^-1]
p = read_MERRA2_profiles('Pressure',year); % [hPa]
O3 = read_MERRA2_profiles('O3',year);% this is mass mixing ratio! [kg/kg]
Theta = read_MERRA2_profiles('Theta',year); % [K]
sPV = read_MERRA2_profiles('sPV',year); % [s^-1]
dT_dZ = read_MERRA2_profiles('dT_dZ',year); % [K km-1]
Qv = read_MERRA2_profiles('Qv',year);% specific humidity, this is mass mixing ratio! [kg/kg]

H_grid = 0.1:0.1:76; % Height from 0.1 km to 76 km on 0.1 km grids
N = size(T);

cd(working_dir);

for i =1:1:N(2)
    T_interp(:,i) = interp1(H(:,i),T(:,i),H_grid);
    PV_interp(:,i) = interp1(H(:,i),PV(:,i),H_grid);
    O3_interp(:,i) = interp1(H(:,i),O3(:,i),H_grid).*(28.8/48);% this convert MMR to VMR
    Theta_interp(:,i) = interp1(H(:,i),Theta(:,i),H_grid);
    sPV_interp(:,i) = interp1(H(:,i),sPV(:,i),H_grid);
    dT_dZ_interp(:,i) = interp1(H(:,i),dT_dZ(:,i),H_grid);
    
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
    if plot_ozone_partial_column == 1
        TF2 = H_grid(:) < 3; % filter out O3 below 3 km
        TF = TF | TF2;
    end
    O3_DU_interp_4int = O3_DU_interp(:,i);
    O3_DU_interp_4int(TF,:) = [];
    O3_DU_int(i) = sum(O3_DU_interp_4int)/2.69e20;% calculate O3 in DU value, 1DU = 2.69e20 molec/m^2
end

%     %%%%% H2O profiles %%%%%%%
%     
%     figure;hold all;
%     color_pool = colormap(jet(48));
%     for i =761:1:768
%         plot(Qv_interp(:,i)*1e6,H_grid,'Color',color_pool(i-720,:));% convert to ppmv
%     end
%     ylabel('GPH [km]');
%     xlabel('H_2O VMR [ppmv]');
%     ylim([0 top_height]);
%     colorbar;
%     title('MERRA2 H_2O Profiles [ppmv]');
%     print_setting(size_fig,save_fig,['MERRA2_H2O_1Dprofiles_' year]);
%     
    

%%%%% H2O total pwv time-series %%%%%%%
figure;
plot(time,pwv_int);
datetick('x','mmm-dd','keeplimits');
ylabel('H_2O pwv [mm]');
xlabel(year);
title('MERRA2 H_2O pwv');
print_setting(size_fig,save_fig,['MERRA2_H2O_pwv_' year]);

%%%%% O3 total DU time-series %%%%%%%
figure;
plot(time,O3_DU_int);
datetick('x','mmm-dd','keeplimits');
ylabel('O_3 [DU]');
xlabel(year);
title('MERRA2 O_3 DU');
print_setting(size_fig,save_fig,['MERRA2_O3_DU_' year]);

if plot_2D_profiles_time_series == 1
    %%%%% Temperature profiles %%%%%%%
    figure;
    imagesc(time,H_grid,T_interp);
    set(gca,'YDir','normal');
    datetick('x','mmm-dd','keeplimits');
    ylabel('GPH [km]');
    xlabel(year);
    ylim([0 top_height]);
    colorbar;
    title('MERRA2 Temperature Profiles [K]');
    
    hold all;
    z = [1.6e-4];% inner edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[0 0 0]);
    z = [1.2e-4];% outer edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[1 1 1]);
    caxis([180 290]);
    print_setting(size_fig,save_fig,['MERRA2_T_profiles_' year]);
    
    %%%%% PV profiles %%%%%%%
    figure;
    imagesc(time,H_grid,PV_interp.*1e6);
    set(gca,'YDir','normal');
    datetick('x','mmm-dd','keeplimits');
    ylabel('GPH [km]');
    xlabel(year);
    ylim([0 top_height]);
    caxis([0 15]);
    colorbar;
    title('MERRA2 PV Profiles [PVU]');
    print_setting(size_fig,save_fig,['MERRA2_PV_profiles_' year]);
    
    %%%%% O3 profiles %%%%%%%
    figure;
    imagesc(time,H_grid,O3_interp.*1e9);% plot VMR in ppb (1e-9) level
    set(gca,'YDir','normal');
    datetick('x','mmm-dd','keeplimits');
    ylabel('GPH [km]');
    xlabel(year);
    ylim([0 top_height]);
    caxis([100 8000]);
    colorbar;
    title('MERRA2 O3 Profiles [ppbv]');
    print_setting(size_fig,save_fig,['MERRA2_O3_profiles_' year]);
    
    hold all;
    z = [1.6e-4];% inner edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[0 0 0]);
    z = [1.2e-4];% outer edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[1 1 1]);
    caxis([0 2e4]);
    print_setting(size_fig,save_fig,['MERRA2_O3_profiles_' year]);
    
    %%%%% Theta profiles %%%%%%%
    figure;
    imagesc(time,H_grid,Theta_interp);
    set(gca,'YDir','normal');
    datetick('x','mmm-dd','keeplimits');
    ylabel('GPH [km]');
    xlabel(year);
    ylim([0 top_height]);
    %caxis([0 1000]);
    colorbar;
    title('MERRA2 Theta Profiles [K]');
    print_setting(size_fig,save_fig,['MERRA2_Theta_profiles_' year]);
    
    %%%%% sPV profiles %%%%%%%
    figure;
    imagesc(time,H_grid,sPV_interp);
    set(gca,'YDir','normal');
    datetick('x','mmm-dd','keeplimits');
    ylabel('GPH [km]');
    xlabel(year);
    ylim([0 top_height]);
    %caxis([0 1000]);
    colorbar;
    title('MERRA2 sPV Profiles [s^-^1]');
    
    hold all;% inner edge of Vortex
    z = [1.6e-4];
    contour(time,H_grid,sPV_interp,z,'LineColor',[1 1 1]);
    
    print_setting(size_fig,save_fig,['MERRA2_sPV_profiles_' year]);
    
    
    %%%%% dT_dZ profiles %%%%%%%
    figure;
    imagesc(time,H_grid,dT_dZ_interp);
    set(gca,'YDir','normal');
    datetick('x','mmm-dd','keeplimits');
    ylabel('GPH [km]');
    xlabel(year);
    ylim([0 top_height]);
    %caxis([0 1000]);
    colorbar;
    title('MERRA2 dT/dZ Profiles [K/km]');
    print_setting(size_fig,save_fig,['MERRA2_dT_dZ_profiles_' year]);
    
    hold all;
    z = [1.6e-4];% inner edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[0 0 0]);
    z = [1.2e-4];% outer edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[1 1 1]);
    caxis([-15 45]);
    print_setting(size_fig,save_fig,['MERRA2_dT_dZ_profiles_' year]);
    
    %%%%% Qv profiles %%%%%%%
    figure;
    imagesc(time,H_grid,Qv_interp.*1e6);% convert to ppmv
    set(gca,'YDir','normal');
    datetick('x','mmm-dd','keeplimits');
    ylabel('GPH [km]');
    xlabel(year);
    ylim([0 top_height]);
    %caxis([0 1000]);
    colorbar;
    title('MERRA2 H_2O Profiles [ppmv]');
    print_setting(size_fig,save_fig,['MERRA2_H2O_profiles_' year]);
    
    hold all;
    z = [1.6e-4];% inner edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[0 0 0]);
    z = [1.2e-4];% outer edge of Vortex
    contour(time,H_grid,sPV_interp,z,'LineColor',[1 1 1]);
    %caxis([-15 45]);
    print_setting(size_fig,save_fig,['MERRA2_H2O_profiles_' year]);
    
end

function combined_profiles = read_MERRA2_profiles(var,year)
data_dir = ['E:\H\work\MERRA\MERRA2_from_Sophie\' year];
cd(data_dir);
files = ls(['GEOS5MERRA2_DynVarsAtEureka_' year '_*.nc4']);
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


