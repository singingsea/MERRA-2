function Data = plot_MERRA2_WMO_v2(year,var)
% this can plot MERRA2 WMO thermal tropopause information
    save_fig = 0;
    %var = 'WMO_Tropopauses';
    %var = 'Temp_at_WMO_Tropopause';
    %var = 'Press_at_WMO_Tropopause';
    %var = 'Theta_at_WMO_Tropopause';
    %var = 'WMOSurf';
    if year == 2012 || year == 2008 || year == 2016
        day =1:1:366;
    elseif year == 2017
        day =1:1:304; % this change is due to we only have 2017 data until end of Oct 2017!
    else
        day = 1:1:365;
    end
    f1 = figure;hold all;
    %path = 'H:\work\MERRA\MERRA2\';
    %path = ['H:\work\MERRA\MERRA2_from_Sophie\' num2str(year) '\'];
    path = ['E:\H\work\MERRA\MERRA2_from_Sophie\' num2str(year) '\'];
    cd(path);
    %Year_2dig = num2str(year);
    lc = dir([path,'GEOS5MERRA2_TropopausesAtEureka*.nc4']); %read in all MERRA2 data
    D = char(lc.name); % build file name table
    N=size(D);
    fid = fopen([var '_' num2str(year) '.txt'],'w+');
    fclose(fid);
    %fprintf(fid,'%s %s %s\n', 'f_day', 'Dyn_Trops', 'PV index');
    
    for j =1:1:N
        %WMO_Trop = ncread('GEOS5MERRA2_TropopausesAtEureka_2013-02-01_2013-04-30_000000.nc4','WMO_Tropopauses');
        %WMO_PV = ncread('GEOS5MERRA2_TropopausesAtEureka_2013-02-01_2013-04-30_000000.nc4','PV_at_WMO_Tropopause');
        Dyn_Trop(:,:,j) = ncread(D(j,:),var);
        time = 0:3:21;
        f_time = time/24;
        f_day = day + f_time(j);
        f_days = f_day';
        UTC = (datenum(year,1,1)+f_time(j)): 1 :( (datenum(year,1,1)+f_time(j)) + max(day)-1);
        UTC = UTC';
        %[date, month] = Julian2Date(year,day);
        %UTC = datenum([year,month,date,time(j),0,0]);
        for k=1:1:1
            Dyn_Trops = Dyn_Trop(1,:,j)';            
            %figure(f1);
            %plot(f_day,Dyn_Trops,'.');
            N_size = size(Dyn_Trops);          
            k_s = repmat(k,N_size(1),1);
            %fprintf(fid,'%f %f %f\n', f_days, Dyn_Trops, k_s);
            Output = [UTC,Dyn_Trops,k_s];
            dlmwrite([var '_' num2str(year) '.txt'],Output,'-append','delimiter',' ','precision',10);
        end
               
    end
    
    Data = dlmread([var '_' num2str(year) '.txt']);
    gscatter(Data(:,1),Data(:,2),Data(:,3));
    %xlabel('Day of the year 2013');
    xlabel(['Day of the year ' num2str(year)]);
    %ylabel('Height (km)');
    ylabel([var ' km']);
    legend(['WMO ' var]);
    %print_setting(1/2,save_fig,var);
    print_setting(1/2,save_fig,[var '_' num2str(year) ]);
    
    
    