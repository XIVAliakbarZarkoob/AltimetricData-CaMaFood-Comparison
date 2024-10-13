%% Aliakbar Zarkoob, AKA "XIV"
%  Gmail: XIV.Aliakbar.Zarkoob@gmail.com
%  Telegram: @XIVAliakbar

clc, clear, close all, beep off, format long g
set(0,'defaultTextInterpreter','latex')

%#ok<*MINV>
%#ok<*AGROW>
%#ok<*NANMEAN>

%% Select Region

REGION = questdlg('Select Region','Region','Niger','Ganges-Brahmaputra','Niger');
if strcmp(REGION,'Niger')
    % Niger ------------------
    REGION_LAT = [2,20];
    REGION_LON = [-15,15];
    % ------------------------
elseif strcmp(REGION,'Ganges-Brahmaputra')
    % ganges-brahmaputra -----
    REGION_LAT = [20,33];
    REGION_LON = [72,100];
    % ------------------------
else 
    error('No regoin was selected!')
end

%% Load Altimetric Data

flag.main = string();
files.main = string();
NUM_SOURCES = 0;
k = 1;

% DAHITI
[file_d,path_d,index_d] = uigetfile('*.nc*','Select Dahiti Virtual Station Data Files','MultiSelect','on');
if index_d == 1
    NUM_SOURCES = NUM_SOURCES + 1;
    if iscell(file_d)
        for i = 1:size(file_d,2)
            ncid = netcdf.open([path_d,char(file_d(i))], 'NC_NOWRITE');
            files.main(k,1) = 'WL' + extractBefore(string(file_d(i)),'.nc');
            data.main.(files.main(k,1)).datetime = datetime(netcdf.getVar(ncid,0));
            data.main.(files.main(k,1)).water_level = double(netcdf.getVar(ncid,1));
            data.main.(files.main(k,1)).error = double(netcdf.getVar(ncid,2));
            data.main.(files.main(k,1)).latitude = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'latitude');
            data.main.(files.main(k,1)).longitude = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'longitude');
            flag.main(k,1) = 'DH';
            netcdf.close(ncid);
            k = k + 1;
        end
    else
        error('Select all the files or skip this data source by closing the input window.')
    end % DAHITI
end

% HYDROWEB
[file_h,path_h,index_h] = uigetfile('*.txt*','Select Hydroweb Virtual Station Data Files','MultiSelect','on');
if index_h == 1
    NUM_SOURCES = NUM_SOURCES + 1;
    if iscell(file_h)
        for i = 1:size(file_h,2)
            files.main(k,1) = extractBefore(string(file_h(i)),'.txt');
            if contains(files.main(k,1),'-')
                files.main(k,1) = strrep(files.main(k,1),'-','_');
            end
            data.main.(files.main(k,1)) = HWread([path_h,char(file_h(i))]);
            flag.main(k,1) = 'HW';
            k = k + 1;
        end
    else
        error('Select all the files or skip this data source by closing the input window.')
    end % HYDROWEB
end

% CLMS
[file_c,path_c,index_c] = uigetfile('*.json*','Select CLMS Virtual Station Data Files','MultiSelect','on');
if index_c == 1
    NUM_SOURCES = NUM_SOURCES + 1;
    if iscell(file_c)
        for i = 1:size(file_c,2)
            files.main(k,1) = extractBetween(string(file_c(i)),'c_gls_','_0');
            data.main.(files.main(k,1)) = jsondecode(fileread([path_c,char(file_c(i))]));
            data.main.(files.main(k,1)).identifier = string(vertcat(data.main.(files.main(k,1)).data.identifier));
            dt = vertcat(data.main.(files.main(k,1)).data.datetime); 
            dt = datetime([double(string(dt(:,1:4))), double(string(dt(:,6:7))), double(string(dt(:,9:10)))]);
            data.main.(files.main(k,1)).datetime = dt;
            data.main.(files.main(k,1)).wse = vertcat(data.main.(files.main(k,1)).data.orthometric_height_of_water_surface_at_reference_position);
            data.main.(files.main(k,1)).uncertainty = vertcat(data.main.(files.main(k,1)).data.associated_uncertainty);
            % data.main.(files.main(k,1)).satellite = string(vertcat(data.main.(files.main(k,1)).data.satellite));
            % data.main.(files.main(k,1)).ground_track_number = vertcat(data.main.(files.main(k,1)).data.ground_track_number);
            data.main.(files.main(k,1)) = rmfield(data.main.(files.main(k,1)),'data');
            flag.main(k,1) = 'CLMS';
            k = k + 1;
        end
    else
        error('Select all the files or skip this data source by closing the input window.')
    end % CLMS
end

if NUM_SOURCES == 0
    error('No data from any resources were selected!')
end

%% Timeseries Range & Virtual Station Coordinates

years.start = zeros(length(files.main),1);
years.end = zeros(length(files.main),1);
for i = 1:length(files.main)
    data_i =  data.main.(files.main(i));
    if strcmp(flag.main(i), 'DH')
        years.start(i) = data_i.datetime(1).Year;
        years.end(i) = data_i.datetime(end).Year;
    elseif strcmp(flag.main(i), 'HW')
        years.start(i) = data_i.FirstDate.Year;
        years.end(i) = data_i.LastDate.Year;
    elseif strcmp(flag.main(i), 'CLMS')
        years.start(i) = data_i.datetime(1).Year;
        years.end(i) = data_i.datetime(end).Year;
    end
end

%% 2002 to 2008 Data Matrix

k = 1;
flag.f2002t2008 = string();
for i = 1:length(files.main)
    data_i =  data.main.(files.main(i));
    if years.start(i)==2002
        files.f2002t2008(k,1) = files.main(i);
        flag.f2002t2008(k,1) = flag.main(i);
        k = k+1;
    end
end

[X.f2002t2008, time.f2002t2008, coords.f2002t2008] = DataMatrix(2002, 2009, files.f2002t2008, data.main, flag.f2002t2008);

rmv_idx = [];
for i = 1:size(X.f2002t2008,1)
    gap_flag = MonthlyGap(X.f2002t2008(i,:),7);
    if gap_flag
        rmv_idx = [rmv_idx; i]; 
    end
end
X.f2002t2008(rmv_idx,:) = [];
coords.f2002t2008(rmv_idx,:) = [];
files.f2002t2008(rmv_idx) = [];
flag.f2002t2008(rmv_idx) = [];

SFill = 7; 
EFill = size(time.f2002t2008,1);
time.f2002t2008 = time.f2002t2008(SFill:EFill);
X.f2002t2008 = X.f2002t2008(:,SFill:EFill);
Xf.f2002t2008 = fillmissing(X.f2002t2008,'linear',2);

% heatmap(X.f2002t2008)

%% 2008 to 2016 Data Matrix

k = 1;
flag.f2008t2016 = string();
for i = 1:length(files.main)
    data_i =  data.main.(files.main(i));
    if years.start(i)==2008
        files.f2008t2016(k,1) = files.main(i);
        flag.f2008t2016(k,1) = flag.main(i);
        k = k+1;
    end
end

[X.f2008t2016, time.f2008t2016, coords.f2008t2016] = DataMatrix(2008, 2016, files.f2008t2016, data.main, flag.f2008t2016);

rmv_idx = [];
for i = 1:size(X.f2008t2016,1)
    gap_flag = MonthlyGap(X.f2008t2016(i,:),6);
    if gap_flag
        rmv_idx = [rmv_idx; i]; 
    end
end
X.f2008t2016(rmv_idx,:) = [];
coords.f2008t2016(rmv_idx,:) = [];
files.f2008t2016(rmv_idx) = [];
flag.f2008t2016(rmv_idx) = [];

SFill = 7; 
EFill = size(time.f2008t2016,1);
time.f2008t2016 = time.f2008t2016(SFill:EFill);
X.f2008t2016 = X.f2008t2016(:,SFill:EFill);
Xf.f2008t2016 = fillmissing(X.f2008t2016,'linear',2);

% heatmap(X.f2008t2016)

%% 2016 to 2019 Data Matrix

k = 1;
flag.f2016t2019 = string();
for i = 1:length(files.main)
    data_i =  data.main.(files.main(i));
    if years.start(i)==2016
        files.f2016t2019(k,1) = files.main(i);
        flag.f2016t2019(k,1) = flag.main(i);
        k = k+1;
    end
end

[X.f2016t2019, time.f2016t2019, coords.f2016t2019] = DataMatrix(2016, 2018, files.f2016t2019, data.main, flag.f2016t2019);

rmv_idx = [];
for i = 1:size(X.f2016t2019,1)
    gap_flag = MonthlyGap(X.f2016t2019(i,:),5);
    if gap_flag
        rmv_idx = [rmv_idx; i]; 
    end
end
X.f2016t2019(rmv_idx,:) = [];
coords.f2016t2019(rmv_idx,:) = [];
files.f2016t2019(rmv_idx) = [];
flag.f2016t2019(rmv_idx) = [];

SFill = 5; 
EFill = size(time.f2016t2019,1);
time.f2016t2019 = time.f2016t2019(SFill:EFill);
X.f2016t2019 = X.f2016t2019(:,SFill:EFill);
Xf.f2016t2019 = fillmissing(X.f2016t2019,'linear',2);

% heatmap(X.f2016t2019);

%% 2019 to 2024 Data Matrix

k = 1;
flag.f2019t2024 = string();
for i = 1:length(files.main)
    data_i =  data.main.(files.main(i));
    if (years.start(i)==2016 || years.start(i)==2019 || years.start(i)==2018) && years.end(i)==2024
        files.f2019t2024(k,1) = files.main(i);
        flag.f2019t2024(k,1) = flag.main(i);
        k = k+1;
    end
end

[X.f2019t2024, time.f2019t2024, coords.f2019t2024] = DataMatrix(2018, 2024, files.f2019t2024, data.main, flag.f2019t2024);

rmv_idx = [];
for i = 1:size(X.f2019t2024,1)
    gap_flag = MonthlyGap(X.f2019t2024(i,:),5);
    if gap_flag
        rmv_idx = [rmv_idx; i]; 
    end
end
X.f2019t2024(rmv_idx,:) = [];
coords.f2019t2024(rmv_idx,:) = [];
files.f2019t2024(rmv_idx) = [];
flag.f2019t2024(rmv_idx) = [];

SFill = 1; 
EFill = 80;
time.f2019t2024 = time.f2019t2024(SFill:EFill);
X.f2019t2024 = X.f2019t2024(:,SFill:EFill);
Xf.f2019t2024 = fillmissing(X.f2019t2024,'linear',2);

% heatmap(X.f2019t2024)

%% Total VS in separated years

[files.total,idx] = unique([files.f2002t2008;files.f2008t2016;files.f2016t2019;files.f2019t2024]);
flag.total = [flag.f2002t2008;flag.f2008t2016;flag.f2016t2019;flag.f2019t2024];
flag.total = flag.total(idx);
coords.total = table(zeros(length(files.total),1),zeros(length(files.total),1));
coords.total.Properties.VariableNames = {'lat','lon'};

for i = 1:length(files.total)
    if strcmp(flag.total(i),'DH')
        data_i = data.main.(files.total(i));
        coords.total.lat(i) = data_i.latitude;
        coords.total.lon(i) = data_i.longitude;
    elseif strcmp(flag.total(i),'HW')
        data_i = data.main.(files.total(i)).MainData;
        tmp = data_i.("LATITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=90);
        coords.total.lat(i) = mean(tmp);
        tmp = data_i.("LONGITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=180);
        coords.total.lon(i) = mean(tmp);
    elseif strcmp(flag.total(i),'CLMS')
        data_i = data.main.(files.total(i));
        coords.total.lat(i) = data_i.geometry.coordinates(2);
        coords.total.lon(i) = data_i.geometry.coordinates(1);
    end
end

%% Load CaMa-Flood Outputs (WSE & Q) 

[file_WSE,path_WSE] = uigetfile('*.nc*','Select CaMa-Flood output files (WSE/sfcelv)','MultiSelect','on');
if iscell(file_WSE)
    fscmf.WSE = string();
    for i = 1:size(file_WSE,2)
        ncid = netcdf.open([path_WSE,char(file_WSE(i))], 'NC_NOWRITE');
        fscmf.WSE(i,1) = extractBefore(string(file_WSE(i)),'.nc');
        cmf.WSE.(fscmf.WSE(i,1)).lat = netcdf.getVar(ncid,0);
        cmf.WSE.(fscmf.WSE(i,1)).lon = netcdf.getVar(ncid,1);
        cmf.WSE.(fscmf.WSE(i,1)).time = netcdf.getVar(ncid,2);
        cmf.WSE.(fscmf.WSE(i,1)).sfcelv = netcdf.getVar(ncid,3);
        netcdf.close(ncid);
    end
else
    error('Select all the WSE/sfcelv files!')
end

[file_Q,path_Q] = uigetfile('*.nc*','Select CaMa-Flood output files (Q/rivout)','MultiSelect','on');
if iscell(file_Q)
    fscmf.Q = string();
    for i = 1:size(file_Q,2)
        ncid = netcdf.open([path_Q,char(file_Q(i))], 'NC_NOWRITE');
        fscmf.Q(i,1) = extractBefore(string(file_Q(i)),'.nc');
        cmf.Q.(fscmf.Q(i,1)).lat = netcdf.getVar(ncid,0);
        cmf.Q.(fscmf.Q(i,1)).lon = netcdf.getVar(ncid,1);
        cmf.Q.(fscmf.Q(i,1)).time = netcdf.getVar(ncid,2);
        cmf.Q.(fscmf.Q(i,1)).rivout = netcdf.getVar(ncid,3);
    end
else
    error('Select all the Q/rivout files!')
end

if length(fscmf.Q) ~= length(fscmf.WSE)
    error('Output files for Q and WSE must be the same number!')
end

%% Mask CaMa-Flood Outputs on Virtual Stations

startDate = datetime(double(extractAfter(fscmf.Q(1),'o_rivout')), 1, 1);
endDate = datetime(double(extractAfter(fscmf.Q(end),'o_rivout')), 12, 31);
time_cmf = (startDate:endDate)'; time_cmf = time_cmf(1:(time_cmf(end).Year-time_cmf(1).Year+1)*365);
X_Q = zeros(length(files.total),length(time_cmf));
X_WSE = zeros(length(files.total),length(time_cmf));
lat = zeros(length(files.total),1); lon = zeros(length(files.total),1);
idx_lat = lat; idx_lon = lon;
for i = 1:length(files.total)
    data_Q = cmf.Q.(fscmf.Q(1));
    data_WSE = cmf.WSE.(fscmf.WSE(1));
    tmp = find(min(abs(data_Q.lat-coords.total.lat(i))) == abs(data_Q.lat-coords.total.lat(i)));
    idx_Q(1) = tmp(1);
    tmp = find(min(abs(data_Q.lon-coords.total.lon(i))) == abs(data_Q.lon-coords.total.lon(i)));
    idx_Q(2) = tmp(1);
    % idx_WSE(1) = find(min(abs(data_WSE.lat-coords.vs.f2019t2024.lat(i))) == abs(data_WSE.lat-coords.vs.f2019t2024.lat(i)));
    % idx_WSE(2) = find(min(abs(data_WSE.lon-coords.vs.f2019t2024.lon(i))) == abs(data_WSE.lon-coords.vs.f2019t2024.lon(i)));
    lat(i) = double(data_Q.lat(idx_Q(1))); lon(i) = double(data_Q.lon(idx_Q(2)));
    idx_lat(i) = idx_Q(1); idx_lon(i) = idx_Q(2);
    for j = 1:length(fscmf.Q)
        tmp = reshape(cmf.Q.(fscmf.Q(j)).rivout(idx_Q(2),idx_Q(1),:),[],1); tmp = tmp(1:365);
        X_Q(i,j*365-364:j*365) = tmp;
        tmp = reshape(cmf.WSE.(fscmf.WSE(j)).sfcelv(idx_Q(2),idx_Q(1),:),[],1); tmp = tmp(1:365);
        X_WSE(i,j*365-364:j*365) = tmp;
    end
end
coords.cmf = table(lat,lon,idx_lat,idx_lon);

if sum(contains(fscmf.Q,'2024'))==1
    X_Q = X_Q(:,200:end-300);
    X_WSE = X_WSE(:,200:end-300);
    time_cmf = time_cmf(200:end-300);
end

%% PCA on Altimetric Data 

[E.f2002t2008, P.f2002t2008, eigvalue.f2002t2008, ~, Pers.f2002t2008, ~] = pca(Xf.f2002t2008',"NumComponents",1);
E.f2002t2008 = repmat(std(P.f2002t2008),size(E.f2002t2008,1),1).*E.f2002t2008;
P.f2002t2008 = P.f2002t2008./repmat(std(P.f2002t2008),size(P.f2002t2008,1),1);

[E.f2008t2016, P.f2008t2016, eigvalue.f2008t2016, ~, Pers.f2008t2016, ~] = pca(Xf.f2008t2016',"NumComponents",1);
E.f2008t2016 = repmat(std(P.f2008t2016),size(E.f2008t2016,1),1).*E.f2008t2016;
P.f2008t2016 = P.f2008t2016./repmat(std(P.f2008t2016),size(P.f2008t2016,1),1);

[E.f2016t2019, P.f2016t2019, eigvalue.f2016t2019, ~, Pers.f2016t2019, ~] = pca(Xf.f2016t2019',"NumComponents",1);
E.f2016t2019 = repmat(std(P.f2016t2019),size(E.f2016t2019,1),1).*E.f2016t2019;
P.f2016t2019 = P.f2016t2019./repmat(std(P.f2016t2019),size(P.f2016t2019,1),1);

[E.f2019t2024, P.f2019t2024, eigvalue.f2019t2024, ~, Pers.f2019t2024, ~] = pca(Xf.f2019t2024',"NumComponents",1);
E.f2019t2024 = repmat(std(P.f2019t2024),size(E.f2019t2024,1),1).*E.f2019t2024;
P.f2019t2024 = P.f2019t2024./repmat(std(P.f2019t2024),size(P.f2019t2024,1),1);

%% CaMa-Flood WSE & Q PCA

[E_Q, P_Q, eigvalue_Q, ~, Pers_Q, ~] = pca(X_Q',"NumComponents",1);
E_Q = repmat(std(P_Q),size(E_Q,1),1).*E_Q;
P_Q = P_Q./repmat(std(P_Q),size(P_Q,1),1);

[E_WSE, P_WSE, eigvalue_WSE, ~, Pers_WSE, ~] = pca(X_WSE',"NumComponents",1);
E_WSE = repmat(std(P_WSE),size(E_WSE,1),1).*E_WSE;
P_WSE = P_WSE./repmat(std(P_WSE),size(P_WSE,1),1);

%% Time Series of Some Virtual Stations

% [~,idx] = ismember(files.f2002t2008,files.total);
% Xm.f2002t2008 = (X.f2002t2008-nanmean(X.f2002t2008,2)) + mean(X_WSE(idx,:),2); 
% [~,idx] = ismember(files.f2008t2016,files.total);
% Xm.f2008t2016 = (X.f2008t2016-nanmean(X.f2008t2016,2)) + mean(X_WSE(idx,:),2); 
% [~,idx] = ismember(files.f2016t2019,files.total);
% Xm.f2016t2019 = (X.f2016t2019-nanmean(X.f2016t2019,2)) + mean(X_WSE(idx,:),2); 
% [~,idx] = ismember(files.f2019t2024,files.total);
% Xm.f2019t2024 = (X.f2019t2024-nanmean(X.f2019t2024,2)) + mean(X_WSE(idx,:),2); 

k = 1; kk = 1; save_path = './Figures/';
vs_range = 1:35:size(files.total,1); vs_range = vs_range(1:end-mod(length(vs_range),4));
coords.tsc = table(zeros(length(vs_range),1),zeros(length(vs_range),1));
coords.tsc.Properties.VariableNames = {'lat','lon'};
for i = vs_range
    
    if mod(k,4) == 1
        figure()
        % figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        sgtitle('Timeseries of Some VS Compared to CaMa-Flood Outputs')
    end
    data_i = data.main.(files.total(i));
    subplot(2,2,kk)
    title(sprintf('VS%g',k))
    hold on; grid on; box on
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.TickLabelInterpreter = 'latex';

    yyaxis left
    ylabel('WSE [$m$]')
    tmp = mean(X_WSE(i,:));
    if strcmp(flag.total(i),'DH')
        wse = data_i.water_level;
        wse_m = wse - mean(wse) + tmp;
        plot(data_i.datetime, wse_m, '.-k', 'LineWidth', 2, 'MarkerSize', 10)
        xlim([data_i.datetime(1) data_i.datetime(end)])
        coords.tsc.lat(k) = data_i.latitude;
        coords.tsc.lon(k) = data_i.longitude;
    elseif strcmp(flag.total(i),'HW')
        wse = data_i.MainData.("ORTHOMETRIC HEIGHT (M) OF WATER SURFACE AT REFERENCE POSITION");
        wse_m = wse - mean(wse) + tmp;
        plot(data_i.MainData.Time, wse_m, '.-k', 'LineWidth', 2, 'MarkerSize', 10)
        xlim([data_i.MainData.Time(1) data_i.MainData.Time(end)])
        tmp = data_i.MainData.("LATITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=90);
        coords.tsc.lat(k) = mean(tmp);
        tmp = data_i.MainData.("LONGITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=180);
        coords.tsc.lon(k) = mean(tmp);
    elseif strcmp(flag.total(i),'CLMS')
        wse = data_i.wse;
        wse_m = wse - mean(wse) + tmp;
        plot(data_i.datetime, wse_m, '.-k', 'LineWidth', 2, 'MarkerSize', 10)
        xlim([data_i.datetime(1) data_i.datetime(end)])
        coords.tsc.lat(k) = data_i.geometry.coordinates(2);
        coords.tsc.lon(k) = data_i.geometry.coordinates(1);
    end
    plot(time_cmf, X_WSE(i,:), '-b', 'LineWidth',1.2)
    yyaxis right
    ylabel('Q [$\frac{m^3}{s}$]')
    plot(time_cmf, X_Q(i,:), '-r', 'LineWidth',1.2)
    xlabel('Time');
    legend('Altimetric WSE', 'CMF WSE', 'CMF Q', 'Interpreter', 'latex')

    k = k+1; kk = kk+1;
    if mod(k,4) == 1
        % filename = sprintf('%s_TimeSeries_Compare_%d', REGION, ceil(k/4)-1); 
        % filename = [save_path,filename];
        % saveas(gcf, [filename '.png']);      
        % saveas(gcf, [filename '.svg']); 
        kk = 1;
    end

end

figure()
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
title('Location of Compared Virtual Stations','FontSize',15)
geoplot(coords.tsc.lat,coords.tsc.lon,'^b','MarkerFaceColor','b','MarkerSize',10)
text(coords.tsc.lat,coords.tsc.lon+0.5,strcat("VS", string(1:length(vs_range))),'BackgroundColor','w')
geolimits(REGION_LAT,REGION_LON)


%% Temporal Components

figure()
hold on; grid on; box on; axis tight
ax = gca;
ax.XAxis.FontSize = 12;
ax.XAxis.TickLabelInterpreter = 'latex';

xlabel('Time','FontSize',20)

yyaxis left
plot(time.f2002t2008,P.f2002t2008,'o-k','LineWidth',1.5,'MarkerSize',10)
plot(time.f2008t2016,P.f2008t2016,'+-k','LineWidth',1.5,'MarkerSize',10)
plot(time.f2016t2019,P.f2016t2019,'x-k','LineWidth',1.5,'MarkerSize',10)
plot(time.f2019t2024,P.f2019t2024,'*-k','LineWidth',1.5,'MarkerSize',10)
plot(time_cmf,P_WSE(:,1),'-b','LineWidth',1.5,'MarkerSize',15)
ylabel('$WSE$','FontSize',20)

yyaxis right
plot(time_cmf,P_Q(:,1),'-r','LineWidth',1.5,'MarkerSize',15)
ylabel('$Q$','FontSize',20)

title('1st Temporal Components of Altimetric Data Compared to CaMa-Flood WSE and Q','FontSize',15)
legend('Altimetric WSE - 2002 to 2008','Altimetric WSE - 2008 to 2016', ...
    'Altimetric WSE - 2016 to 2019','Altimetric WSE - 2019 to 2024', ...
    'CMF WSE','CMF Q', 'Interpreter', 'latex','FontSize',10)
xticks(datetime(2001, 1, 1):calyears(1):datetime(2024, 1, 1))


%% ------------------------------ Functions -----------------------------------

% Create Data Matrix
function [X, time, coords] = DataMatrix(StartYear, EndYear, files, data, flag)

    time = [repelem((StartYear:EndYear)',12), ...
            repmat((1:12)',length(StartYear:EndYear)',1), ...
            ones((EndYear-StartYear+1)*12,1)*15];
    time = datetime(time);
    X = NaN(length(files),length(time));
    
    lat = zeros(length(files),1);
    lon = zeros(length(files),1);
    coords = table(lat,lon);
    
    for i = 1:length(files)
        if strcmp(flag(i),'DH')
            data_i = data.(files(i));
            for j = 1:length(time)
                idx = (time.Year(j) == data_i.datetime.Year) + (time.Month(j) == data_i.datetime.Month);
                idx = idx == 2;
                X(i,j) = mean(data_i.water_level(idx));
            end
            coords.lat(i) = data_i.latitude;
            coords.lon(i) = data_i.longitude;
        elseif strcmp(flag(i),'HW')
            data_i = data.(files(i)).MainData;
            for j = 1:length(time)
                idx = (time.Year(j) == data_i.Time.Year) + (time.Month(j) == data_i.Time.Month);
                idx = idx == 2;
                X(i,j) = mean(data_i.("ORTHOMETRIC HEIGHT (M) OF WATER SURFACE AT REFERENCE POSITION")(idx));
            end
            tmp = data_i.("LATITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=90);
            coords.lat(i) = mean(tmp);
            tmp = data_i.("LONGITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=180);
            coords.lon(i) = mean(tmp);
        elseif strcmp(flag(i),'CLMS')
            data_i = data.(files(i));
            for j = 1:length(time)
                idx = (time.Year(j) == data_i.datetime.Year) + (time.Month(j) == data_i.datetime.Month);
                idx = idx == 2;
                X(i,j) = mean(data_i.wse(idx));
            end
            coords.lat(i) = data_i.geometry.coordinates(2);
            coords.lon(i) = data_i.geometry.coordinates(1);
        end
    end

end

% Monthly Gap
function gap_flag = MonthlyGap(data,threshold)
    nan_idx = find(isnan(data));
    count = 1;
    gap_flag = false;
    for i = 2:length(nan_idx)
        if nan_idx(i) == nan_idx(i-1) + 1
            count = count + 1;
        else
            count = 1;
        end
        if count > threshold
            gap_flag = true;
            break;
        end
    end
end

% Read Hydroweb Data
function data = HWread(filepath)

data = rheader(filepath);
MainData = readtimetable(filepath,'NumHeaderLines',data.EOH,'TextType','string');
MainData.Var4 = [];
time = split(MainData.Var1,':');
dt = datetime(double([MainData.Time.Year,MainData.Time.Month,MainData.Time.Day,time(:,1),time(:,2),zeros(length(time),1)]));
MainData.Time = dt;
MainData.Var1 =[];
MainData.Properties.VariableNames = data.Content(3:end);
data.MainData = MainData;

function header = rheader(filepath)
        FileID = fopen(filepath);
        header.EOH = 0;
        while true
            header.EOH = header.EOH + 1;
            line = fgetl(FileID);
            if contains(line,'#############')
                fclose(FileID);
                break
            end
            if contains(line,'#BASIN')
                header.Basin = strtrim(line(9:end));
            end
            if contains(line,'#RIVER')
                header.River = strtrim(line(9:end));
            end
            if contains(line,'#ID')
                header.ID = strtrim(line(6:end));
            end
            if contains(line,'#TRIBUTARY OF')
                header.Tributary = strtrim(line(16:end));
            end
            if contains(line,'#APPROX. WIDTH OF REACH (m)')
                header.ApproxWidthOfReach = str2double(line(30:end));
            end
            if contains(line,'SURFACE OF UPSTREAM WATERSHED (km2)')
                header.SurfaceOfUpstreamWatershed = str2double(line(39:end));
            end
            if contains(line,'#RATING CURVE PARAMETERS')
                header.RatingCurveParameters = sscanf(line(67:end),'%f');
                if isempty(header.RatingCurveParameters)
                    header.RatingCurveParameters = NaN;
                end
            end
            if contains(line,'#REFERENCE ELLIPSOID')
                header.Refrence.Ellipsoid = strtrim(line(23:end));
            end
            if contains(line,'#REFERENCE LONGITUDE')
                header.Refrence.Longitude = str2double(line(23:end));
            end
            if contains(line,'#REFERENCE LATITUDE')
                header.Refrence.Latitude = str2double(line(22:end));
            end
            if contains(line,'#REFERENCE DISTANCE (km)')
                header.Refrence.Distance = str2double(line(27:end));
            end
            if contains(line,'#GEOID MODEL')
                header.GeoidModel = strtrim(line(15:end));
            end
            if contains(line,'#GEOID ONDULATION AT REF POSITION(M.mm)')
                header.GeoidOnulation = str2double(line(42:end));
            end
            if contains(line,'##MISSION(S)-TRACK(S)')
                header.MissionTrack = strtrim(line(23:end));
            end
            if contains(line,'#STATUS')
                header.Status = strtrim(line(10:end));
            end
            if contains(line,'#VALIDATION CRITERIA')
                header.ValidationCriteria = strtrim(line(23:end));
            end
            if contains(line,'#MEAN ALTITUDE(M.mm)')
                header.MeanAltitude = str2double(line(23:end));
            end
            if contains(line,'#MEAN SLOPE (mm/km)')
                header.MeanSlope = str2double(line(22:end));
            end
            if contains(line,'#NUMBER OF MEASUREMENTS IN DATASET')
                header.MeasurmentsNum = str2double(line(37:end));
            end
            if contains(line,'#FIRST DATE IN DATASET')
                header.FirstDate = datetime(sscanf(strrep(line(25:end),'-',' '),'%f')');
            end
            if contains(line,'#LAST DATE IN DATASET')
                header.LastDate = datetime(sscanf(strrep(line(24:end),'-',' '),'%f')');
            end
            if contains(line,'#DISTANCE MIN IN DATASET (km)')
                header.MinDistance = str2double(line(32:end));
            end
            if contains(line,'#DISTANCE MAX IN DATASET (km)')
                header.MaxDistance = str2double(line(32:end));
            end
            if contains(line,'#PRODUCTION DATE')
                header.ProductionDate = datetime(sscanf(strrep(line(19:end),'-',' '),'%f')');
            end
            if contains(line,'#PRODUCT VERSION')
                header.ProductionVersion = strtrim(line(19:end));
            end
            if contains(line,'#PRODUCT CITATION')
                header.ProductCitation = strtrim(line(20:end));
            end
            if contains(line,'#SOURCES')
                header.Sources = strtrim(line(11:end));
                if isempty(header.Sources)
                    header.Sources = NaN;
                end
            end
            if contains(line,'#COL')
                cnum = str2double(extractBefore(line(5:end),' : '));
                if cnum == 1
                    header.Content = "";
                end
                header.Content(cnum,1) = extractAfter(line,' : ');
            end
        end
    end
end
