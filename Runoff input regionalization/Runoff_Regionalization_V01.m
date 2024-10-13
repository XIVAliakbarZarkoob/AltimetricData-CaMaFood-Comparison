%% Aliakbar Zarkoob, AKA "XIV"
%  Gmail: XIV.Aliakbar.Zarkoob@gmail.com
%  Telegram: @XIVAliakbar

clc, clear, close all, beep off, format long g
set(0,'defaultTextInterpreter','latex')


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
    error('No Region Was Selected!')
end

%%

SAVE_PATH = ['./',REGION,'/'];
if ~exist(SAVE_PATH, 'dir')
    mkdir(SAVE_PATH);
end

[file,path,index] = uigetfile('*.nc*','Select Runoff Files','MultiSelect','on');
files = string();
if index == 1
    if iscell(file)
        for i = 1:size(file,2)

            file_i = char(file(i));
            data_i = rncdf([path,file_i]);
            if i == 1
                tmp = abs(data_i.lat-REGION_LAT(1));
                lat_idx(1) = find(min(tmp) == tmp,1,'last');
                tmp = abs(data_i.lat-REGION_LAT(2));
                lat_idx(2) = find(min(tmp) == tmp,1,'last');
                tmp = abs(data_i.lon-REGION_LON(1));
                lon_idx(1) = find(min(tmp) == tmp,1,'last');
                tmp = abs(data_i.lon-REGION_LON(2));
                lon_idx(2) = find(min(tmp) == tmp,1,'last');
                lat_num = max(lat_idx)-min(lat_idx)+1;
                lon_num = max(lon_idx)-min(lon_idx)+1;
            end
            time_num = size(data_i.time,1);

            ncFileName = [SAVE_PATH,file_i];
            nccreate(ncFileName, 'lat', 'Dimensions', {'lat', lat_num}, 'Datatype', 'single','Format','netcdf4');
            nccreate(ncFileName, 'lon', 'Dimensions', {'lon', lon_num}, 'Datatype', 'single','Format','netcdf4');
            nccreate(ncFileName, 'time', 'Dimensions', {'time', time_num}, 'Datatype', 'int64','Format','netcdf4');
            nccreate(ncFileName, 'Runoff', ...
                'Dimensions', {'lon', lon_num, 'lat', lat_num, 'time', time_num}, ...
                'Datatype', 'single');
            
            ncwrite(ncFileName, 'lat', data_i.lat(min(lat_idx):max(lat_idx)));
            ncwrite(ncFileName, 'lon', data_i.lon(min(lon_idx):max(lon_idx)));
            ncwrite(ncFileName, 'time', data_i.time);
            ncwrite(ncFileName, 'Runoff', data_i.Runoff(min(lon_idx):max(lon_idx),min(lat_idx):max(lat_idx),:));

            fprintf('%s Was Written.\n',file_i)
            clear data_i tmp

        end
    end 
end

% global attributes
% ncwriteatt(ncFileName, '/', 'title', 'Runoff Data'); % Example global attribute
