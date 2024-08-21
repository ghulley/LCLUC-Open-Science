%#ok<*UNRCH,*ALIGN,*NOPRT,*NASGU,*AGROW,*NBRAK,*NBRAK2>
% % inputs example: yearStr = '2021'; startDayStr = '04/01'; endDayStr = '08/01';
% % outputLoc=['./']
% % datLocationTemplate=['./Input_Data/yearStr,/Month/Day'];
% % this code assumes that the input filenames have this format: 'VNP21IMG_NRT.AYYYYDOY.HHmm.002.<processing_date_stamp>.nc' and the files are stored like the line above
function [] = VNP21IMG_LST_Cutouts(yearStr, startDayStr, endDayStr)

    addpath('./findMinPointIn2D/')

    narginchk(3,3);
    
    % % use: 
    % % #1 for distance() % MATLAB map toolbox required    
    % % #2 halversine % base MatLab no extra toolbox 
    % % #3 vicenty % base MatLab no extra toolbox
    % % #1 is prefered even though extensive testing of #2 and #3 showed an equivalent number of matches for 72 scenes
    % % #3 is more expensive but more accurate than #2 
    geoLocateMinimumDistance = 3;

    if geoLocateMinimumDistance == 3
        maxIter = 21;
	else
		earthRadiusInMeters = 6371000;
    end

    percentZeroThreshold = 40;

% %   Only one use_<TemperatureUnit> must be true
% %   Anything else will default to Kelvin 

    use_Celsius = false;
%     use_Celsius = true;

    use_Fahrenheit = false;
%     use_Fahrenheit = true;

    use_Kelvin = false;
%     use_Kelvin = true; 

    ch = [use_Celsius, use_Fahrenheit, use_Kelvin];
    ch = ch == true;
    ch_cnt = sum(ch);    

    if ch_cnt ~= 1
        use_Celsius = false;
        use_Fahrenheit = false;
        use_Kelvin = true;
    end

% %   Required Boolean Flags 
    writeTifs = false;
%     writeTifs = true;

    writeImages = false;
%     writeImages = true;

%     displayImages = false; 
    displayImages = true;    

    datLocationTemplate=['./Input_Data/',yearStr,'/'];  
	
	outputLoc=['./'];

    set(0,'DefaultFigureColormap',jet);    
    
    forwardSlashPos = strfind(startDayStr, '/');
    startDayLogStr = startDayStr;
    startDayLogStr(forwardSlashPos) = [];
    forwardSlashPos = strfind(endDayStr, '/');
    endDayLogStr = endDayStr;
    endDayLogStr(forwardSlashPos) = [];
    fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'w');
    fclose(fid);

    fprintf(['./VNP21IMG_LST_Cutouts(yearStr=''', yearStr,''', startDayStr=''' startDayStr,''', endDayStr=''', endDayStr,''')', newline]);
    fprintf(['yearStr = ', yearStr, newline]); 
    fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
    fwrite(fid, ['./VNP21IMG_LST_Cutouts(yearStr=''', yearStr,''', startDayStr=''' startDayStr,''', endDayStr=''', endDayStr,''')', newline]);
    fclose(fid);

    if ch_cnt ~= 1
        fprintf(['./LST can only be in one temperature unit ... defaulting to Kelvin ... ', newline]);
        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
        fwrite(fid, ['./LST can only be in one temperature unit ... defaulting to Kelvin ... ', newline]);
        fclose(fid);
    end

	fprintf(['use_Celsius = ', num2str(use_Celsius), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['use_Celsius = ', num2str(use_Celsius), newline]);
    fclose(fid);
    
	fprintf(['use_Fahrenheit = ', num2str(use_Fahrenheit), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['use_Fahrenheit = ', num2str(use_Fahrenheit), newline]);
    fclose(fid);
    
	fprintf(['use_Kelvin = ', num2str(use_Kelvin), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['use_Kelvin = ', num2str(use_Kelvin), newline]);
    fclose(fid);

	fprintf(['writeTifs = ', num2str(writeTifs), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['writeTifs = ', num2str(writeTifs), newline]);
    fclose(fid);

	fprintf(['writeImages = ', num2str(writeImages), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['writeImages = ', num2str(writeImages), newline]);
    fclose(fid);
    
	fprintf(['displayImages = ', num2str(displayImages), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['displayImages = ', num2str(displayImages), newline]);
    fclose(fid);

    cnt = 0; 
       
    prevFile = '';    
    
    % number of cities
    citiesInfoBoolArrTemplate = zeros(1,6); 
    citiesInfoBoolArr = citiesInfoBoolArrTemplate;
    citiesInfoArr = cell(1,6);
    
    % set each city's settings 
    Delhi.Name='Delhi';
    Delhi.Center=[28.689006,77.100125];
    Delhi.NE=[29.149104,77.633174];
    Delhi.SW=[28.228908,76.567076];
    Delhi.DayHour=[03,12];
    Delhi.DayMinute=[00,59];
    Delhi.NightHour=[18,00];
    Delhi.NightMinute=[00,59];
    citiesInfoArr{1} = Delhi;
    
    Ahmedabad.Name='Ahmedabad';
    Ahmedabad.Center=[23.021781,72.571285];
    Ahmedabad.NE=[23.481879,73.104334];
    Ahmedabad.SW=[22.561683,72.038236];
    Ahmedabad.DayHour=[03,12];
    Ahmedabad.DayMinute=[00,59];
    Ahmedabad.NightHour=[17,00];
    Ahmedabad.NightMinute=[00,59];
    citiesInfoArr{2} = Ahmedabad;
    
    Lucknow.Name='Lucknow';
    Lucknow.Center=[26.845325,80.945618];
    Lucknow.NE=[27.305423,81.478667];
    Lucknow.SW=[26.385227,80.412569];
    Lucknow.DayHour=[03,12];
    Lucknow.DayMinute=[00,59];
    Lucknow.NightHour=[17,00];
    Lucknow.NightMinute=[00,59];
    citiesInfoArr{3} = Lucknow;
    
    Bangalore.Name='Bangalore';
    Bangalore.Center=[12.971599,77.593906];
    Bangalore.NE=[13.431697,78.126955];
    Bangalore.SW=[12.511501,77.060857];
    Bangalore.DayHour=[03,12];
    Bangalore.DayMinute=[00,59];
    Bangalore.NightHour=[17,00];
    Bangalore.NightMinute=[00,59];
    citiesInfoArr{4} = Bangalore;
    
    Kolkata.Name='Kolkata';
    Kolkata.Center=[22.495881,88.337883];
    Kolkata.NE=[22.955979,88.870932];
    Kolkata.SW=[22.035783,87.804834];
    Kolkata.DayHour=[03,12];
    Kolkata.DayMinute=[00,59];
    Kolkata.NightHour=[17,00];
    Kolkata.NightMinute=[00,59];
    citiesInfoArr{5} = Kolkata;
    
    Pune.Name='Pune';
    Pune.Center=[18.510062,73.852529];
    Pune.NE=[18.97016,74.385578];
    Pune.SW=[18.049964,73.31948];
    Pune.DayHour=[03,12];
    Pune.DayMinute=[00,59];
    Pune.NightHour=[17,00];
    Pune.NightMinute=[00,59];
    citiesInfoArr{6} = Pune;
    
    citiesInfoArrTemplate = citiesInfoArr;

    % don't need the city variables if they are stored in citiesInfoArr
    clearvars Delhi Ahmedabad Luknow Bangalore Kolkata Pune;
        
    startDoyStr = sprintf('%03d', day(datetime([startDayStr,'/',yearStr],'InputFormat','MM/dd/uuuu'), 'dayofyear'));
    startDoy = str2double(startDoyStr); 
	fprintf(['startDoy = ', num2str(startDoy), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['startDoy = ', num2str(startDoy), newline]);
    fclose(fid);
    endmonth = endDayStr(1:2);
    endmonth = str2double(endmonth);
    endday = str2double(endDayStr(4:5)); 
    if endday == 1
	    endmonth = endmonth-1;
        endday = eomday(str2double(yearStr),endmonth);    
    end
    endMonthStr = sprintf('%02d',endmonth);
    endDayStr = sprintf('%02d',endday);
    endDayStr = [endMonthStr,'/',endDayStr];
    endDoyStr = sprintf('%03d', day(datetime([endDayStr,'/',yearStr],'InputFormat','MM/dd/uuuu'), 'dayofyear'));
    endDoy = str2double(endDoyStr);
	fprintf(['endDoy = ', num2str(endDoy), newline]);
	fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
	fwrite(fid, ['endDoy = ', num2str(endDoy), newline, newline]);
    fclose(fid);
	      
    year = str2double(yearStr);

    % Prevent reading from same file repeatedly
    LatTemplate = [];
    LonTemplate = [];
    LSTTemplate = [];
    QCTemplate = [];
    lastFileIO = '';

    for dateIdx = startDoy : (endDoy + 1)
        dateIdxStr = sprintf('%03d', dateIdx); 
        citiesInfoBoolArr = citiesInfoBoolArrTemplate;
        citiesInfoArr = citiesInfoArrTemplate;
        [~, mth, dd] = datevec(datenum(year,1,dateIdx));
        if mth < 10
            mthStr = ['0', num2str(mth)];
        end
        if dd < 10
            dayStr = ['0', num2str(dd)];
        end     
        targetDayStr = [mthStr,'/',dayStr];
        if ~isempty(lastFileIO)
            fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');                            
            fwrite(fid, newline);
            fclose(fid);

            fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
            fwrite(fid, ['Month/Day/Year: ', targetDayStr, '/', yearStr, newline]);
            fclose(fid);
        else
            fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
            fwrite(fid, ['Month/Day/Year: ', targetDayStr, '/', yearStr, newline]);
            fclose(fid);
        end
        
        datLocation = [datLocationTemplate, targetDayStr,'/'];         
        files_dir = dir([datLocation,'*',dateIdxStr,'*.nc']);                                   
        for b = 1:length(files_dir)        
            file = [files_dir(1).folder,'/',files_dir(b).name];             
            filename = files_dir(b).name;                     
            pos = strfind(file, '.A');
            doy_str = file(pos+2:pos + 13);            
            if strcmp(prevFile, doy_str)
                fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                fwrite(fid, [filename, ': Duplicate file detected ... skipping ... ', newline]);
                fclose(fid);
                fprintf([filename, ': Duplicate file detected ... skipping ... ', newline])
                continue;
            end               
            
            time = file(pos+10:pos+13);
            hr = time(1:2);
            hr = str2double(hr);                    
    
            mn = time(3:4);
            mn = str2double(mn);     
    
            try
                info = h5info(file);
            catch ME
                fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                fwrite(fid, ['Unable to read ', filename, ': ', cityName, 'skipping ', doy_str, newline]);
                fclose(fid);
                fprintf(['Unable to read ', filename, ': ', cityName, 'skipping ', doy_str, newline]);
                continue;
            end
            attribute_num = info.Attributes;           
            sz = size(attribute_num);
            attributes = zeros(1,5);
            % Populate 
            for i = 1 : sz(1)
                attribute_name = attribute_num(i).Name;
                if strcmp(attribute_name, 'NorthBoundingCoordinate')
                    maxLat = info.Attributes(i).Value;
                    attributes(1) = true;
                elseif strcmp(attribute_name, 'SouthBoundingCoordinate')
                    minLat = info.Attributes(i).Value;
                    attributes(2) = true;
                elseif strcmp(attribute_name, 'WestBoundingCoordinate')
                    minLon = info.Attributes(i).Value;
                    attributes(3) = true;
                elseif strcmp(attribute_name, 'EastBoundingCoordinate')
                    maxLon = info.Attributes(i).Value;  
                    attributes(4) = true;
                elseif strcmp(attribute_name, 'startDirection')
                    startDirection = info.Attributes(i).Value;  
                    attributes(5) = true;
                end                         
            end
            if ~all(attributes)
                fprintf([filename, ': missing an attribute ... skipping ', doy_str, newline]);
                continue;
            end
            prevFile = doy_str;
	        citiesInfoBoolArr = citiesInfoBoolArrTemplate;
            citiesInfoArr = citiesInfoArrTemplate;
            fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
            fwrite(fid, [newline]);            %#ok<*NBRAK2> 
            fclose(fid);
            fprintf([newline]);
            % Early For Loop to filter out cities that are out of range
            for cityIdx = 1 : size(citiesInfoArr,2)
                city = citiesInfoArr{cityIdx};
                cityName = city.Name;                
                cityNE = city.NE;
                citySW = city.SW;  
                if citySW(1) < minLat || cityNE(1) > maxLat 
                    fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                    fwrite(fid, [filename,': ',cityName, ' crop Lat out of range ', newline]);
                    fclose(fid);
                    fprintf([filename,': ',cityName, ' crop Lat out of range ', newline]);
		            citiesInfoBoolArr(cityIdx) = false;
                    continue;
                elseif citySW(2) < minLon || cityNE(2) > maxLon
                    fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                    fwrite(fid, [filename,': ',cityName, ' crop Lon is out of range ', newline]);
                    fclose(fid);
                    fprintf([filename,': ',cityName, ' crop Lon is out of range ', newline]);
		            citiesInfoBoolArr(cityIdx) = false;
                    continue;
                end
                citiesInfoBoolArr(cityIdx) = true;                
            end            	        
            citiesInfoArr(~citiesInfoBoolArr) = [];
            % only iterate if a city was found within the granule
            if any(citiesInfoBoolArr)
                for cityIdx = 1 : size(citiesInfoArr,2) 
                    city = citiesInfoArr{cityIdx};
                    cityName = city.Name;
					fprintf([newline, 'cityName = ', cityName, newline]);
                    cityCenter = city.Center;
                    cityNE = city.NE;
                    citySW = city.SW;
                    cityDayHour = city.DayHour;
                    cityDayMinute = city.DayMinute;
                    cityNightHour = city.NightHour;
                    cityNightMinute = city.NightMinute;                    
                    curName = files_dir(b).name;
					fprintf(['curName = ', curName, newline]);
                    if strcmp(curName,lastFileIO)                        
                        Lat = LatTemplate;
                        Lon = LonTemplate;
                        LST = LSTTemplate;
                        QC = QCTemplate;
                    else
                        if ~isempty(lastFileIO)
                            fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');                            
                            fwrite(fid, newline);
                            fclose(fid);
                        end
                        % Read file data
                        Lat = hdf5read(file, '/VIIRS_I5_LST/Geolocation Fields/Latitude', 'V71Dimensions', true);Lat = double(Lat); %#ok<*HDFR>
                        Lon = hdf5read(file, '/VIIRS_I5_LST/Geolocation Fields/Longitude', 'V71Dimensions', true);Lon = double(Lon);            
                        LST = hdf5read(file, '/VIIRS_I5_LST/Data Fields/LST', 'V71Dimensions', true);
                        LST = double(LST) .* 0.02;        
                        flip = false;                         
                        if strcmp(startDirection, 'Descending')
                            flip = false;
                        elseif strcmp(startDirection, 'Ascending')
                            flip = true;
                        end                            
                        QC = hdf5read(file, '/VIIRS_I5_LST/Data Fields/QC', 'V71Dimensions', true);     
                        % if flip is set rotate data
                        if flip
                            Lat = rot90(Lat,2);
                            Lon = rot90(Lon,2);
                            LST = rot90(LST,2);
                            QC = rot90(QC,2);
                        end
                        lastFileIO = curName; 
                        LatTemplate = Lat;
                        LonTemplate = Lon;
                        LSTTemplate = LST;
                        QCTemplate = QC;
                    end            
                    
                    test1 = bitget(QC,1);
                    test2 = bitget(QC,2);
                    pixtrim = (test1==1 & test2==1);       
            
                    fill = LST==0;
                    if use_Celsius
                        LST = LST - 273;
                    elseif use_Fahrenheit
                        LST = (LST-273)*9/5 + 32; 
                    elseif use_Kelvin
                        % code defaults to Kelvin
                    end
                    LST(fill) = 0;                   
                    
                    minlat = citySW(1); maxlat = cityNE(1);
                    minlon = min(citySW(2), cityNE(2));
                    maxlon = max(citySW(2), cityNE(2));
                            
                    VIIRS_sz = size(Lat);        
                    if geoLocateMinimumDistance == 1
                        [dist,~] = distance(Lat,Lon,maxlat,maxlon,earthRadiusInMeters,'degrees');                    
                        [ch,IND] = min(dist(:));
					    [I,J] = ind2sub(size(Lon),IND);    				    
                    elseif geoLocateMinimumDistance == 2                         
                        [I, J, ch] = haversin2dFindMin(Lat, Lon, maxlat, maxlon, earthRadiusInMeters);                           
                    else                                           
                        [I, J, ch, ~] = vicentyInv2dFindMinInMeters(Lat, Lon, maxlat, maxlon, maxIter); 
                    end
                    if I == 1 || J == 1 || I == VIIRS_sz(1) || J == VIIRS_sz(2)
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': NE VIIRS PT for ',cityName, ' is outside of or on scene edge', newline]);
                        fclose(fid);
                        fprintf([filename, ': NE VIIRS PT for ',cityName, ' is outside of or on scene edge', newline]);
                        continue;
                    end                    
                    if ch > 375
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': closest NE VIIRS PT for ',cityName,' is ', num2str(ch),' meters from actual point', newline]);
                        fclose(fid);
                        fprintf([filename, ': closest NE VIIRS PT for ',cityName,' is ', num2str(ch),' meters from actual point', newline]);
                        continue;
                    end     
                    if geoLocateMinimumDistance == 1
                        [dist,~] = distance(Lat,Lon,minlat,minlon,earthRadiusInMeters,'degrees'); 
                        [ch,IND] = min(dist(:));
                        [I2,J2] = ind2sub(size(Lon),IND);	
                    elseif geoLocateMinimumDistance == 2                         
                        [I2, J2, ch] = haversin2dFindMin(Lat, Lon, minlat, minlon, earthRadiusInMeters);                           
                    else                        
                        [I2, J2, ch] = vicentyInv2dFindMinInMeters(Lat, Lon, minlat, minlon, maxIter); 
                    end
                    if I2 == 1 || J2 == 1 || I2 == VIIRS_sz(1) || J2 == VIIRS_sz(2)
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': SW VIIRS PT for ',cityName, ' is outside of or on scene edge', newline]);
                        fclose(fid);
                        fprintf([filename, ': SW VIIRS PT for ',cityName, ' is outside of or on scene edge', newline]);
                        continue;
                    end        
                    if ch > 375
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': closest SW VIIRS PT for ',cityName,' is ', num2str(ch),' meters from actual point', newline]);
                        fclose(fid);
                        fprintf([filename, ': closest SW VIIRS PT for ',cityName,' is ', num2str(ch),' meters from actual point', newline]);
                        continue;
                    end
                    
                    % 375m grid
                    spacepix = 0.004;        
                    lats = minlat:spacepix:maxlat;
                    lons = minlon:spacepix:maxlon;
	                
	                szlats = size(lats);
	                szlons = size(lons);
	                
	                if szlats(1) <= 0 || szlats(2) <= 0
		                fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': ', cityName, ' lats dimensions are incorrect for ndgrid.', newline]);
                        fclose(fid);
                        fprintf(['lats dimensions are incorrect for meshgrat.', newline]);
		                continue;
	                end
                    if szlons(1) <= 0 || szlons(2) <= 0 
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': ', cityName, ' lons dimensions are incorrect for ndgrid.', newline]);
                        fclose(fid);
                        fprintf(['lons dimensions are incorrect for meshgrat.', newline]);
                        continue;
                    end
                    
                    [latg,long] = ndgrid(lats,lons);
                    latg = flipud(latg);
                    sz = size(latg);
            
                    first_cut = find(Lat(:,:)< min(latg(:)) |  Lat(:,:)> max(latg(:)));
                    Lat(first_cut) = [];
                    Lon(first_cut) = [];
                    
                    if isempty(Lat)
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': ', cityName, ' Lat is empty after first crop.', newline]);
                        fclose(fid);
                        fprintf([filename, ': ', cityName, ' Lat is empty after first crop.', newline]);
                        continue;
                    end
                    if isempty(Lon)
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': ', cityName, ' Lon is empty after first crop.', newline]);
                        fclose(fid);
                        fprintf([filename, ': ', cityName, ' Lon is empty after first crop.', newline]);
                        continue;
                    end
            
                    second_cut = find(Lon(:,:)< min(long(:)) |  Lon(:,:)> max(long(:)));
                    Lon(second_cut) = [];
                    Lat(second_cut) = [];
                    
                    if isempty(Lat)
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': ', cityName, ' Lat is empty after second crop.', newline]);
                        fclose(fid);
                        fprintf([filename, ': ', cityName, ' Lat is empty after second crop.', newline]);
                        continue;
                    end
                    if isempty(Lon)
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': ', cityName, ' Lon is empty after second crop.', newline]);
                        fclose(fid);
                        fprintf([filename,': ', cityName,' Lon is empty after second crop.', newline]);
                        continue;
                    end
                    QC(first_cut) = [];
                    QC(second_cut) = [];        
                    LST(first_cut) = [];
                    LST(second_cut) = [];
                    pixtrim(first_cut) = [];
                    pixtrim(second_cut) = [];
            
                    Lat(pixtrim) = [];
                    Lon(pixtrim) = [];
                    LST(pixtrim) = [];
                    QC(pixtrim) = [];
                    
                    test1 = bitget(QC,1);
                    test2 = bitget(QC,2);
                    ocean_mask = (test1==1 & test2==1); 
                    cloud_mask = (test1 == 0 & test2 == 1);
                    
                    LST(ocean_mask) = NaN;
                    LST(cloud_mask) = 0;
                    
                    if isempty(LST)
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': for ', cityName, ' LST is empty', newline]);
                        fclose(fid);
						fprintf([filename, ': for ', cityName,' has fill values ', num2str(percent_zero) ,'%% before griddata', newline]);
                        continue;
                    end
                    
                    total_sz = size(LST);        
                    total_pixels = total_sz(1) * total_sz(2);
                    szero = (LST(:) == 0);
                    szero = sum(szero(:),'omitnan');
                    percent_zero = (szero/total_pixels) * 100; 
                    if percent_zero == 100
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': for ', cityName,' has fill values ', num2str(percent_zero) ,'% before griddata', newline]);
                        fclose(fid);
						fprintf([filename, ': for ', cityName,' has fill values ', num2str(percent_zero) ,'%% before griddata', newline]);
                        continue;
                    end        
                    
                    [lat_grid,lon_grid,LST_grid] = griddata(Lat,Lon,LST,latg,long,'natural');         
                    [~,~,cloud_mask_grid] = griddata(Lat,Lon,double(cloud_mask),latg,long,'natural'); 
                    [~,~,QC_grid] = griddata(Lat,Lon,double(QC),latg,long,'natural'); 
                    QC_grid = uint8(QC_grid);
                    cloud_mask_grid(isnan(cloud_mask_grid(:))) = 0;
                    cloud_mask_grid = logical(cloud_mask_grid);                    
                    LST_grid(cloud_mask_grid) = 0;            
					
                    total_sz = size(LST_grid);        
                    total_pixels = total_sz(1) * total_sz(2);
                    szero = (LST_grid(:) == 0);
                    szero = sum(szero(:),'omitnan');
                    percent_zero = (szero/total_pixels) * 100;         
					
                    if percent_zero < percentZeroThreshold
						fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': percent_zero = ', num2str(percent_zero), '%', newline]);
                        fclose(fid);
						fprintf([filename, ': percent_zero = ', num2str(percent_zero), '%%', newline]);
					else
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': for ', cityName,' has fill values ', num2str(percent_zero) ,' > ',num2str(percentZeroThreshold),'% after griddata', newline]);
                        fclose(fid);
                        fprintf([filename, ': for ', cityName,' has fill values ', num2str(percent_zero) ,' > ',num2str(percentZeroThreshold),'%% after griddata', newline]);
                        continue;
                    end
                    if percent_zero == 100
                        fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
                        fwrite(fid, [filename, ': for ', cityName,' is all fill values', newline]);
                        fclose(fid);
                        fprintf([filename, ': for ', cityName,' is all fill values', newline]);
                        continue;
                    end
                    out_rows = all(isnan(LST_grid),2);    % rows that are all nan
                    out_cols = all(isnan(LST_grid),1);   % cols that are all nan
                    LST_grid(out_rows, :) = [];
                    LST_grid(:, out_cols) = [];
                    lat_grid(out_rows, :) = [];
                    lat_grid(:, out_cols) = [];
                    lon_grid(out_rows, :) = [];
                    lon_grid(:, out_cols) = [];
                    % Output Geotiff
                    latlim = [min(lat_grid(:)) max(lat_grid(:))];
                    lonlim = [min(lon_grid(:)) max(lon_grid(:))];
                    rasterSize = size(lat_grid);
                    R = georefcells(latlim,lonlim,rasterSize,'ColumnsStartFrom','north');        
                    
                    fileswrite = [file(pos-12:end-3),'_LST.tif'];                      
                    fileswrite = [outputLoc, cityName,'/', yearStr, '/', fileswrite];
                    imgwrite = [outputLoc, cityName,'/', yearStr, '/images/'];
                    if ~displayImages
                        set(0,'DefaultFigureVisible','off');
                    end
                    if writeImages
                        figure;
                        imagesc([LST_grid])                        
                        hold on;   
                        uniqueValues = unique(LST_grid(:));
                        uniqueValues(uniqueValues == 0) = [];
                        uniqueValues(isnan(uniqueValues(:))) = [];
                        if length(uniqueValues(:)) > 1                                                                                       
                            clim([uniqueValues(1) uniqueValues(end)]);
                        else
                            error([filename, ' has only one unique value'])
                        end
                        if use_Celsius                
                            clim([10 prctile(LST_grid(:),99.9)]);
                        elseif use_Fahrenheit
                            clim([35 max(LST_grid(:))]); 
                        elseif use_Kelvin
                            
                        end
                        colorbar
                        titlenameColPos = strfind(file,'.002') - 1;
                        title([cityName,' ', file(pos+2:titlenameColPos(end))],'FontSize',18); 
                        set(gca,'Xticklabel',[],'YTicklabel',[],'Xtick',[],'YTick',[]); 
                        set(gca,'YGrid','on','Fontsize',14,'XGrid','on');
                        set(gcf,'Color','White');
                        exportgraphics(gca,[imgwrite,file(pos-12:titlenameColPos(end)),'_LST.jpg']); 
                        if ~displayImages
                            set(0,'DefaultFigureVisible','on');        
                        end
	                    close(gcf); clc;        		
                    end
                    if writeTifs
                        geotiffwrite(fileswrite,LST_grid,R);               
                    end
                    cnt = cnt + 1;
                    fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
					if writeImages && writeTifs
						fprintf([filename, ': ',cityName,' good; tif and jpg files generated', newline]);
						fwrite(fid, [filename, ': ',cityName,' good; tif and jpg files generated', newline]);
					elseif writeImages && ~writeTifs
						fprintf([filename, ': ',cityName,' good; jpg file generated', newline]);
						fwrite(fid, [filename, ': ',cityName,' good; jpg file generated', newline]);
					elseif ~writeImages && writeTifs
						fprintf([filename, ': ',cityName,' good; tif file generated', newline]);
						fwrite(fid, [filename, ': ',cityName,' good; tif file generated', newline]);
                    else
                        fprintf([filename, ': ',cityName,' good', newline]);
						fwrite(fid, [filename, ': ',cityName,' good', newline]);
                    end                    
                    fclose(fid);                    
                end
            else
                fprintf(newline);
                continue;
            end
        end
    end	
    fprintf(['\n', num2str(cnt), ' tif files created.', newline])
    fid = fopen(['log_', startDayLogStr, '_through_',endDayLogStr ,'.txt'], 'a');
    fwrite(fid, [num2str(cnt), ' scenes identified.']);
    fclose(fid);
end