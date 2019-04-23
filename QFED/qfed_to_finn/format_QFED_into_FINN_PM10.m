% USE WITH MATLAB R2010b OR NEWER FOR NETCDF4 SUPPORT
% This is a modified version of a script given to me by Xinxin Ye in Pablo
% Saide's group to convert QFED data into FINN format for use in the
% wrffire files. 
% 
% The files necessary for this script to run are:
% deg2utm_vec_nozone.m
% read_netcdf_vars_1hour.m
% redhc_CO_gas_factors_QFED.csv

close all
clear

date_i=[2013 08 22 0 0 0];
date_f=[2013 08 23 0 0 0];

path_qfed='/met1/WRF/DATA/QFED/2013/M08/';
spc_read_qfed={'bc','ch4','co','co2','nh3','no','oc','pm25','so2'};
mol_weight=[0 16.04 28.01 44.01 17.03 30.01 0 0 64.07]; %g/mol
prefix_qfed='qfed2.emis_';
sufix1_qfed='.006.';
sufix2_qfed='.nc4';
var_extract={'biomass_xf','biomass_tf','biomass_sv','biomass_gl'}; %same order as in FINN

dt_files=24;%in hours
file_out_prefix='../QFED_in_FINN_format_pm10_';
%file_out_prefix='QFED_in_FINN_format_MOZART_pm10_';
file_out_sufix='.txt';
%geos5 0.1 center degree grid
delta_latlon=0.1; %in degrees
lat_i=-89.95;
lat_f= 89.95;
lon_i=-179.95;
lon_f= 179.95;
lat_vec=(lat_i:delta_latlon:lat_f);
lon_vec=(lon_i:delta_latlon:lon_f);
[lon_grid lat_grid]=ndgrid(lon_vec,lat_vec);
[Nx,Ny]=size(lon_grid);

%Read SAPRC table
%file=fopen('saprc_CO_gas_factors_QFED.csv','r');
%Read MOZART
%file=fopen('mozaic_CO_gas_factors_QFED.csv','r');
%Read REDHC
file=fopen('redhc_CO_gas_factors_QFED.csv','r');
data=textscan(file,'%s%f%f%f%f','HeaderLines',1,'Delimiter',',');
fclose(file);
saprc_spc=data{1};
saprc_co_to_voc=[data{2} data{3} data{4} data{5}]; %same order as in "var_extract"

% Compute area (m2) per each gridcell
x_pol=NaN*zeros(4,Nx,Ny);
y_pol=NaN*zeros(4,Nx,Ny);
[aux_x aux_y]=deg2utm_vec_nozone(reshape(lat_grid,numel(lat_grid),1)-0.05, ...
                                 reshape(lon_grid,numel(lon_grid),1)-0.05);
x_pol(1,:,:)=reshape(aux_x,Nx,Ny);
y_pol(1,:,:)=reshape(aux_y,Nx,Ny);
[aux_x aux_y]=deg2utm_vec_nozone(reshape(lat_grid,numel(lat_grid),1)-0.05, ...
                                 reshape(lon_grid,numel(lon_grid),1)+0.05);
x_pol(2,:,:)=reshape(aux_x,Nx,Ny);
y_pol(2,:,:)=reshape(aux_y,Nx,Ny);
[aux_x aux_y]=deg2utm_vec_nozone(reshape(lat_grid,numel(lat_grid),1)+0.05, ...
                                 reshape(lon_grid,numel(lon_grid),1)+0.05);
x_pol(3,:,:)=reshape(aux_x,Nx,Ny);
y_pol(3,:,:)=reshape(aux_y,Nx,Ny);
[aux_x aux_y]=deg2utm_vec_nozone(reshape(lat_grid,numel(lat_grid),1)+0.05, ...
                                 reshape(lon_grid,numel(lon_grid),1)-0.05);
x_pol(4,:,:)=reshape(aux_x,Nx,Ny);
y_pol(4,:,:)=reshape(aux_y,Nx,Ny);
area_geos5 = squeeze(polyarea(x_pol,y_pol));
%fix problems
%fix horizontal maximums
area_geos5(:,899)=area_geos5(:,902);
area_geos5(:,900)=area_geos5(:,902);
area_geos5(:,901)=area_geos5(:,902);
%fix vertical maximums
[aux_x aux_y]=find(area_geos5>2e+8);
x_index_fix=union(aux_x,aux_x);
area_geos5(x_index_fix,:)=area_geos5(x_index_fix-2,:);
clear x_pol y_pol aux_x aux_y;

%build format
format='%6d,%6d,%6d,%20.10E,%20.10E,'; %first 3 are integers, then lat and lon
for i=1:numel(spc_read_qfed)
 format=[format '%20.10E,'];
end
for i=1:numel(saprc_spc)
 format=[format '%20.10E,'];
end
format=[format '%20.10E,']; %pm10
format=[format '\n'];

%header
header_out='DAY,TIME,GENVEG,LATI,LONGI';
for i=1:numel(spc_read_qfed)
 header_out=[header_out ',' upper(spc_read_qfed{i})];
end
for i=1:numel(saprc_spc)
 header_out=[header_out ',' upper(saprc_spc{i})];
end
header_out=[header_out ',' upper('PM10')];

%open file
file=fopen([file_out_prefix datestr(date_i,'yyyymmdd') ...
     '_' datestr(date_f,'yyyymmdd') file_out_sufix],'w');
fprintf(file,'%s\n',header_out);

Ndays_extract=ceil(datenum(date_f)-datenum(date_i));
datenum_aux=datenum(date_i);
for j=1:Ndays_extract
 data_qfed=NaN*zeros(Nx,Ny,4,numel(spc_read_qfed));
 datenum_list(j)=datenum_aux;
 datestr_aux=[datestr(datenum_aux,'yyyymmdd')]
 index_t=1;
 for i=1:numel(spc_read_qfed)
  [names_spc data] = read_netcdf_vars_1hour([path_qfed ...
    prefix_qfed spc_read_qfed{i} sufix1_qfed datestr_aux sufix2_qfed], ...
    var_extract,index_t);
  for k=1:4
   data_qfed(:,:,k,i)=data{k};
  end
 end
 
 % Change units
 % QFED units: kg s-1 m-2;
 % FINN units: gases: mole/day, aerosols= kg/day
 for i=1:numel(spc_read_qfed)
  if(mol_weight(i)>0) %gas
   for k=1:4 
    data_qfed(:,:,k,i)=(1000.0*3600*24/(mol_weight(i)))*squeeze(data_qfed(:,:,k,i)).*area_geos5;
   end
  else %aerosol
   for k=1:4
    data_qfed(:,:,k,i)=(3600*24)*squeeze(data_qfed(:,:,k,i)).*area_geos5;
   end  
  end
 end

 % Generate VOC
 co_index=find(strcmp('co',spc_read_qfed));
 data_saprc_voc=NaN*zeros(Nx,Ny,4,numel(saprc_spc));
 for i=1:numel(saprc_spc)
  for k=1:4
   data_saprc_voc(:,:,k,i)=data_qfed(:,:,k,co_index).*saprc_co_to_voc(i,k);
  end
 end

 %Add PM10
 pm10_to_pm25_ratio=2.0; % based on Akagi et al 2011 emission factors
 pm25_index=find(strcmp('pm25',spc_read_qfed));
 data_pm10=NaN*zeros(Nx,Ny,4);
 for k=1:4
  data_pm10(:,:,k)=data_qfed(:,:,k,pm25_index).*pm10_to_pm25_ratio; %PM10 in kg/day
 end

 % Write on FINN format 
 for k=1:4
  index_fire=(data_qfed(:,:,k,1)>0);
  date_aux=datevec(datenum_aux);
  if(sum(sum(index_fire))>0)
   Ndata=numel(find(index_fire));
   DAY=(datenum_aux-datenum([date_aux(1) 01 0 0 0 0]))*ones(Ndata,1);
   TIME=0*ones(Ndata,1); %all 0 UTC
   GENVEG=k*ones(Ndata,1);
   LONI=lon_grid(index_fire);
   LATI=lat_grid(index_fire);
   emiss_write=NaN*zeros(Ndata,numel(spc_read_qfed));
   for i=1:numel(spc_read_qfed) 
    aux=squeeze(data_qfed(:,:,k,i));
    emiss_write(:,i)=aux(index_fire);
   end
   emiss_write_saprc=NaN*zeros(Ndata,numel(saprc_spc));
   for i=1:numel(saprc_spc)
    aux=squeeze(data_saprc_voc(:,:,k,i));
    emiss_write_saprc(:,i)=aux(index_fire);
   end
   %pm10
   emiss_write_pm10=NaN*zeros(Ndata,1);
   aux=squeeze(data_pm10(:,:,k));
   emiss_write_pm10=aux(index_fire);

   matrix_write=[DAY TIME GENVEG LATI LONI emiss_write emiss_write_saprc emiss_write_pm10];
   str_data=sprintf(format,matrix_write');
   str_data=strrep(str_data,'E+','D+');
   str_data=strrep(str_data,'E-','D-');
   fprintf(file,'%s',str_data);
  end
 end 
  
 datenum_aux=datenum(date_i)+j*dt_files/(24); %increase time step 

end

fclose(file);

