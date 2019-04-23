% Coded by Pablo Saide, pablo-saide@uiowa.edu

function [varname data] = read_netcdf_vars_1hour(name,vars_extract,index_t) 

file = netcdf.open(name, 'NC_NOWRITE');
[ndims_file,nvars_file,ngatts_file,unlimdimid_file] = netcdf.inq(file);

%vars and local atributes
for i=0:(nvars_file-1)
 [varname{i+1},xtype(i+1),dimids{i+1},natts(i+1)] = netcdf.inqVar(file,i);
% data_aux{i+1}= netcdf.getVar(file,i);
% attributes:
% for j=0:(natts(i+1)-1)
%  attname{i+1,j+1}=netcdf.inqAttName(file,i,j);
%  [xtype_att(i+1,j+1),attlen(i+1,j+1)] = netcdf.inqAtt(file,i,attname{i+1,j+1});
%  attrvalue{i+1,j+1} = netcdf.getAtt(file,i,attname{i+1,j+1});
% end
end

%dimensions
for i=0:(ndims_file-1)
 [dimname{i+1}, dimlen(i+1)] = netcdf.inqDim(file,i);
end

if(strcmp(lower(vars_extract{1}),'all')) %extrac all
 data=cell(nvars_file,1);
 for i=1:nvars_file
  data{i}= netcdf.getVar(file,i-1);
 end

else %extract selection
 indx_extract=zeros(numel(vars_extract),1);
 data=cell(numel(vars_extract),1);
 for i=1:numel(vars_extract)
%  vars_extract{i}
  logical=strcmp(vars_extract{i},varname);
  indx_extract(i)=find(logical);
  
  dimids_aux=dimids{indx_extract(i)};
%  varname{indx_extract(i)}
%  dimids{indx_extract(i)}
  if(numel(dimids{indx_extract(i)})==4)
   aux=dimlen(dimids_aux+1);
   Nx=aux(1);
   Ny=aux(2);
   Nz=aux(3);
   Nt=aux(4);
   start=[0 , 0 , 0 , index_t-1];
   count=[Nx , Ny , Nz , 1];
  elseif(numel(dimids{indx_extract(i)})==3)
   aux=dimlen(dimids_aux+1);
   Nx=aux(1);
   Ny=aux(2);
   Nt=aux(3);
   start=[0 , 0 , index_t-1];
   count=[Nx , Ny , 1];
  elseif(numel(dimids{indx_extract(i)})==2)
   aux=dimlen(dimids_aux+1);
   Nx=aux(1);
   Ny=aux(2);
   start=[0 , 0];
   count=[Nx , Ny];
  else
   fprintf(1,'%s\n','Dimension unknown');
   return
  end

  data{i}= netcdf.getVar(file,indx_extract(i)-1,start,count);

 end
 varname=varname(indx_extract);

end

netcdf.close(file)
