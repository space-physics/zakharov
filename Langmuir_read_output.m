function data = Langmuir_read_output(fn)

fid = fopen(fn,'r');
data = fread(fid,'float64');
fclose(fid);

if nargout==0
clear data
end
end