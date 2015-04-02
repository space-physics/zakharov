function params = Langmuir_read_params(fnParam)

fidparam = fopen(fnParam,'r');
params = fread(fidparam,'float64');
fclose(fidparam);

%if nargout==0
    display([fnParam,': vbeam_ev=',num2str(params(15)),' vbeam=',num2str(params(16)),' tetabeam=',num2str(params(17))])
%clear params
%end
end