function [nnComplex,EEComplex] = Langmuir_read_data(fnNN,fnEE,useSinglePrec)
%% nn
display(['reading file ',fnNN])
fid=fopen(fnNN);
    if useSinglePrec
        nn=fread(fid, 'float64=>float32');
    else %double
        nn=fread(fid,'float64');
    end
    
    fclose(fid);
%% EE    
display(['reading file ',fnEE])
    fid2=fopen(fnEE);
    if useSinglePrec
        EE=fread(fid2,'float64=>float32');
    else %double
    EE=fread(fid2,'float64');
    end

    fclose(fid2);
    
%% give only needed format
% creates complex number
% note: odd indices are "real", even indicies are "imag" -- odd,even pair

EEComplex(1,:) = EE(1:2:end)+1j*EE(2:2:end);
nnComplex(1,:) = nn(1:2:end)+1j*nn(2:2:end);

end
