function P = Langmuir_read_params(fn)

assert(exist(fn,'file')==2,[fn,' not found'])

f = fopen(fn,'r');

params = fread(f,'float64');
fclose(f);

P.pi=params(1);
P.me=params(2);
P.electroncharge=params(3);
P.mi=params(4);
P.Kb=params(5);
P.eV=params(6);
P.epsilon0=params(7);
P.Z=params(8);
P.Te=params(9);
P.Ti=params(10);
P.nuic=params(11);
P.nuec=params(12);
P.n0=params(13);
P.nbeam=params(14);
P.vbeam_ev=params(15);
P.vbeam=params(16);
P.tetabeam=params(17);
P.endTime=params(18);
P.Tstep=params(19);
P.TT=params(20);
P.res=params(21);
P.TT_res=params(22);
P.L=params(23);
P.N=params(24);
P.Xstep=params(25);
P.QW=params(26);
P.eta=params(27);
vve=params(28);
P.Cs=params(29);
P.omegae=params(30);
P.lambdaD=params(31);
P.SEED=params(32);

disp(P)


end