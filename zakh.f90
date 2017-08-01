! 1D Langmuir turbulence simulation according to Guio and Forme 2006
! Last update: 9/5/2013
! For Kappa on: 01/05/2016

! Michael Hirsch Oct 2013 -- updated vbeam,tetabeam to use C++ vector format

program zakharov1d

use, intrinsic:: iso_fortran_env, only: wp=>real64

implicit none

real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

real(wp), parameter :: me=9.10938356e-31_wp ! kg
real(wp), parameter :: electroncharge=1.60217662e-19_wp ! coulombs
real(wp), parameter :: mi=16*1.66e-27_wp ! atomic oxygen
real(wp), parameter :: Kb=1.38064852e-23_wp    ! Boltzmann cte
real(wp), parameter :: eV=1.602176565e-19_wp, epsilon0=8.854187817e-12_wp

integer, parameter :: Z=1
real(wp), parameter :: Te=3000.0_wp, Ti=1000.0_wp
real(wp), parameter :: nuic=1.0_wp ! ion collision freq
real(wp), parameter :: nuec=100.0_wp ! electron collision freq
real(wp), parameter :: n0=5.0e11_wp ! background density

real(wp), parameter :: vbeam_ev(*) = [500.0_wp]
integer, parameter :: Nvbeam=size(vbeam_ev)
real(wp) :: vbeam(Nvbeam), tetabeam(Nvbeam)


real(wp), parameter :: power_n_cte=7.7735e6_wp/sqrt(53.5_wp)
real(wp), parameter :: power_E_cte=0.5033_wp*0.5033_wp/sqrt(3.0_wp)

!    Kappa parameters
! look at Broughton et al., Modeling of MF wave mode conversion
real(wp), parameter :: se_percent=0.001_wp
real(wp), parameter :: kappa=1.584_wp
real(wp), parameter :: T_se=18.2_wp*eV
real(wp), parameter :: theta_se=sqrt((kappa-1.5_wp)/kappa*2.0_wp*T_se/me)
real(wp), parameter :: se_cte= 0.7397_wp / theta_se**3.0_wp

!    Simulation parameters

real(wp), parameter :: endTime=100.0e-3_wp ! simulation ends (seconds)
real(wp), parameter :: Tstep=0.5e-7_wp ! simulation time steps
integer, parameter :: TT= floor(endTime / Tstep)              ! floor(endTime/Tstep)+2;
integer, parameter :: res=20
integer, parameter :: TT_res = floor(endTime/Tstep/res)  !floor(endTime/Tstep/res);
real(wp), parameter :: L=70.0_wp           ! simulation box length (meter)j
integer, parameter :: N=2046          ! number of samples in L; should be devidable by 6
real(wp), parameter :: Xstep= L / N
integer, parameter :: QW=1     ! number of realizations

integer, allocatable :: SEED(:)
integer :: Nseed, clock
real(wp) :: rdist


real(wp), parameter :: eta=(Te+3.0_wp*Ti) / Te
real(wp), parameter :: ve=sqrt(Kb*Te/me)
real(wp), parameter :: Cs=sqrt(eta*me/mi)*ve
real(wp), parameter :: omegae=sqrt(n0*electroncharge**2.0_wp/me/epsilon0)
real(wp), parameter :: lambdaD=ve/omegae

real(wp), allocatable :: beamev(:), nbeam(:)
integer :: Nnbeam

character(:), allocatable :: odir, ofn
character(256) :: argv
integer :: argc,i, ii,iij1,iij2,counter1, q,tt1,c1,c2,beami, beamj, realization,u,uEE,uNN
real(wp) :: tic,toc

!---- main loop variables

integer :: p(N)
real(wp) :: k(N), Xsection_ion=0.0_wp, Xsection_pl=0.0_wp, E_thermal_k_squared, n_thermal_k_squared, Source_factor_E(N), &
  Source_factor_n(N), omegaL(N), gamas, nui(N), gamal1, gamal2, gamal3, gamal, nue(N), output1(N,12)

type params

  real(wp) :: pi, me,electroncharge,mi, Kb,eV,epsilon0, Te, Ti, nuic, nuec, n0, nbeam, &
            vbeam_ev, vbeam, tetabeam, endTime, Tstep, TT, res, TT_res, L, N, Xstep, QW, &
            eta, ve, Cs, omegae, lambdaD

  integer :: Z
!  	//have to include other parameters regarding the Kappa distribution

end type params

type(params) :: parameters


real(wp):: EE(3,N,2), nn(3,N,2), vv(3,N,2), CC(2), SSE(N,2), SSn(2), cte1, cte2, &
    k1(N,2), k2(N,2), k3(N,2), k4(N,2), kn1(2), kn2(2), kn3(2), kn4(2), kv1(2), kv2(2), kv3(2), kv4(2), &
    total_EE(20000*N*2),total_nn(20000*N*2)

integer :: LL,UU,pp


!----

! argparse
argc = command_argument_count()

call get_command_argument(1,argv)
odir = trim(argv)
print *,'writing output to', odir
call execute_command_line('mkdir -p '//odir)

allocate(beamev(argc-1))
do i = 2,argc
   call get_command_argument(i,argv)
    read(argv,*) beamev(i-1)
enddo

Nnbeam = size(beamev)

nbeam = beamev*n0

print *, "Nnbeam=", Nnbeam
print *, "Nvbeam=",Nvbeam
print *,"TT=",TT,"time steps"

call cpu_time(tic)

! initialization

vbeam=sqrt(eV*vbeam_ev*2/me)
 tetabeam=0.3*vbeam

print *,"vbeam=",vbeam
print *,"tetabeam=",tetabeam

call random_seed(size=nseed)
allocate(seed(nseed))

do beami=1,Nnbeam
do beamj=1,Nvbeam

  parameters%pi = pi
  parameters%me = me
  parameters%electroncharge = electroncharge
  parameters%mi = mi
  parameters%Z = Z
  parameters%Kb = Kb
  parameters%eV = eV
  parameters%epsilon0= epsilon0
  parameters%Te = Te
  parameters%Ti = Ti
  parameters%nuic = nuic
  parameters%nuec= nuec
  parameters%n0 = n0
  parameters%nbeam = nbeam(beami)
  parameters%vbeam_ev = vbeam_ev(beamj)
  parameters%vbeam = vbeam(beamj)
  parameters%tetabeam = tetabeam(beamj)
  parameters%endTime = endTime
  parameters%Tstep = Tstep
  parameters%TT = TT
  parameters%res = res
  parameters%TT_res = TT_res
  parameters%L = L
  parameters%N = N
  parameters%Xstep = Xstep
  parameters%QW = QW
  parameters%eta = eta
  parameters%ve = ve
  parameters%Cs = Cs
  parameters%omegae = omegae
  parameters%lambdaD = lambdaD

  write(argv,'(A,I0.3,A,I0.3)')  odir//"/parameters_n" , beami, "_v" , beamj
  ofn = trim(argv)

  open(newunit=u, file=ofn, form='formatted',status='unknown',action='write')
  write(u,*) parameters
  write(u,*) SEED
  close(u)
  print *, "Wrote",parameters,SEED,"to",ofn

  do ii=1,N
    p(ii)=ii-N/2

    if (ii==N/2) then
      k(ii) = 0
      Xsection_ion=0.0_wp
      Xsection_pl=0.0_wp
      n_thermal_k_squared=0.0_wp
      E_thermal_k_squared=0.0_wp
    else
      k(ii)=2*pi*p(ii)/N/Xstep
      call Xsection(Xsection_ion,Xsection_pl,k(ii))
      Xsection_ion=Xsection_ion/N/N
      Xsection_pl=Xsection_pl/N/N
      n_thermal_k_squared=Xsection_ion*n0
      E_thermal_k_squared=Xsection_pl *n0* electroncharge/epsilon0/k(ii)**2.0_wp
    end if

    omegaL(ii)=sqrt(omegae**2.0_wp+3*k(ii)*ve**2.0_wp)
    gamas= (-1)*sqrt(pi/8)*(sqrt(me/mi)+ Te/Ti**2.0_wp/sqrt(Te/Ti)*exp((-1)*(Te/2.0/Ti)-1.5))*abs(k(ii))*Cs
		!gamas= (-1)*sqrt(pi/2)*(sqrt(me/mi)+4*pow(Te/2/Ti,2)/sqrt(Te/2/Ti)*exp((-1)*(Te*4/Ti)))*abs(k(ii))*Cs*10;   //based on Robinson 2002
		!gamas= (-1)*sqrt(pi/8)*pow(1/(1+k(ii)*k(ii)*lambdaD*lambdaD)+3*Ti/Te,2)/sqrt(1/(1+k(ii)*k(ii)*lambdaD*lambdaD)+3*Ti/Te)*(sqrt(me/mi)+pow(Te/Ti,2)/sqrt(Te/Ti)*exp((-1)*(Te/2.0/Ti)/(1+k(ii)*k(ii)*lambdaD*lambdaD)-1.5))*abs(k(ii))*Cs;   //Based on some Chinese paper!!
    nui(ii)=(nuic/2-gamas)

    if (ii==N/2) then
      gamal=0.0_wp           ! this one is Nan due to division by zero
      gamal1=0.0_wp
      nue(ii)=nuec/2-gamal1
      Source_factor_n(ii)=0.0_wp
      Source_factor_E(ii)=0.0_wp
    else
      gamal1=(-1)*sqrt(pi/8)* omegae/k(ii)/ve**2.0_wp * sign(1.0_wp,k(ii)) * omegaL(ii)**2.0_wp / &
              (k(ii)*ve)*exp((-1)*omegaL(ii)/k(ii)/ve**2.0_wp/2)  !Landau damping due to the thermal electrons

      gamal2=(-1)*sqrt(pi/8)* omegae / k(ii) / tetabeam(beamj)**2.0_wp * sign(1.0_wp,k(ii)) * nbeam(beami) / &
              n0*omegaL(ii)*(omegaL(ii)-k(ii)*vbeam(beamj)) / (k(ii)*tetabeam(beamj)) * &
              exp((-1)* omegaL(ii)-k(ii)*vbeam(beamj) / k(ii) / tetabeam(beamj)**2.0_wp/2) !Landau damping due to the beam
!gamal2=(-1)*sqrt(pi/8)*pow(omegae/k(ii)/tetabeam.at(beamj),2)* sign(1,k(ii)) *nbeam.at(beami)/n0*omegaL(ii)*(omegaL(ii)-k(ii)*vbeam.at(beamj))/(k(ii)*tetabeam.at(beamj))*exp((-1)*pow((omegaL(ii)-k(ii)*vbeam.at(beamj))/k(ii)/tetabeam.at(beamj),2)/2);  //Landau damping due to the beam
      gamal3=(-1)*sqrt(pi)*omegae*omegaL(ii)**2.0_wp / k(ii)**3.0_wp * sign(1.0_wp,k(ii)) *se_cte * &
              (1+ omegaL(ii)**2.0_wp/kappa/ k(ii)*theta_se**2.0_wp)**((-1)*(kappa+1))

      gamal=gamal1*(1-se_percent)+gamal2+se_percent*gamal3 ! here decide to include the beam and Kappa distribution
      nue(ii) = nuec/2-gamal1

      Source_factor_n(ii)=2*nui(ii)*sqrt(4*nui(ii)*k(ii)*k(ii)/(4*nui(ii)*nui(ii)+k(ii)*k(ii))* &
                          n_thermal_k_squared*power_n_cte)
      ! source factor is the factor by which we balance the thermal source intensity
      Source_factor_E(ii)=sqrt(2*nue(ii)*E_thermal_k_squared*power_E_cte)
      nue(ii)=nuec/2-gamal
    end if


    output1(ii,1) = p(ii)
    output1(ii,2) = k(ii)
    output1(ii,3) = Xsection_ion
    output1(ii,4) = Xsection_pl
    output1(ii,5) = E_thermal_k_squared
    output1(ii,6) = n_thermal_k_squared
    output1(ii,7) = omegaL(ii)
    output1(ii,8) = gamas
    output1(ii,9) = nui(ii)
    output1(ii,10) = nue(ii)
    output1(ii,11) = Source_factor_E(ii)
    output1(ii,12) = Source_factor_n(ii)
  end do ! ii

  write(argv,'(A,I0.3,A,I0.3)') odir//"/output1_n",beami,"_v", beamj
  open(newunit=u,file=trim(argv),status='unknown',action='write')

  write(u,*) output1

  close(u)
  print *, "Wrote to", trim(argv)


do realization=1,QW
  cte2=omegae/2.0_wp/n0/1

  call system_clock(clock)
  seed = clock + 37 * [ (i - 1, i = 1, nseed) ]
  call random_seed(put=seed)

  write(argv,'(A,I0.3,A,I0.3,A,I0.3)') odir//"/EE",realization,"_n",beami,"_v",beamj
  open(newunit=uEE,file=trim(argv), status='unknown',action='write')
  print *,'writing to',trim(argv)

  write(argv,'(A,I0.3,A,I0.3,A,I0.3)') odir//"/nn",realization, "_n" ,beami, "_v",beamj
  open(newunit=uNN,file=trim(argv), status='unknown',action='write')
  print *,'writing to',trim(argv)

!   main loops

  do iij1=1,4
    do iij2=1,N
      call random_number(rdist)
      EE (iij1,iij2,1)=sqrt(output1(iij2,5)/2.0)*rdist
      call random_number(rdist)
      EE (iij1,iij2,2)=sqrt(output1(iij2,5)/2.0)*rdist
      call random_number(rdist)
      nn (iij1,iij2,1)=sqrt(output1(iij2,6)/2.0)*rdist
      call random_number(rdist)
      nn (iij1,iij2,2)=sqrt(output1(iij2,6)/2.0)*rdist

      vv(iij1,iij2,1)=0.0_wp
      vv(iij1,iij2,2)=0.0_wp

    enddo ! iij2 N

    do iij2=1,N/2
      nn (iij1,N-iij2,1)=  nn(iij1,iij2,1)
      nn (iij1,N-iij2,2)= -nn(iij1,iij2,2)
    enddo ! iij2 N/2

    nn (iij1,N/2,1)=0
    nn (iij1,N/2,2)=0

  end do ! iij1 4

  counter1=0
  do tt1=1,TT

!		int c0=(tt1-1) % 3;
    c1= mod(tt1, 3)
    c2= mod(tt1+1, 3)
!		long double omega_off=omegae+2*pi*300000;

		! update display every 50th iteration
    if (mod(tt1,50) == 0) print '(A,I0.3,F5.2,A,I0.3,A,I0.3)',"Realization: ",&
        realization,tt1*100.0/TT,"% complete.  n",beami," v",beamj

    do pp=1,N

      LL = max(p(pp)-N/3,-N/3)
      UU = min(N/3,p(pp)+N/3)
      CC(1)=0.0;
      CC(2)=0.0;

      do q=LL,UU
        CC(1)=CC(1)+EE(c1,q+N/2,1)*nn(c1,p(pp)-q+N/2,1)-EE(c1,q+N/2,2)*nn(c1,p(pp)-q+N/2,2)
        CC(2)=CC(2)+EE(c1,q+N/2,1)*nn(c1,p(pp)-q+N/2,2)+EE(c1,q+N/2,2)*nn(c1,p(pp)-q+N/2,1);
      end do

      call random_number(rdist)
      SSE(pp,1) = rdist*Source_factor_E(pp)/sqrt(Tstep)
      call random_number(rdist)
      SSE(pp,2)= rdist*Source_factor_E(pp)/sqrt(Tstep)

      cte1=1.5*omegae*(lambdaD*k(pp))*(lambdaD*k(pp))
			!cte1=1.5*Kb*Te/me/omega_off*k(pp)*k(pp)-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
      k1(pp,1)=Tstep*(cte1*EE(c1,pp,2)-nuE(pp)*EE(c1,pp,1)+cte2*CC(2))
      k1(pp,2)=Tstep*((-1)*cte1*EE(c1,pp,1)-nuE(pp)*EE(c1,pp,2)-cte2*CC(1))
    end do ! pp N


    do pp=1,N
      LL= max(p(pp)-N/3,-N/3)
      UU= min(N/3,p(pp)+N/3)
      CC(1)=0.0
      CC(2)=0.0

      do q=LL,UU
        CC(1)=CC(1)+(EE(c1,q+N/2,1)+k1(q+N/2,1)/2)*nn(c1,p(pp)-q+N/2,1)-(EE(c1,q+N/2,2)+k1(q+N/2,2)/2)*nn(c1,p(pp)-q+N/2,2)
        CC(2)=CC(2)+(EE(c1,q+N/2,1)+k1(q+N/2,1)/2)*nn(c1,p(pp)-q+N/2,2)+(EE(c1,q+N/2,2)+k1(q+N/2,2)/2)*nn(c1,p(pp)-q+N/2,1)
      end do

      cte1=1.5*omegae*(lambdaD*k(pp))*(lambdaD*k(pp))
			!cte1=1.5*Kb*Te/me/omega_off*k(pp)*k(pp)-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
      k2(pp,1)=Tstep*(cte1*(EE(c1,pp,2)+k1(pp,2)/2.0-SSE(pp,1)/2.0*Tstep) - &
                nuE(pp) * (EE(c1,pp,1)+k1(pp,1)/2.0+SSE(pp,2)/2.0*Tstep)+cte2*CC(2))
      k2(pp,2)=Tstep*((-1)*cte1*(EE(c1,pp,1)+k1(pp,1)/2.0+SSE(pp,2)/2.0*Tstep) - &
                nuE(pp)*(EE(c1,pp,2)+k1(pp,2)/2.0-SSE(pp,1)/2.0*Tstep)-cte2*CC(1))

    end do ! pp N


    do pp=1,N
      LL= max(p(pp)-N/3,-N/3)
      UU= min(N/3,p(pp)+N/3)
      CC(1)=0.0
      CC(2)=0.0

      do q=LL,UU
        CC(1)=CC(1)+(EE(c1,q+N/2,1)+k2(q+N/2,1)/2)*nn(c1,p(pp)-q+N/2,1)-(EE(c1,q+N/2,2)+k2(q+N/2,2)/2)*nn(c1,p(pp)-q+N/2,2);
        CC(2)=CC(2)+(EE(c1,q+N/2,1)+k2(q+N/2,1)/2)*nn(c1,p(pp)-q+N/2,2)+(EE(c1,q+N/2,2)+k2(q+N/2,2)/2)*nn(c1,p(pp)-q+N/2,1);
      end do

      cte1=1.5*omegae*(lambdaD*k(pp))*(lambdaD*k(pp))
			!cte1=1.5*Kb*Te/me/omega_off*k(pp)*k(pp)-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
      k3(pp,1)=Tstep*(cte1*(EE(c1,pp,2)+k2(pp,2)/2.0-SSE(pp,1)/2.0*Tstep) - &
              nuE(pp)*(EE(c1,pp,1)+k2(pp,1)/2.0+SSE(pp,2)/2.0*Tstep)+cte2*CC(2));
      k3(pp,2)=Tstep*((-1)*cte1*(EE(c1,pp,1)+k2(pp,1)/2.0+SSE(pp,2)/2.0*Tstep) - &
              nuE(pp)*(EE(c1,pp,2)+k2(pp,2)/2.0-SSE(pp,1)/2.0*Tstep)-cte2*CC(1));

    end do ! pp N


    do pp=1,N
      LL= max(p(pp)-N/3,-N/3);
      UU= min(N/3,p(pp)+N/3);
      CC(1)=0.0;
      CC(2)=0.0;

      do q=LL,UU
        CC(1)=CC(1)+(EE(c1,q+N/2,1)+k3(q+N/2,1))*nn(c1,p(pp)-q+N/2,1)-(EE(c1,q+N/2,2)+k3(q+N/2,2))*nn(c1,p(pp)-q+N/2,2);
        CC(2)=CC(2)+(EE(c1,q+N/2,1)+k3(q+N/2,1))*nn(c1,p(pp)-q+N/2,2)+(EE(c1,q+N/2,2)+k3(q+N/2,2))*nn(c1,p(pp)-q+N/2,1);
      end do


      cte1=1.5*omegae*(lambdaD*k(pp))*(lambdaD*k(pp))
			!cte1=1.5*Kb*Te/me/omega_off*k(pp)*k(pp)-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
      k4(pp,1)=Tstep*(cte1*(EE(c1,pp,2)+k3(pp,2)-SSE(pp,1)*Tstep)-nuE(pp)*(EE(c1,pp,1)+k3(pp,1)+SSE(pp,2)*Tstep)+cte2*CC(2));
      k4(pp,2)=Tstep*((-1)*cte1*(EE(c1,pp,1)+k3(pp,1)+SSE(pp,2)*Tstep)-nuE(pp)*(EE(c1,pp,2)+k3(pp,2)-SSE(pp,1)*Tstep)-cte2*CC(1));

      EE(c2,pp,1)=EE(c1,pp,1)+(k1(pp,1)+2.0*k2(pp,1)+2.0*k3(pp,1)+k4(pp,1))/6.0+SSE(pp,2)*Tstep;
      EE(c2,pp,2)=EE(c1,pp,2)+(k1(pp,2)+2.0*k2(pp,2)+2.0*k3(pp,2)+k4(pp,2))/6.0-SSE(pp,1)*Tstep;
      EE(c2,N/2,1)=0.0
      EE(c2,N/2,2)=0.0


    end do ! pp N


    do pp=1,N/2
      call random_number(rdist)
      SSn(1)= rdist*Source_factor_n(pp)/sqrt(Tstep);
      call random_number(rdist)
      SSn(2)= rdist*Source_factor_n(pp)/sqrt(Tstep);


      LL= max(p(pp)-N/3,-N/3);
      UU= min(N/3,p(pp)+N/3);
      CC(1)=0.0;
      CC(2)=0.0;

      do q=LL,UU
        CC(1)=CC(1)+EE(c1,q+N/2,1)*EE(c1,q-p(pp)+N/2,1)+EE(c1,q+N/2,2)*EE(c1,q-p(pp)+N/2,2);
        CC(2)=CC(2)+EE(c1,q+N/2,2)*EE(c1,q-p(pp)+N/2,1)-EE(c1,q+N/2,1)*EE(c1,q-p(pp)+N/2,2);
      end do

      kn1(1)=Tstep*(vv(c1,pp,1));
      kn1(2)=Tstep*(vv(c1,pp,2));
      kv1(1)=Tstep*((-2.0)*nui(pp)*vv(c1,pp,1)- Cs*k(pp)**2 *nn(c1,pp,1)-k(pp)*k(pp)*epsilon0/4/mi*CC(1));
      kv1(2)=Tstep*((-2.0)*nui(pp)*vv(c1,pp,2)- Cs*k(pp)**2 *nn(c1,pp,2)-k(pp)*k(pp)*epsilon0/4/mi*CC(2));

      CC(1)=0.0;
      CC(2)=0.0;

      do q=LL,UU
        CC(1)=CC(1)+(EE(c1,q+N/2,1)+k1(q+N/2,1)/2)*(EE(c1,q-p(pp)+N/2,1)+k1(q-p(pp)+N/2,1)/2)+ &
              (EE(c1,q+N/2,2)+k1(q+N/2,2)/2)*(EE(c1,q-p(pp)+N/2,2)+k1(q-p(pp)+N/2,2)/2);
        CC(2)=CC(2)+(EE(c1,q+N/2,2)+k1(q+N/2,2)/2)*(EE(c1,q-p(pp)+N/2,1)+k1(q-p(pp)+N/2,1)/2)- &
              (EE(c1,q+N/2,1)+k1(q+N/2,1)/2)*(EE(c1,q-p(pp)+N/2,2)+k1(q-p(pp)+N/2,2)/2);
      end do

      kn2(1)=Tstep*(vv(c1,pp,1)+kv1(1)/2+SSn(1)/2*Tstep);
      kn2(2)=Tstep*(vv(c1,pp,2)+kv1(2)/2+SSn(2)/2*Tstep);
      kv2(1)=Tstep*((-2.0)*nui(pp)*(vv(c1,pp,1)+kv1(1)/2+SSn(1)/2*Tstep)- &
              Cs*k(pp)**2 *(nn(c1,pp,1)+kn1(1)/2)-k(pp)*k(pp)*epsilon0/4/mi*CC(1))
      kv2(2)=Tstep*((-2.0)*nui(pp)*(vv(c1,pp,2)+kv1(2)/2+SSn(2)/2*Tstep)- &
              Cs*k(pp)**2 *(nn(c1,pp,2)+kn1(2)/2)-k(pp)*k(pp)*epsilon0/4/mi*CC(2))

      CC(1)=0.0;
      CC(2)=0.0;

      do q=LL,UU
        CC(1)=CC(1)+(EE(c1,q+N/2,1)+k2(q+N/2,1)/2)*(EE(c1,q-p(pp)+N/2,1)+k2(q-p(pp)+N/2,1)/2) + &
              (EE(c1,q+N/2,2)+k2(q+N/2,2)/2)*(EE(c1,q-p(pp)+N/2,2)+k2(q-p(pp)+N/2,2)/2);
        CC(2)=CC(2)+(EE(c1,q+N/2,2)+k2(q+N/2,2)/2)*(EE(c1,q-p(pp)+N/2,1)+k2(q-p(pp)+N/2,1)/2) - &
              (EE(c1,q+N/2,1)+k2(q+N/2,1)/2)*(EE(c1,q-p(pp)+N/2,2)+k2(q-p(pp)+N/2,2)/2);
      end do

      kn3(1)=Tstep*(vv(c1,pp,1)+kv2(1)/2+SSn(1)/2*Tstep)
      kn3(2)=Tstep*(vv(c1,pp,2)+kv2(2)/2+SSn(2)/2*Tstep)
      kv3(1)=Tstep*((-2.0)*nui(pp)*(vv(c1,pp,1)+kv2(1)/2+SSn(1)/2*Tstep) - &
             Cs*k(pp)**2 *(nn(c1,pp,1)+kn2(1)/2)-k(pp)*k(pp)*epsilon0/4/mi*CC(1))
      kv3(2)=Tstep*((-2.0)*nui(pp)*(vv(c1,pp,2)+kv2(2)/2+SSn(2)/2*Tstep) - &
             Cs*k(pp)**2 *(nn(c1,pp,2)+kn2(2)/2)-k(pp)*k(pp)*epsilon0/4/mi*CC(2))


      CC(1)=0.0;
      CC(2)=0.0;
      do q=LL,UU
        CC(1)=CC(1)+(EE(c1,q+N/2,1)+k3(q+N/2,1)/2)*(EE(c1,q-p(pp)+N/2,1)+k3(q-p(pp)+N/2,1)/2) + &
              (EE(c1,q+N/2,2)+k3(q+N/2,2)/2)*(EE(c1,q-p(pp)+N/2,2)+k3(q-p(pp)+N/2,2)/2);
        CC(2)=CC(2)+(EE(c1,q+N/2,2)+k3(q+N/2,2)/2)*(EE(c1,q-p(pp)+N/2,1)+k3(q-p(pp)+N/2,1)/2) - &
              (EE(c1,q+N/2,1)+k3(q+N/2,1)/2)*(EE(c1,q-p(pp)+N/2,2)+k3(q-p(pp)+N/2,2)/2);
      end do

      kn4(1)=Tstep*(vv(c1,pp,1)+kv3(1)+SSn(1)*Tstep);
      kn4(2)=Tstep*(vv(c1,pp,2)+kv3(2)+SSn(2)*Tstep);
      kv4(1)=Tstep*((-2.0)*nui(pp)*(vv(c1,pp,1)+kv3(1)+SSn(1)*Tstep) - &
             Cs*k(pp)**2 *(nn(c1,pp,1)+kn3(1))-k(pp)*k(pp)*epsilon0/4/mi*CC(1));
      kv4(2)=Tstep*((-2.0)*nui(pp)*(vv(c1,pp,2)+kv3(2)+SSn(2)*Tstep) - &
             Cs*k(pp)**2 *(nn(c1,pp,2)+kn3(2))-k(pp)*k(pp)*epsilon0/4/mi*CC(2));


      vv(c2,pp,1)=vv(c1,pp,1)+(kv1(1)+2*kv2(1)+2*kv3(1)+kv4(1))/6+SSn(1)*Tstep;
      vv(c2,pp,2)=vv(c1,pp,2)+(kv1(2)+2*kv2(2)+2*kv3(2)+kv4(2))/6+SSn(2)*Tstep;
      nn(c2,pp,1)=nn(c1,pp,1)+(kn1(1)+2*kn2(1)+2*kn3(1)+kn4(1))/6;
      nn(c2,pp,2)=nn(c1,pp,2)+(kn1(2)+2*kn2(2)+2*kn3(2)+kn4(2))/6;
      nn(c2,N/2,1)=0.0;
      nn(c2,N/2,2)=0.0;

      if (pp>=1) then
        nn(c2,N-pp,1)=nn(c2,pp,1);
        nn(c2,N-pp,2)=-nn(c2,pp,2);
      end if

    end do ! pp N/2

    if ( mod(tt1,res) == 1) then
      do pp=1,N
        total_EE(counter1*N*2+pp*2+1)=EE(c2,pp,1);
        total_EE(counter1*N*2+pp*2+2)=EE(c2,pp,2);
        total_nn(counter1*N*2+pp*2+1)=nn(c2,pp,1);
        total_nn(counter1*N*2+pp*2+2)=nn(c2,pp,2);
      end do ! pp N
      counter1 = counter1 + 1
    end if

    if (counter1==20000) then
      write(uEE,*) total_EE
      write(unn,*) total_nn
      counter1=0
    end if

  end do ! tt1

  if (counter1>0) then
    write(uEE,*) total_EE
    write(uNN,*) total_nn
  end if


  close(uEE)
  close(unn)
end do !relizations
end do !Nvbeam
end do !Nnbeam


call cpu_time(toc)
print *,"Elapsed Time: ", toc-tic


contains


pure elemental subroutine Xsection(Xsec_ion, Xsec_pl, k)

  real(wp), intent(out) :: Xsec_ion, Xsec_pl
  real(wp), intent(in) :: k

  real(wp) :: alpha,XX

  alpha=1.0_wp / (k*lambdaD)
  Xsec_ion= 2.0_wp * pi/(1.0_wp + alpha**2.0_wp)*(Z* alpha**4.0_wp/(1.0_wp+alpha**2.0_wp + alpha**2.0_wp * (Z*Te/Ti) ))

  XX=2.0_wp*pi*(1.0_wp + alpha**2.0_wp*Z*Te/Ti)/(1.0_wp+alpha**2.0_wp+alpha**2.0_wp*(Z*Te/Ti))

  Xsec_pl=XX-Te/Ti*Xsec_ion

end subroutine Xsection


end program
