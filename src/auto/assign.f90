! Automatically generated file.  Do not edit!

  subroutine assign(nline,nvalues,var,value,values)

  use param

  use procinfo

  implicit none

  logical contains
  integer i,nvalues
  integer(kind=2) nline
  character(50) type
  character(1000) range
  character(len=*) var,value
  character(len=*) values(1:nvalues)

  if (var=='parfile') then
     type = 'CHARACTER'
     parfile = value
     if (nvalues>1) goto 200
     return
  end if

  if (var=='dr0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dr0
     return
  end if

  if (var=='rbound') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) rbound
     return
  end if

  if (var=='Nrtotal') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Nrtotal
     return
  end if

  if (var=='Nr') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Nr
     return
  end if

  if (var=='Nl') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Nl
     return
  end if

  if (var=='ghost') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) ghost
     return
  end if

  if (var=='intorder') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) intorder
     return
  end if

  if (var=='nproctot') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) nproctot
     return
  end if

  if (var=='adjuststep') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) adjuststep
     if (nvalues>1) goto 200
     return
  end if

  if (var=='adjuststepmax') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) adjuststepmax
     if (nvalues>1) goto 200
     return
  end if

  if (var=='dt0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dt0
     return
  end if

  if (var=='dtfac') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dtfac
     return
  end if

  if (var=='Nt') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Nt
     return
  end if

  if (var=='directory') then
     type = 'CHARACTER'
     directory = value
     if (nvalues>1) goto 200
     return
  end if

  if (var=='checkpointfile') then
     type = 'CHARACTER'
     checkpointfile = value
     if (nvalues>1) goto 200
     return
  end if

  if (var=='checkpointinitial') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) checkpointinitial
     if (nvalues>1) goto 200
     return
  end if

  if (var=='checkpoint') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) checkpoint
     if (nvalues>1) goto 200
     return
  end if

  if (var=='closefiles') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) closefiles
     if (nvalues>1) goto 200
     return
  end if

  if (var=='Ninfo') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Ninfo
     return
  end if

  if (var=='Noutput0D') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Noutput0D
     return
  end if

  if (var=='Noutput1D') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Noutput1D
     return
  end if

  if (var=='Ncheckpoint') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) Ncheckpoint
     return
  end if

  if (var=='nvars0D') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) nvars0D
     return
  end if

  if (var=='nvars0D') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nvars0D
     return
  end if

  if (var=='nvars1D') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) nvars1D
     return
  end if

  if (var=='nvars1D') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nvars1D
     return
  end if

  if (var=='norm_rmin') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) norm_rmin
     return
  end if

  if (var=='output_r1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) output_r1
     return
  end if

  if (var=='outvars0D') then
     type = 'CHARACTER'
     outvars0D = value
     return
  end if

  if (var=='outvars1D') then
     type = 'CHARACTER'
     outvars1D = value
     return
  end if

  if (var=='checkvars') then
     type = 'CHARACTER'
     checkvars = value
     return
  end if

  if (var=='commenttype') then
     type = 'CHARACTER'
     commenttype = value
     if (nvalues>1) goto 200
     range = &
     '(xgraph,gnuplot)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='boundtype') then
     type = 'CHARACTER'
     boundtype = value
     if (nvalues>1) goto 200
     range = &
     '(none,static,flat,radiative,constraint)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='dorigin') then
     type = 'CHARACTER'
     dorigin = value
     if (nvalues>1) goto 200
     range = &
     '(centered,onesided)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='slicing') then
     type = 'CHARACTER'
     slicing = value
     if (nvalues>1) goto 200
     range = &
     '(static,maximal,harmonic,1+log,shockavoid,alphaminus2,cosmocf-harmonic,cosmocf-1' // &
     '+log,cosmocf-shockavoid,cosmosync-harmonic,cosmosync-1+log,cosmosync-shockavoid)' // &
     ''
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='ilapse') then
     type = 'CHARACTER'
     ilapse = value
     if (nvalues>1) goto 200
     range = &
     '(none,one,isotropic,psiminus2,psiminus4,psi2,psi4,maximal)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='maximalbound') then
     type = 'CHARACTER'
     maximalbound = value
     if (nvalues>1) goto 200
     range = &
     '(robin,dirichlet,conformal)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='gauge_f') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) gauge_f
     return
  end if

  if (var=='lapsepert') then
     type = 'CHARACTER'
     lapsepert = value
     if (nvalues>1) goto 200
     range = &
     '(none,gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='lapse_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) lapse_a0
     return
  end if

  if (var=='lapse_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) lapse_r0
     return
  end if

  if (var=='lapse_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) lapse_s0
     return
  end if

  if (var=='lapse_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) lapse_t0
     return
  end if

  if (var=='lapseeta') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) lapseeta
     return
  end if

  if (var=='lapsediss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) lapsediss
     return
  end if

  if (var=='shift') then
     type = 'CHARACTER'
     shift = value
     if (nvalues>1) goto 200
     range = &
     '(none,zero,static,Gammadriver0,Gammadriver1,Gammadriver2,Gammadriver3,Gammadrive' // &
     'rshock1,Gammadrivershock2)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='driverD0') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) driverD0
     if (nvalues>1) goto 200
     return
  end if

  if (var=='drivercsi') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) drivercsi
     return
  end if

  if (var=='drivereta') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) drivereta
     return
  end if

  if (var=='shiftpert') then
     type = 'CHARACTER'
     shiftpert = value
     if (nvalues>1) goto 200
     range = &
     '(none,gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='shift_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) shift_a0
     return
  end if

  if (var=='shift_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) shift_r0
     return
  end if

  if (var=='shift_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) shift_s0
     return
  end if

  if (var=='shift_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) shift_t0
     return
  end if

  if (var=='rescaledata') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) rescaledata
     return
  end if

  if (var=='idata') then
     type = 'CHARACTER'
     idata = value
     if (nvalues>1) goto 200
     range = &
     '(checkpoint,minkowski,schwarzschild,schwarzschildKS,trumpetBH,reissnernordstrom,' // &
     'desitter,scalarpulse,ghostpulse,nonminpulse,complexpulse,complexghostpulse,boson' // &
     'star,chargedboson,procapulse,procastar,l-procastar,chargedproca,diracpulse,dirac' // &
     'star,dustshell,fluidshell,TOVstar,TOVcomplex,blastwave,scalarDM,complexDM,ghostw' // &
     'ormhole,duststep)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='spacetime') then
     type = 'CHARACTER'
     spacetime = value
     if (nvalues>1) goto 200
     range = &
     '(dynamic,background,minkowski)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='formulation') then
     type = 'CHARACTER'
     formulation = value
     if (nvalues>1) goto 200
     range = &
     '(bssn,z4c)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='bssnflavor') then
     type = 'CHARACTER'
     bssnflavor = value
     if (nvalues>1) goto 200
     range = &
     '(eulerian,lagrangian)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='integrator') then
     type = 'CHARACTER'
     integrator = value
     if (nvalues>1) goto 200
     range = &
     '(icn,rk4)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='order') then
     type = 'CHARACTER'
     order = value
     if (nvalues>1) goto 200
     range = &
     '(two,four,six,eight)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='nolambda') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) nolambda
     if (nvalues>1) goto 200
     return
  end if

  if (var=='noDeltar') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) noDeltar
     if (nvalues>1) goto 200
     return
  end if

  if (var=='noKTA') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) noKTA
     if (nvalues>1) goto 200
     return
  end if

  if (var=='chimethod') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) chimethod
     if (nvalues>1) goto 200
     return
  end if

  if (var=='regular2') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) regular2
     if (nvalues>1) goto 200
     return
  end if

  if (var=='icniter') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) icniter
     return
  end if

  if (var=='chipower') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) chipower
     return
  end if

  if (var=='lambdapower') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) lambdapower
     return
  end if

  if (var=='eta') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) eta
     return
  end if

  if (var=='kappa1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) kappa1
     return
  end if

  if (var=='kappa2') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) kappa2
     return
  end if

  if (var=='geodiss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) geodiss
     return
  end if

  if (var=='scalardiss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalardiss
     return
  end if

  if (var=='nonmindiss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmindiss
     return
  end if

  if (var=='elecdiss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) elecdiss
     return
  end if

  if (var=='procadiss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) procadiss
     return
  end if

  if (var=='diracdiss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracdiss
     return
  end if

  if (var=='fluiddiss') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluiddiss
     return
  end if

  if (var=='ahfound') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) ahfound
     if (nvalues>1) goto 200
     return
  end if

  if (var=='ahfind') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) ahfind
     if (nvalues>1) goto 200
     return
  end if

  if (var=='ahfind_every') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) ahfind_every
     return
  end if

  if (var=='ahafter') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ahafter
     return
  end if

  if (var=='ahmass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ahmass
     return
  end if

  if (var=='mattertype') then
     type = 'CHARACTER'
     mattertype = value
     range = &
     '(vacuum,cosmo,scalar,ghost,complex,complexghost,nonmin,electric,proca,complexpro' // &
     'ca,dirac,fluid,dust)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='cosmic_run') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) cosmic_run
     if (nvalues>1) goto 200
     return
  end if

  if (var=='lambda_cosmo') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) lambda_cosmo
     return
  end if

  if (var=='BHmass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) BHmass
     return
  end if

  if (var=='BHcharge') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) BHcharge
     return
  end if

  if (var=='scalar_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_a0
     return
  end if

  if (var=='scalar_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_r0
     return
  end if

  if (var=='scalar_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_s0
     return
  end if

  if (var=='scalar_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_t0
     return
  end if

  if (var=='scalar_mass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_mass
     return
  end if

  if (var=='scalar_lambda') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_lambda
     return
  end if

  if (var=='scalarprofile') then
     type = 'CHARACTER'
     scalarprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,r2gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='scalarpotential') then
     type = 'CHARACTER'
     scalarpotential = value
     if (nvalues>1) goto 200
     range = &
     '(none,phi2,phi4)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='scalarmethod') then
     type = 'CHARACTER'
     scalarmethod = value
     if (nvalues>1) goto 200
     range = &
     '(first,second)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='scalar_relax') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) scalar_relax
     if (nvalues>1) goto 200
     return
  end if

  if (var=='scalar_bg_pert') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) scalar_bg_pert
     if (nvalues>1) goto 200
     return
  end if

  if (var=='scalar_bg_phi0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_bg_phi0
     return
  end if

  if (var=='scalar_bg_pi0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) scalar_bg_pi0
     return
  end if

  if (var=='complexR_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_a0
     return
  end if

  if (var=='complexR_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_r0
     return
  end if

  if (var=='complexR_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_s0
     return
  end if

  if (var=='complexR_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_t0
     return
  end if

  if (var=='complexR_a1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_a1
     return
  end if

  if (var=='complexR_r1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_r1
     return
  end if

  if (var=='complexR_s1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_s1
     return
  end if

  if (var=='complexR_t1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexR_t1
     return
  end if

  if (var=='complexI_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_a0
     return
  end if

  if (var=='complexI_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_r0
     return
  end if

  if (var=='complexI_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_s0
     return
  end if

  if (var=='complexI_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_t0
     return
  end if

  if (var=='complexI_a1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_a1
     return
  end if

  if (var=='complexI_r1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_r1
     return
  end if

  if (var=='complexI_s1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_s1
     return
  end if

  if (var=='complexI_t1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexI_t1
     return
  end if

  if (var=='complex_mass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complex_mass
     return
  end if

  if (var=='complex_lambda') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complex_lambda
     return
  end if

  if (var=='complex_l') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) complex_l
     return
  end if

  if (var=='complex_lmax') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) complex_lmax
     return
  end if

  if (var=='k_parameter') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) k_parameter
     return
  end if

  if (var=='complexprofile') then
     type = 'CHARACTER'
     complexprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tophat,comp-gaussian,comp-tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='complexpotential') then
     type = 'CHARACTER'
     complexpotential = value
     if (nvalues>1) goto 200
     range = &
     '(none,phi2,phi4)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='complexmethod') then
     type = 'CHARACTER'
     complexmethod = value
     if (nvalues>1) goto 200
     range = &
     '(first,second)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='complexDM_type') then
     type = 'CHARACTER'
     complexDM_type = value
     if (nvalues>1) goto 200
     range = &
     '(harmonicTS,harmonicMOM,growing)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='boson_phi0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_phi0
     return
  end if

  if (var=='boson_omega') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_omega
     return
  end if

  if (var=='omega_left') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) omega_left
     return
  end if

  if (var=='omega_right') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) omega_right
     return
  end if

  if (var=='boson_factor') then
     type = 'CHARACTER'
     boson_factor = value
     if (nvalues>1) goto 200
     range = &
     '(physical,harmonic)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='boson_gauge') then
     type = 'CHARACTER'
     boson_gauge = value
     if (nvalues>1) goto 200
     range = &
     '(PA,CF)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='boson_relax') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) boson_relax
     if (nvalues>1) goto 200
     return
  end if

  if (var=='bosongauss') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) bosongauss
     if (nvalues>1) goto 200
     return
  end if

  if (var=='boson_phiR_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_phiR_a0
     return
  end if

  if (var=='boson_phiR_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_phiR_r0
     return
  end if

  if (var=='boson_phiR_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_phiR_s0
     return
  end if

  if (var=='boson_piI_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_piI_a0
     return
  end if

  if (var=='boson_piI_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_piI_r0
     return
  end if

  if (var=='boson_piI_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) boson_piI_s0
     return
  end if

  if (var=='complex_bg_pert') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) complex_bg_pert
     if (nvalues>1) goto 200
     return
  end if

  if (var=='complex_bg_phiR0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complex_bg_phiR0
     return
  end if

  if (var=='complex_bg_phiI0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complex_bg_phiI0
     return
  end if

  if (var=='complex_initialcharged') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) complex_initialcharged
     if (nvalues>1) goto 200
     return
  end if

  if (var=='complex_q') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complex_q
     return
  end if

  if (var=='charge_factor') then
     type = 'CHARACTER'
     charge_factor = value
     if (nvalues>1) goto 200
     range = &
     '(standard,jetzer)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='ghost_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ghost_a0
     return
  end if

  if (var=='ghost_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ghost_r0
     return
  end if

  if (var=='ghost_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ghost_s0
     return
  end if

  if (var=='ghost_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ghost_t0
     return
  end if

  if (var=='ghost_mass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ghost_mass
     return
  end if

  if (var=='ghost_lambda') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) ghost_lambda
     return
  end if

  if (var=='ghostprofile') then
     type = 'CHARACTER'
     ghostprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='ghostpotential') then
     type = 'CHARACTER'
     ghostpotential = value
     if (nvalues>1) goto 200
     range = &
     '(none,phi2,phi4)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='ghostmethod') then
     type = 'CHARACTER'
     ghostmethod = value
     if (nvalues>1) goto 200
     range = &
     '(first,second)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='complexghostR_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostR_a0
     return
  end if

  if (var=='complexghostR_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostR_r0
     return
  end if

  if (var=='complexghostR_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostR_s0
     return
  end if

  if (var=='complexghostR_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostR_t0
     return
  end if

  if (var=='complexghostI_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostI_a0
     return
  end if

  if (var=='complexghostI_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostI_r0
     return
  end if

  if (var=='complexghostI_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostI_s0
     return
  end if

  if (var=='complexghostI_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghostI_t0
     return
  end if

  if (var=='complexghost_mass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghost_mass
     return
  end if

  if (var=='complexghost_lambda') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) complexghost_lambda
     return
  end if

  if (var=='complexghost_l') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) complexghost_l
     return
  end if

  if (var=='complexghost_lmax') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) complexghost_lmax
     return
  end if

  if (var=='complexghostprofile') then
     type = 'CHARACTER'
     complexghostprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='complexghostpotential') then
     type = 'CHARACTER'
     complexghostpotential = value
     if (nvalues>1) goto 200
     range = &
     '(none,phi2,phi4)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='complexghostmethod') then
     type = 'CHARACTER'
     complexghostmethod = value
     if (nvalues>1) goto 200
     range = &
     '(first,second)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='nonmin_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmin_a0
     return
  end if

  if (var=='nonmin_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmin_r0
     return
  end if

  if (var=='nonmin_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmin_s0
     return
  end if

  if (var=='nonmin_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmin_t0
     return
  end if

  if (var=='nonmin_fxi') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmin_fxi
     return
  end if

  if (var=='nonmin_theta') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmin_theta
     return
  end if

  if (var=='nonmin_phi0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) nonmin_phi0
     return
  end if

  if (var=='nonminf') then
     type = 'CHARACTER'
     nonminf = value
     if (nvalues>1) goto 200
     range = &
     '(trivial,quadratic,BransDicke)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='nonminprofile') then
     type = 'CHARACTER'
     nonminprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,r2gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='nonminpotential') then
     type = 'CHARACTER'
     nonminpotential = value
     if (nvalues>1) goto 200
     range = &
     '(none)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='nonminmethod') then
     type = 'CHARACTER'
     nonminmethod = value
     if (nvalues>1) goto 200
     range = &
     '(first,second)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='nonmin_pulse') then
     type = 'CHARACTER'
     nonmin_pulse = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tanh)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='proca_mass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_mass
     return
  end if

  if (var=='proca_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_a0
     return
  end if

  if (var=='proca_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_r0
     return
  end if

  if (var=='proca_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_s0
     return
  end if

  if (var=='proca_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_t0
     return
  end if

  if (var=='procaprofile') then
     type = 'CHARACTER'
     procaprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='cproca_mass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) cproca_mass
     return
  end if

  if (var=='cproca_q') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) cproca_q
     return
  end if

  if (var=='cproca_l') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) cproca_l
     return
  end if

  if (var=='proca_phi0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_phi0
     return
  end if

  if (var=='proca_omega') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_omega
     return
  end if

  if (var=='proca_factor') then
     type = 'CHARACTER'
     proca_factor = value
     if (nvalues>1) goto 200
     range = &
     '(physical,harmonic)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='procagauss') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) procagauss
     if (nvalues>1) goto 200
     return
  end if

  if (var=='proca_PhiR_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_PhiR_a0
     return
  end if

  if (var=='proca_PhiR_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_PhiR_r0
     return
  end if

  if (var=='proca_PhiR_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_PhiR_s0
     return
  end if

  if (var=='proca_AI_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_AI_a0
     return
  end if

  if (var=='proca_AI_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_AI_r0
     return
  end if

  if (var=='proca_AI_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) proca_AI_s0
     return
  end if

  if (var=='dirac_mass') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dirac_mass
     return
  end if

  if (var=='diracprofile') then
     type = 'CHARACTER'
     diracprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='diractype') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) diractype
     return
  end if

  if (var=='dirac_k') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dirac_k
     return
  end if

  if (var=='diracFR_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFR_a0
     return
  end if

  if (var=='diracFR_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFR_r0
     return
  end if

  if (var=='diracFR_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFR_s0
     return
  end if

  if (var=='diracFR_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFR_t0
     return
  end if

  if (var=='diracFI_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFI_a0
     return
  end if

  if (var=='diracFI_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFI_r0
     return
  end if

  if (var=='diracFI_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFI_s0
     return
  end if

  if (var=='diracFI_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracFI_t0
     return
  end if

  if (var=='diracGR_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGR_a0
     return
  end if

  if (var=='diracGR_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGR_r0
     return
  end if

  if (var=='diracGR_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGR_s0
     return
  end if

  if (var=='diracGR_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGR_t0
     return
  end if

  if (var=='diracGI_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGI_a0
     return
  end if

  if (var=='diracGI_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGI_r0
     return
  end if

  if (var=='diracGI_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGI_s0
     return
  end if

  if (var=='diracGI_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) diracGI_t0
     return
  end if

  if (var=='dirac_f0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dirac_f0
     return
  end if

  if (var=='dirac_omega') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dirac_omega
     return
  end if

  if (var=='diracgauss') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) diracgauss
     if (nvalues>1) goto 200
     return
  end if

  if (var=='dust_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dust_a0
     return
  end if

  if (var=='dust_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dust_r0
     return
  end if

  if (var=='dust_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dust_s0
     return
  end if

  if (var=='dust_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dust_t0
     return
  end if

  if (var=='dust_atmos') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) dust_atmos
     return
  end if

  if (var=='dustprofile') then
     type = 'CHARACTER'
     dustprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='dust_method') then
     type = 'CHARACTER'
     dust_method = value
     if (nvalues>1) goto 200
     range = &
     '(center,upwind,limiter,mp5)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='fluid_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_a0
     return
  end if

  if (var=='fluid_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_r0
     return
  end if

  if (var=='fluid_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_s0
     return
  end if

  if (var=='fluid_t0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_t0
     return
  end if

  if (var=='fluid_N') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_N
     return
  end if

  if (var=='fluid_gamma') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_gamma
     return
  end if

  if (var=='fluid_kappa') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_kappa
     return
  end if

  if (var=='fluid_atmos') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_atmos
     return
  end if

  if (var=='fluid_q1') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_q1
     return
  end if

  if (var=='fluid_q2') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) fluid_q2
     return
  end if

  if (var=='fluidprofile') then
     type = 'CHARACTER'
     fluidprofile = value
     if (nvalues>1) goto 200
     range = &
     '(gaussian,tophat)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='fluid_EOS') then
     type = 'CHARACTER'
     fluid_EOS = value
     if (nvalues>1) goto 200
     range = &
     '(none,ideal)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='fluid_method') then
     type = 'CHARACTER'
     fluid_method = value
     if (nvalues>1) goto 200
     range = &
     '(llf,hlle,mp5)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='fluid_limiter') then
     type = 'CHARACTER'
     fluid_limiter = value
     if (nvalues>1) goto 200
     range = &
     '(minmod,vanleer,superbee,mc,koren,ospre,sweby)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='fluid_usesoundspeed') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) fluid_usesoundspeed
     if (nvalues>1) goto 200
     return
  end if

  if (var=='TOV_rho0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) TOV_rho0
     return
  end if

  if (var=='TOV_rad') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) TOV_rad
     return
  end if

  if (var=='TOVgauss') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) TOVgauss
     if (nvalues>1) goto 200
     return
  end if

  if (var=='TOV_a0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) TOV_a0
     return
  end if

  if (var=='TOV_r0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) TOV_r0
     return
  end if

  if (var=='TOV_s0') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) TOV_s0
     return
  end if

  if (var=='blast_R') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) blast_R
     return
  end if

  if (var=='blast_rhol') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) blast_rhol
     return
  end if

  if (var=='blast_rhor') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) blast_rhor
     return
  end if

  if (var=='blast_pl') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) blast_pl
     return
  end if

  if (var=='blast_pr') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) blast_pr
     return
  end if

  if (var=='newr') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) newr
     if (nvalues>1) goto 200
     return
  end if

  if (var=='rtype') then
     type = 'CHARACTER'
     rtype = value
     if (nvalues>1) goto 200
     range = &
     '(sigmoid,choptuik,smoothstep,polynomal,sinh)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='drinf') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) drinf
     return
  end if

  if (var=='beta_transf') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) beta_transf
     return
  end if

  if (var=='gamma_transf') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) gamma_transf
     return
  end if

  if (var=='delta_transf') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) delta_transf
     return
  end if

  if (var=='integralmethod') then
     type = 'CHARACTER'
     integralmethod = value
     if (nvalues>1) goto 200
     range = &
     '(quadrature,fourth)'
     do i=1,nvalues
        if ((.not.contains(range,'('//trim(values(i))//')')).and. &
            (.not.contains(range,'('//trim(values(i)))).and. &
            (.not.contains(range,trim(values(i))//')')).and. &
            (.not.contains(range,trim(values(i))))) then
           if (rank==0) then
              print *
              print *, 'Parfile error.'
              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline
              print *
              print *, 'Aborting! (subroutine assign.f90)'
              print *
           end if
           call die
        end if
     end do
     return
  end if

  if (var=='NCheb') then
     if (nvalues>1) goto 200
     type = 'INTEGER'
     read(value,*,ERR=100) NCheb
     return
  end if

  if (var=='rtmin') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) rtmin
     return
  end if

  if (var=='rtmax') then
     if (nvalues>1) goto 200
     type = 'REAL'
     read(value,*,ERR=100) rtmax
     return
  end if

  if (var=='TrackMinkowski') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) TrackMinkowski
     if (nvalues>1) goto 200
     return
  end if

  if (var=='TrackSchwarzschild') then
     if (nvalues>1) goto 200
     type = 'LOGICAL'
     read(value,*,ERR=100) TrackSchwarzschild
     if (nvalues>1) goto 200
     return
  end if

  print *
  if (rank==0) then
     print *, 'Parfile error, non-existent parameter in line:',nline
     print *
     print *, 'Aborting! (subroutine assign.f90)'
     print *
  end if
  call die
  100 continue
  if (rank==0) then
     print *
     print *, 'There was an error assigning the variable ''',trim(var),''' in line:',nline
     print *, 'Are you sure you gave a ',trim(type),' value?'
     print *
     print *, 'Aborting! (subroutine assign.f90)'
     print *
  end if
  call die

  200 continue
  if (rank==0) then
     print *
     print *, 'Multiple values not allowed for variable ''',trim(var),''' in line:',nline
     print *
     print *, 'Aborting! (subroutine assign.f90)'
     print *
  end if
  call die

  end subroutine assign
