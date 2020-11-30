      program diishf
!
!     This program perfoms a DIIS Hartree-Fock calculation.
!
!     L. M. Thompson, 2018
!
      use mqc_gaussian
      use iso_fortran_env

      implicit none
      character(len=:),allocatable::command,fileName
      type(mqc_gaussian_unformatted_matrix_file)::matFile
      integer(kind=int64)::iOut=output_unit,iIn=input_unit,iPrint=0,i,j,nElec,nAlpha,nBeta,nBasis, &
        iter=0,max_iter=256,nBasUse,multi,nStor=20
      type(mqc_scalar)::conver,thresh
      type(mqc_wavefunction)::wavefunc
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs)::eris
      type(mqc_scalar)::Vnn,Energy,half
      type(mqc_vector)::tmp_vector
      type(mqc_matrix)::residual_matrix
      type(mqc_scf_integral)::Gmat,Ext_Fock,Xmat,old_density
      type(mqc_scf_integral),dimension(:),allocatable::error_matrix,Fock
!
      Write(IOut,*)
      Write(IOut,*) 'Hartree-Fock Density Matrix and Energy Calculator'
      Write(IOut,*)
      Write(IOut,*) 'L. M. Thompson 2018'
      Write(IOut,*)
      half = 0.5
      thresh = 1.0e-8
!
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f') then
!
!*      -f matrix_file                   Input matrix file with initial set of molecular orbitals.
!*
          call mqc_get_command_argument(i+1,fileName)
          j = i+2
        else
          call mqc_error_A('Unrecognised input flag',6,'command',command)
        endIf
        deallocate(command)
      endDo

      call matFile%load(fileName)
      call matFile%getESTObj('wavefunction',wavefunc)
      call wavefunc%print(iOut,'all')
      call matFile%getMolData(moleculeInfo)
      call matFile%get2ERIs('regular',eris)
      call eris%print(iOut,'AO 2ERIs')
!
!     Compute the nuclear-nuclear repulsion energy.
!
      Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
      call moleculeInfo%print(iOut)
      call Vnn%print(iOut,'Nuclear Repulsion Energy (au)')
!
      nElec = wavefunc%nElectrons%ival()
      nAlpha = wavefunc%nAlpha%ival()
      nBeta = wavefunc%nBeta%ival()
      nBasis = wavefunc%nbasis%ival()
      multi = wavefunc%multiplicity%ival()
      allocate(error_matrix(nStor),Fock(nStor))
!
!     initialize density matrix
      call wavefunc%density_matrix%identity(nBasis,nBasis)
      call wavefunc%density_matrix%print(iOut,'Initial density matrix')
!
!     determine orthogonalization matrix
      call mqc_scf_transformation_matrix(wavefunc%overlap_matrix,Xmat,nBasUse)
      call Xmat%print(iOut,'Orthogonalization matrix')
      if(nBasUse.ne.nBasis) call mqc_error('Linear dependencies not implemented')
!
!     enter iterations
      do while (iter+1.le.max_iter)
        write(iOut,'(A,I3)') 'Iteration: ',iter+1
!
!       form G matrix
        Gmat = contraction([eris],wavefunc%density_matrix)
        call Gmat%print(iOut,'G(P)')
!
!       form Fock matrix and save
        Fock(mod(iter,nStor)+1) = wavefunc%core_hamiltonian + Gmat
        call Fock(mod(iter,nStor)+1)%print(iOut,'Fock matrix')
!
!       compute residual and add to residual list
        error_matrix(mod(iter,nStor)+1) = matmul(matmul(Fock(mod(iter,nStor)+1),wavefunc%density_matrix),wavefunc%overlap_matrix) & 
          - matmul(wavefunc%overlap_matrix,matmul(wavefunc%density_matrix,Fock(mod(iter,nStor)+1))) 
!
!       construct system of equations
        call residual_matrix%init(min(nStor+1,iter+2),min(nStor+1,iter+2))
        tmp_vector = [(-1,i=1,min(nStor,iter+1)),0]
        call residual_matrix%vput(tmp_vector,[0],[-1])
        call residual_matrix%vput(tmp_vector,[-1],[0])
        do i = 1, min(nStor,iter+1)
          do j = 1, min(nStor,iter+1)
            call residual_matrix%put(mqc_scf_integral_trace(matmul(error_matrix(i),error_matrix(j))),i,j)
          endDo
        endDo
        call residual_matrix%print(iOut,'DIIS Residual_matrix')
!
!       solve for coefficients
        tmp_vector = [(0,i=1,min(nStor,iter+1)),-1]
        tmp_vector = matmul(residual_matrix%inv(),tmp_vector)
        call tmp_vector%print(iout,'Fock coefficient vector')
!
!       construct new Fock matrix 
        call Ext_Fock%init(nBasis,nBasis)
        do i = 1, min(nStor,iter+1)
          Ext_Fock = Ext_Fock + tmp_vector%at(i)*Fock(i)
        endDo
        call Ext_Fock%print(iout,'Extrapolated Fock matrix')
!
!       orthogonalize extrapolated Fock basis
        Ext_Fock = matmul(dagger(Xmat),matmul(Ext_Fock,Xmat))
        call Ext_Fock%print(iOut,'Orthogonalized Extrapolated Fock matrix')
!
!       diagonalize extrapolated Fock matrix
        call Ext_Fock%diag(wavefunc%MO_energies,wavefunc%MO_coefficients)
!
!       back-transform MO coefficients
        call wavefunc%MO_coefficients%print(iOut,' Orthogonal MO coefficients')
        wavefunc%MO_coefficients = matmul(Xmat,wavefunc%MO_coefficients)
        if(wavefunc%wf_type.eq.'G') then
          wavefunc%MO_coefficients = wavefunc%MO_coefficients%orbitals(alphaOrbsIn=[(i,i=1,nBasis,2),(i,i=-2,-1*nBasis,-2)], &
            betaOrbsIn=[(i,i=-2,-1*nBasis,-2),(i,i=1,nBasis,2)])
        endIf
        call wavefunc%MO_energies%print(iOut,'Orbital energies')
        call wavefunc%MO_coefficients%print(iOut,'MO coefficients')
!        
!       form density matrix 
        old_density = wavefunc%density_matrix
        wavefunc%density_matrix = matmul(wavefunc%MO_coefficients%orbitals('occupied',[nAlpha],[nBeta]), &
          dagger(wavefunc%MO_coefficients%orbitals('occupied',[nAlpha],[nBeta])))
        call wavefunc%density_matrix%print(iOut,'Density matrix')
!
!       compute energy at iteration
        Energy = contraction(wavefunc%core_hamiltonian,wavefunc%density_matrix)
        Energy = Energy + half*contraction(wavefunc%density_matrix,Gmat)
        Energy = Energy + Vnn
        call Energy%print(iOut,'Energy',FormatStr='F20.8')
!
!       test convergence and exit or iterate
        conver = mqc_integral_norm((wavefunc%density_matrix-old_density),'F')
        call conver%print(iOut,'Convergence on density matrix',FormatStr='F20.8')
        if(conver.le.thresh) then
          write(iOut,'(1x,A,I5,1x,A)') 'Convergence achieved in',iter,'iterations'
          exit
        endIf
!        
        iter = iter+1
      endDo
!
      end program diishf
