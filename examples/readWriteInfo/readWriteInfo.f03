      program readWriteInfo
!
!     This program contains examples of the matrix file read and write routines  
!
      use mqc_gaussian
!
!
!     Variable Declarations...
!
      implicit none
      character(len=:),allocatable::fileName,newMatFile
 
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      integer::iOut=6,iPrint=2,nBasis,nAlpha,nBeta
      type(mqc_scf_integral)::mo_coefficients,overlap
      type(mqc_matrix)::sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor
!
      write(iOut,*)
      write(iOut,*) 'MQC Read-Write Examples'
      write(iOut,*)
      write(iOut,*) 'L. M. Thompson, 2017.'
      write(iOut,*)
!
!     Command line arguments can be returned as strings of correctly allocated length using 
!     mqc_get_command_argument
!
      call mqc_get_command_argument(1,fileName)
!
!     The header can be loaded from the file using the load procedure.
!
      call fileInfo%load(filename)
!
!     Using the getArray procedure, matrix file arrays can be read in exactly as found in the 
!     matrix file. 
!
      call fileInfo%getArray('SHELL TO ATOM MAP',sh2AtMp)
      call fileInfo%getArray('SHELL TYPES',shlTyp)
      call fileInfo%getArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
      call fileInfo%getArray('PRIMITIVE EXPONENTS',prmExp)
      call fileInfo%getArray('CONTRACTION COEFFICIENTS',conCoef)
      call fileInfo%getArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
      call fileInfo%getArray('COORDINATES OF EACH SHELL',shCoor)
!
!     The getESTObj procedure allows certain matrix file arrays to be simultaneously read and 
!     organized into spin blocks. Arrays are automatically be ordered regardless of whether the 
!     matrix file contains RHF, UHF, GHF and/or complex integrals.
!
      call fileInfo%getESTObj('overlap',est_integral=overlap)
      call fileInfo%getESTObj('mo coefficients',est_integral=mo_coefficients)
!  
!     EST objects can be printed in exactly the same way as MQC objects. The output will be 
!     automatically adjusted depending on the wavefunction/density type.
! 
      call overlap%print(iOut,'Overlap')
      call mo_coefficients%print(iOut,'MO coefficients')
!
!     The getVal procedure allows certain scalar values to be retrieved from the matrix file.
!
      nBasis = fileInfo%getVal('nBasis')
      nAlpha = fileInfo%getVal('nAlpha')
      nBeta = fileInfo%getVal('nBeta')
!
!     EST objects can be multiplied and transposed with each other, as well as with MQC variables.
!
      call mqc_print_integral(matmul(transpose(mo_coefficients),matmul(overlap,mo_coefficients)),iOut,'MO overlap')
!
!     A new matrix file can be created with the header information stored in fileinfo.
!
      call mqc_get_command_argument(2,newMatFile)
      call fileinfo%create(newMatFile)
!
!     Matrix file arrays can be output to the new matrix file using the writeArray procedure. 
!
      call fileInfo%writeArray('SHELL TO ATOM MAP',sh2AtMp)
      call fileInfo%writeArray('SHELL TYPES',shlTyp)
      call fileInfo%writeArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
      call fileInfo%writeArray('PRIMITIVE EXPONENTS',prmExp)
      call fileInfo%writeArray('CONTRACTION COEFFICIENTS',conCoef)
      call fileInfo%writeArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
      call fileInfo%writeArray('COORDINATES OF EACH SHELL',shCoor)
!
!     EST objects can be written using the writeESTObj procedure. This is automatically written 
!     to the matrix file depending on the spin blocks present, but can be specified manually with
!     the optional override argument.
!
      call fileinfo%writeESTObj('density',est_integral=(matmul(mo_coefficients,transpose(mo_coefficients))), &
        override='spin')
      call fileinfo%writeESTObj('mo coefficients',est_integral=mo_coefficients,override='general')
!
 999  end program readWriteInfo
