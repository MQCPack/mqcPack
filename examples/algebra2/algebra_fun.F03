      program algebra_fun
      use MQC_Algebra2
      use iso_fortran_env
!
!     This program tests the key functionality of MQC module MQC_Algebra2. 
!
!
      implicit none
      type(MQC_Variable)::var1,var2,vec1,vec2,vec3,mat1,mat2
      integer(kind=int64)::intrinsicInt1
      integer(kind=int64),dimension(5)::aaa=[1,2,3,4,5]
      integer(kind=int64)::iout=6
      real(kind=real64)::intrinsicReal1
!
!
!     Test scalars...
!
      write(Iout,1010)
      write(Iout,1000)' Testing scalars...'
!
!     Test initialization and reporting functions/routines.
      write(Iout,1000)'    Testing initialization and reporting.'
! Need format statements.  Different compilers produce slightly different
! Output when using * format
 1000 format( A )
 1010 format( ' ' )
 1020 format( A, G21.15 )
 1030 format( A, I15 )
 1040 format( 2A )
 1050 format( 3A )
      call MQC_Variable_initialize(var1,'real')
      call MQC_Variable_initialize(var2,'int')
      write(Iout,1000)'      Details for var1:'
      write(Iout,1040)'        Type = ',TRIM(MQC_Variable_getType(var1))
      write(Iout,1030)'        Rank = ',MQC_Variable_getRank(var1)
      write(Iout,1030)'        Len  = ',MQC_Variable_getLength(var1)
      write(Iout,1000)'      Details for var2:'
      write(Iout,1040)'        Type = ',TRIM(MQC_Variable_getType(var2))
      write(Iout,1030)'        Rank = ',MQC_Variable_getRank(var2)
      write(Iout,1030)'        Len  = ',MQC_Variable_getLength(var2)
      write(Iout,1010)
      write(Iout,1000)' Hrant - Printing var1 and var2.'
      call mqc_variable_fillVal(var1,999)
      var1 = 678
      call mqc_variable_fillVal(var2,8.88)
      var2 = 1.23
      call mqc_variable_fillVal(var2,1024.5)

      call MQC_Print_mqcVariable(var1,header='var1')
      call MQC_Print_mqcVariable(var2,header='var2')

      intrinsicInt1 = INT(var1)
      write(Iout,1030)' intrinsicInt1(var1) = ',intrinsicInt1

      intrinsicInt1 = INT(var2)
      write(Iout,1030)' intrinsicInt1(var2) = ',intrinsicInt1

      intrinsicInt1 = var1
      write(Iout,1030)' intrinsicInt1(var1), using assignment, = ',intrinsicInt1

      intrinsicReal1 = FLOAT(var1)
      write(Iout,1020)' intrinsicReal1(var1) = ',intrinsicReal1

      intrinsicReal1 = FLOAT(var2)
      write(Iout,1020)' intrinsicReal1(var2) = ',intrinsicReal1

      intrinsicInt1 = var1
      write(Iout,1020)' intrinsicReal1(var1), using assignment, = ',intrinsicReal1
!
!
!     Test vectors...
!
      write(Iout,1010)
      write(Iout,1000)' Testing vectors...'
!
!     Test initialization and reporting functions/routines.
      write(Iout,1000)'    Testing initialization and reporting.'
      call MQC_Variable_initialize(vec1,'int',rank=1,dimensions=[5])
      call MQC_Variable_initialize(vec2,'real',rank=1,dimensions=[32])
      write(Iout,1000)'      Details for vec1:'
      write(Iout,1040)'        Type    = ',TRIM(MQC_Variable_getType(vec1))
      write(Iout,1030)'        Rank    = ',MQC_Variable_getRank(vec1)
      write(Iout,1030)'        Size(1) = ',MQC_Variable_getSize(vec1,1)
      write(Iout,1030)'        Len     = ',MQC_Variable_getLength(vec1)
      write(Iout,1000)'      Details for vec2:'
      write(Iout,1040)'        Type    = ',TRIM(MQC_Variable_getType(vec2))
      write(Iout,1030)'        Rank    = ',MQC_Variable_getRank(vec2)
      write(Iout,1030)'        Size(1) = ',MQC_Variable_getSize(vec2,1)
      write(Iout,1030)'        Len     = ',MQC_Variable_getLength(vec2)

      call mqc_variable_fillVal(vec1,[1,2,3,4,5])
      write(Iout,1010)
      write(Iout,1000)' About to attempt vector printing...'
      call MQC_Print_mqcVariable(vec1,header='vec1')

      write(Iout,1010)
      write(Iout,1000)' Resetting vec1 and printing again...'
      vec1 = [6,4,2,0,-2,-4]
      call vec1%print()
      goto 999

      write(Iout,1010)
      write(Iout,1000)' Resetting vec1 as a real and printing again...'
#ifdef PGI
      vec1 = dfloat([6,4,2,0,-2,-4])
#else
      vec1 = float([6,4,2,0,-2,-4])
#endif
      call vec1%print(header='vec1:')

      vec3 = [1.3,2.4,3.5,3.2,-1.43,-12.0]
      call vec3%print(header='vec3=')

      write(Iout,1010)
      write(Iout,1020)' vec1.conformable.vec2 = ',vec1%isConformable(vec2)
      write(Iout,1020)' vec1.conformable.vec3 = ',vec1%isConformable(vec3)
      write(Iout,1020)' vec2.conformable.vec3 = ',vec2%isConformable(vec3)
      write(Iout,1010)
      write(Iout,1020)' vec1.conformable.vec2 = ',vec1.conformable.vec2
      write(Iout,1020)' vec1.conformable.vec3 = ',vec1.conformable.vec3
      write(Iout,1020)' vec2.conformable.vec3 = ',vec2.conformable.vec3
      write(Iout,1010)

      vec1 = [6,4,2,0,-2,-4]
      vec3 = [16.5,14.,12.,10.,-12.,-14.]
      call vec1%print(header='vec1')
      call vec3%print(header='vec3')
      vec2 = MQC_Variable_Addition(vec1,vec3)
      call mqc_print(vec2,header='vec1+vec3')
      call mqc_print(MQC_Variable_Addition(vec1,vec3),header='vec1+vec3 ... AGAIN:')
      call mqc_print(vec1+vec3,header='vec1+vec3 ... YET AGAIN:')

      call mqc_print(vec1-vec3,header='vec1-vec3:')

      call mqc_print(vec1*vec3,header='vec1*vec3:')

      call mqc_print(vec1/vec3,header='vec1/vec3:')

!
!     Test dot_product.
!
      var1 = mqc_variable_contraction_full(vec1,vec3)
      call mqc_print(var1,header='vec1.vec3')
      call mqc_print(mqc_variable_contraction_full(vec1,vec3),header='vec1.vec3')
!
!
!     Test matrices...
!
      write(Iout,1010)
      write(Iout,1000)' Testing matrices...'
!
!     Test initialization and reporting functions/routines.
      write(Iout,1000)'    Testing initialization and reporting.'
      call MQC_Variable_initialize(mat1,'real',dimensions=[5,2])
      call MQC_Variable_initialize(mat2,'int',rank=2,dimensions=[10,32])
      write(Iout,1000)'      Details for mat1:'
      write(Iout,1040)'        Type    = ',TRIM(MQC_Variable_getType(mat1))
      write(Iout,1030)'        Rank    = ',MQC_Variable_getRank(mat1)
      write(Iout,1030)'        Size(1) = ',MQC_Variable_getSize(mat1,1)
      write(Iout,1030)'        Size(2) = ',MQC_Variable_getSize(mat1,2)
      write(Iout,1030)'        Len     = ',MQC_Variable_getLength(mat1)
      write(Iout,1010)
      write(Iout,1000)'      ** AGAIN USING INTRINSIC INTERFACES...'
      write(Iout,1030)'        Rank    = ',RANK(mat1)
      write(Iout,1030)'        Size(1) = ',SIZE(mat1,1)
      write(Iout,1030)'        Size(2) = ',SIZE(mat1,2)
      write(Iout,1030)'        Len     = ',SIZE(mat1)
      write(Iout,1010)

      write(Iout,1000)'      Details for mat2:'
      write(Iout,1040)'        Type    = ',TRIM(MQC_Variable_getType(mat2))
      write(Iout,1030)'        Rank    = ',MQC_Variable_getRank(mat2)
      write(Iout,1030)'        Size(1) = ',MQC_Variable_getSize(mat2,1)
      write(Iout,1030)'        Size(2) = ',MQC_Variable_getSize(mat2,2)
      write(Iout,1030)'        Len     = ',MQC_Variable_getLength(mat2)


!
!     Test out printing...
!
      call mqc_print_Vector_Array_Integer(aaa,6)
!
!     Test int2char routine...
!
      write(Iout,1010)
      write(Iout,1010)
      write(Iout,1000)' Hrant - This is a test of the int  -> char function...'
      write(Iout,1050)'    here is the number 12345 --->',TRIM(integer2character(12345)),'<-----'
      write(Iout,1010)
      write(Iout,1010)
      write(Iout,1010)
      write(Iout,1000)' Hrant - This is a test of the real -> char function...'
      write(Iout,1050)'    here is the number 1.2345 --->',TRIM(real2character(1.2345e-6,'d50.15')),'<-----'
      write(Iout,1010)
      call MQC_Print_Scalar_Integer(98765,6,'This is my INTEGER number:',.true., &
        .true.)
      call MQC_Print_Scalar_Real(9.8765,6,'This is my REAL number:',.true., &
        .true.)
      
      
!
 999  continue
      end program algebra_fun
