# To Do List for MQC
This is a list of Hrant's ongoing to-do list for MQC -- focusing primarily on the mqc_algebra2 module. There are two sub-lists. The first sub-list shows developments that are planned for implementation; the second list shows the to-do items that have been fully implemented and relatively well-validated/tested. In many cases, code snippets are included to demonstrated how one can use MQC to carry out the listed implemented features.

## MQC_Algebra2

#### TO DO

* Matrix-Matrix multiply
* Matrix-Vector multiply
* Matrix diagonalization
* Matrix inversion

#### DONE

* Matrix transpose (for full storage matrices only)
  ```
  type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
  type(MQC_Variable)::moCoefficientsAlpha

  call fileInfo%load('test.mat')
  call fileInfo%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=moCoefficientsAlpha)
  call mqc_print_mqcVariable(Transpose(moCoefficientsAlpha),header='Transpose')
  ```



* Convert from intrinsic (INT and FLOAT) $\rightarrow$ mqc_var type
  ```
  type(MQC_Variable)::var1

  vec1 = [1.3,2.4,3.5,3.2,-1.43,-12.0]
  call vec1%print(header='vec1=')
  ```

* Convert from mqc_var type $\rightarrow$ intrinsic (INT and FLOAT)
  ```
  type(MQC_Variable)::overlapMatrixAO,densityMatrixAlpha

  write(*,*)' <PS> = ',float(2)*Float(contraction(densityMatrixAlpha,overlapMatrixAO))
  ```

* (Full) contraction
  ```
  type(MQC_Variable)::overlapMatrixAO,densityMatrixAlpha

  write(*,*)' <PS> = ',float(2)*Float(contraction(densityMatrixAlpha,overlapMatrixAO))
  ```

* dot_product
  ```
  type(MQC_Variable)::var1,var2,vec1,vec2,vec3,mat1,mat2

  vec1 = float([6,4,2,0,-2,-4])
  call vec1%print(header='vec1:')

  var1 = mqc_variable_contraction_full(vec1,vec3)
  call mqc_print(var1,header='vec1.vec3')
  ```
