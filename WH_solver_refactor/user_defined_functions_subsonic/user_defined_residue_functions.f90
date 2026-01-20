MODULE user_defined_residue_functions

  USE input_params
  USE bessel_utils
  USE io_utils
  USE user_defined_functions

  IMPLICIT NONE

  PRIVATE ::  residue_pr_instab_hard_duct_mode, &
              residue_pr_instab_guided_jet_mode,&
              compute_F_GJ, compute_d_ds_Trs_lambda

  PUBLIC  ::  residue_pr_instab, residue_pot_instab, &
              residue_k_plus_duct

  CONTAINS

  FUNCTION residue_pr_instab(r,z,input_data)  result(residuepr_op)
    real(dpk)            :: r, z
    integer              :: ii, jj
    complex(dpk)         :: residuepr_op, res, fn
    type(input_params_t) :: input_data

    if (input_data%solution_mode == 'guided_jet') then
       residuepr_op =  residue_pr_instab_guided_jet_mode(r,z,input_data)
    else
       residuepr_op =  residue_pr_instab_hard_duct_mode(r,z,input_data)
    end if

  END FUNCTION residue_pr_instab

  FUNCTION residue_pr_instab_hard_duct_mode(r,z,input_data) result(residuepr_op)
    real(dpk)            :: r, z
    integer              :: ii, jj
    complex(dpk)         :: residuepr_op, res, fn
    type(input_params_t) :: input_data


    if (input_data%vs_param_gamma .EQ. 0) then
       residuepr_op = 0.
    else 
      res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M1)/((input_data%mu_plus - input_data%KH_zero_1)* & 
                                                               input_data%k_minus_at_mu_plus*input_data%k_plus_sz1)

      if (r <= 1) then
         fn = ((1._dpk -(input_data%KH_zero_1*input_data%M1))**2)*compute_Trs_lambda(r,input_data%KH_zero_1,1,input_data)
      else
         fn = ((1._dpk -(input_data%KH_zero_1*input_data%M2))**2)*compute_Trs_lambda(r,input_data%KH_zero_1,2,input_data)
      end if
       
      residuepr_op = ((input_data%omega_r)**2)*res*fn*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%KH_zero_1*z)
      
    end if      
    
  END FUNCTION residue_pr_instab_hard_duct_mode

  FUNCTION residue_k_plus_duct(r,z,input_data) result(residuepr_op)

   real(dpk),    intent(in) :: r, z
   type(input_params_t), intent(in) :: input_data
   complex(dpk) :: residuepr_op
   complex(dpk) :: res, fn, fn_num, fn_den

   if ( z < 0 ) then
      residuepr_op = 0._dpk
      return
   end if 

   if ( r .LE. 1._dpk) then
      fn_num =  ((1._dpk-(input_data%k_d_plus*input_data%M1))**2)
      fn_num =  fn_num*((1._dpk-(input_data%mu_plus*input_data%M1)))
      fn_num =  fn_num*compute_Trs_lambda(r,input_data%k_d_plus,1,input_data)
     

      fn_den =  (((input_data%mu_plus - input_data%KH_zero_1)*input_data%k_minus_at_mu_plus*&
                input_data%k_plus_at_k_d_plus*(input_data%mu_plus-input_data%k_d_plus  ))) 


      fn = fn_num/fn_den

      residuepr_op = ((input_data%omega_r)**2)*fn*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%k_d_plus*z)
   else  
      residuepr_op = 0._dpk
   end if


  END FUNCTION residue_k_plus_duct

  FUNCTION residue_pr_instab_guided_jet_mode(r,z,input_data) result(residuepr_op)

    real(dpk)            :: r, z
    integer              :: ii, jj, stream_idx
    complex(dpk)         :: residuepr_op, pre_factor, fn
    type(input_params_t) :: input_data

    if (r .LE. 1) then
       stream_idx = 1
    else
       stream_idx = 2
    end if

    pre_factor = (CMPLX(0._dpk,1._dpk,kind=dpk)) *&
                 ((input_data%omega_r)**2)*input_data%C0*(1 - (input_data%mu_plus*input_data%M1))
 
    pre_factor = pre_factor*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r * input_data%KH_zero_1 * z)
  
    pre_factor = pre_factor/&
                (input_data%k_minus_at_mu_plus*input_data%k_plus_sz1*((input_data%KH_zero_1-input_data%mu_plus)**2) )


    fn = compute_Trs_lambda(r,input_data%KH_zero_1,stream_idx,input_data)

    if (stream_idx == 1) then
        residuepr_op = pre_factor * fn * (1._dpk - (input_data%KH_zero_1 * input_data%M1))**2
    else if (stream_idx == 2) then
        residuepr_op = pre_factor * fn * (1._dpk - (input_data%KH_zero_1 * input_data%M2))**2
    end if 
    
  END FUNCTION residue_pr_instab_guided_jet_mode

  FUNCTION residue_pr_GJ_guided_jet_mode(r,z,input_data) result(residuepr_op)

    real(dpk)            :: r, z
    integer              :: ii, jj, stream_idx
    complex(dpk)         :: residuepr_op, pre_factor, fn
    type(input_params_t) :: input_data

    if (r .LE. 1) then
       stream_idx = 1
    else
       stream_idx = 2
    end if

    if ( r .LE. 1 .AND. z .LE. 0) then
        residuepr_op = 0._dpk
        return
    end if


    pre_factor = (CMPLX(0._dpk,1._dpk,kind=dpk)) *&
                 ((input_data%omega_r)**2)*input_data%C0*(1 - (input_data%mu_plus*input_data%M1))
 
    pre_factor = pre_factor*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r * input_data%mu_plus * z)
  
    pre_factor = pre_factor/&
                (input_data%k_minus_at_mu_plus*input_data%k_plus_at_mu_plus*((input_data%mu_plus-input_data%KH_zero_1)) )


    fn = compute_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)

    if (stream_idx == 1) then
        residuepr_op = pre_factor * fn * (1._dpk - (input_data%mu_plus * input_data%M1))**2
    else if (stream_idx == 2) then
        residuepr_op = pre_factor * fn * (1._dpk - (input_data%mu_plus * input_data%M2))**2
    end if 
 
  END FUNCTION residue_pr_GJ_guided_jet_mode


!  FUNCTION residue_pr_GJ_guided_jet_mode(r,z,input_data) result(residue_pr_GJ_op)
!
!    real(dpk)            :: r, z
!    integer              :: ii, jj, stream_idx
!    complex(dpk)         :: residue_pr_GJ_op, pre_factor, fn
!    complex(dpk)         :: res1, res2
!    type(input_params_t) :: input_data
!
!    if (r .LE. 1) then
!       stream_idx = 1
!    else
!       stream_idx = 2
!    end if
!
!    pre_factor = CMPLX(0._dpk,1._dpk,kind=dpk)
!    pre_factor = pre_factor*((input_data%omega_r)**2) * input_data%C0  
!    pre_factor = pre_factor/&
!                (input_data%k_minus_at_mu_plus )
!    pre_factor = pre_factor*(1- (input_data%s_GJ*input_data%M1))      
!
!    if (stream_idx == 1) then
!       res1 = (1- (input_data%s_GJ*input_data%M1))*compute_d_ds_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)
!       res2 = input_data%M1*compute_F_GJ(stream_idx,r,z,input_data)
!    elseif (stream_idx == 2) then
!       res1 = (1- (input_data%s_GJ*input_data%M2))*compute_d_ds_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)
!       res2 = input_data%M2*compute_F_GJ(stream_idx,r,z,input_data)
!    end if 
!
!    residue_pr_GJ_op = pre_factor * (res1 - res2)
!    
!  END FUNCTION residue_pr_GJ_guided_jet_mode


  FUNCTION compute_F_GJ(stream_idx,r,z,input_data) result(F_value)

    integer              :: stream_idx
    real(dpk)            :: r, z, mach_no
    type(input_params_t) :: input_data
    complex(dpk)         :: F_value


    if (stream_idx == 1) then
       mach_no = input_data%M1
    else if (stream_idx == 2) then
       mach_no = input_data%M2
    end if
  
    F_value = (1 - input_data%mu_plus*mach_no)
    F_value = F_value*compute_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)
    F_value = F_value*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r *input_data%mu_plus * z)
    F_value = F_value/(input_data%k_plus_at_mu_plus * (input_data%mu_plus - input_data%KH_zero_1))


  END FUNCTION compute_F_GJ

  FUNCTION compute_d_ds_F_GJ(stream_idx,r,z,input_data) result(dF_ds_value)

    integer              :: stream_idx
    real(dpk)            :: r, z, mach_no
    type(input_params_t) :: input_data
    complex(dpk)         :: dF_ds_value, pre_factor
    complex(dpk)         :: term_1, term_2, term_3, term_4,&
                            term_4_1_num, term_4_1_den
  
    if (stream_idx == 1) then
       mach_no = input_data%M1
    else if (stream_idx == 2) then
       mach_no = input_data%M2
    end if

    pre_factor = EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus * z)
    pre_factor = pre_factor/(input_data%k_plus_at_mu_plus * (input_data%mu_plus -input_data%KH_zero_1))


    term_1 = -mach_no*compute_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)

    term_2 = (1 -(input_data%mu_plus*mach_no))*compute_d_ds_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)

    term_3 = (1 - (input_data%mu_plus*mach_no))*compute_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)
    term_3 = term_3*CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r * z

    term_4 = -(1- (input_data%mu_plus*mach_no))
    term_4 = term_4*compute_Trs_lambda(r,input_data%mu_plus,stream_idx,input_data)

    term_4_1_num = input_data%k_plus_prime_at_mu_plus*(input_data%mu_plus - input_data%KH_zero_1)
    term_4_1_num = term_4_1_num + input_data%k_plus_at_mu_plus
    term_4_1_den = input_data%k_plus_at_mu_plus*(input_data%mu_plus - input_data%KH_zero_1)

    term_4 = term_4*term_4_1_num/term_4_1_den

    dF_ds_value = pre_factor*(term_1 + term_2 + term_3 + term_4)

  END FUNCTION compute_d_ds_F_GJ

!
  FUNCTION residue_pot_instab(r,z,ss,input_data) result(residue_pot_instab_op)
 !
 !!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
 !!1. Same as above but for computing velocity potential
 !!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
 !
    real(dpk)            :: r, z
    integer              :: ss, ii, jj
    complex(dpk)         :: residue_pot_instab_op, res, fn
   
    type(input_params_t) :: input_data
 !
 !
    if (input_data%vs_param_gamma .EQ. 0) then
 !
       residue_pot_instab_op = 0.
 !
    else 
 !
        res = input_data%psi*(1._dpk - input_data%mu_plus*input_data%M1)/((input_data%mu_plus - input_data%KH_zero_1)* & 
                                                                 input_data%k_minus_at_mu_plus*input_data%k_plus_sz1)
 !
         
        if (r <= 1) then
         fn = (1._dpk - (input_data%KH_zero_1*input_data%M1))*compute_Trs_lambda(r,input_data%KH_zero_1,1,input_data)
        else
         fn = (1._dpk - (input_data%KH_zero_1*input_data%M2))*compute_Trs_lambda(r,input_data%KH_zero_1,2,input_data)
        end if
         
        residue_pot_instab_op = input_data%omega_r*res*fn*EXP(CMPLX(0.,1._dpk,kind=dpk)*input_data%omega_r*input_data%KH_zero_1*z)
        
    end if      
    
  END FUNCTION residue_pot_instab

END MODULE user_defined_residue_functions

