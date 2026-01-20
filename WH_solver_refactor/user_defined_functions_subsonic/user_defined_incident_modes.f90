Module user_defined_incident_modes

  USE input_params
  USE bessel_utils
  USE io_utils  

  IMPLICIT NONE

  PRIVATE :: compute_psi_incident_hard_duct_mode,&
             compute_psi_incident_guided_jet_mode

  PUBLIC  :: compute_psi_incident

  CONTAINS 

  FUNCTION compute_psi_incident(r,z,input_data) result(psi0)

!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!
!! 1. Compute the \Psi_{mn}(r) of (2.14) or (4.11)
!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!=!!

    complex(dpk)        :: psi0, psimn, resp
    real(dpk)           :: r, z

    type(input_params_t) :: input_data

    if (input_data%solution_mode == 'guided_jet') then
       psi0 = compute_psi_incident_guided_jet_mode(r,z,input_data)
    else
       psi0 = compute_psi_incident_hard_duct_mode(r,z,input_data)
    end if

  END FUNCTION compute_psi_incident


 FUNCTION compute_psi_incident_hard_duct_mode(r,z,input_data) result(psi0)

    complex(dpk)        :: psi0, psimn, resp
    complex(dpk)        :: f1, f2, f3
    real(dpk)           :: r, z
    integer             :: stream_idx

    type(input_params_t) :: input_data

    if (r .LE. 1) then 
       stream_idx = 1
    else
       stream_idx = 2
    end if

    if (stream_idx == 1) then
         psimn = bessj(input_data%alpha1*r,input_data%azim_mode,1)*EXP(ABS(AIMAG(input_data%alpha1*r)))
       
         psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*(1._dpk - input_data%M1*input_data%mu_plus)* &
              psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*z)

     else
          psi0 = 0.
    end if

  END FUNCTION compute_psi_incident_hard_duct_mode

 FUNCTION compute_psi_incident_guided_jet_mode(r,z,input_data) result(psi0)

    complex(dpk)        :: psi0, psimn, resp
    complex(dpk)        :: f1, f2, f3
    real(dpk)           :: r, z
    integer             :: stream_idx

    type(input_params_t) :: input_data

 
    if (r .LE. 1) then 
       stream_idx = 1
    else
       stream_idx = 2
    end if

    if (stream_idx == 1) then
         psimn = bessj(input_data%alpha1*r,input_data%azim_mode,1)*EXP(ABS(AIMAG(input_data%alpha1*r)))
       
         psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*(1._dpk - input_data%M1*input_data%mu_plus)* &
              psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*z)
    else

          f1 = ((1._dpk - (input_data%mu_plus*input_data%M1) )/(1._dpk -(input_data%mu_plus*input_data%M2))) &
                                                                 *bessj(input_data%alpha1,input_data%azim_mode,1)
          f2 =  bessj(input_data%alpha2*r,input_data%azim_mode,1)
          f3 =  hank1(input_data%alpha2,input_data%azim_mode,1)      

          psimn = (input_data%kap_rho)*f1*f2/f3

          psi0 = CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*(1._dpk -input_data%M2*input_data%mu_plus)* &
                 psimn*EXP(CMPLX(0._dpk,1._dpk,kind=dpk)*input_data%omega_r*input_data%mu_plus*z)

    !      psi0 = 0.0

    end if
  
 END FUNCTION compute_psi_incident_guided_jet_mode

END MODULE user_defined_incident_modes
