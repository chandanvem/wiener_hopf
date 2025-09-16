module input_params 
     
     IMPLICIT none

     integer, parameter :: dpk = kind(1.0d0)
   
     type :: input_params_t
 
          real(dpk)                                  :: tol  !! specified tolerance level
          integer                                    :: num_ker_pts_loop !! number of starting pts for kernel contour
          integer                                    :: total_ker_points, tot_IFT_pts
          integer                                    :: prswitch, restart, reflswitch, farswitch, vortswitch
          integer                                    :: num_zeros_s1_s2, num_poles_s1_s2  !! num  zeros & poles in between s_1- s_2
          integer                                    :: num_sup_zeros, num_sup_poles  !! number of supersonic zeros & poles
          integer                                    :: Nphi  !! polar mesh resolution for directivity computation
          complex(dpk), allocatable, dimension(:)    :: zeros_list_bw_s1_s2, poles_list_bw_s1_s2, sup_zeros_list, sup_poles_list
          complex(dpk), allocatable, dimension(:)    :: kzsp
          integer                                    :: num_of_streams, St_flag
          real(dpk)                                  :: M1, M2, M3
          real(dpk)                                  :: h
          real(dpk)                                  :: kapT  !! sqrt(T1/T0)
          real(dpk)                                  :: kap_rho  !! rho1/rho0
          real(dpk)                                  :: theta,  offset
          real(dpk)                                  :: omega_st, w0
          real(dpk)                                  :: Rmax, Rmin, Zmin, Zmax !! dimensions of the physical domain
          real(dpk), allocatable, dimension(:)       :: R, Z
          integer                                    :: Nmeshr, Nmeshz !! 400X400 -> Seg Fault!!
          real(dpk)                                  :: azim_mode
          real(dpk)                                  :: Zo  !! start point of incident vorticity (negative)
          real(dpk)                                  :: vs_param_gamma  !! Vortex shedding parameter
          complex(dpk)                               :: KH_zero_1, KH_zero_2, KH_pole_1  !! instability zero 1
          complex(dpk)                               :: mu_plus !! axial wave no (of incident wave)
          complex(dpk)                               :: mu0 !! incident axial wave no (for incident vort)
          complex(dpk)                               :: psi
          integer                                    :: num_IFT_pts_loop !! number of points for IFT contour
          real(dpk)                                  :: asymplim, asymplim1 !! asymptotic limit for the kernel functions
          real(dpk)                                  :: cnangle 
          real(dpk), allocatable, dimension(:)       :: cnanglei
          complex(dpk)                               :: omega_r
          complex(dpk)                               :: k_minus_at_mu_plus, k_plus_sz1,k_plus_sz2, kpsp1
          real(dpk), allocatable, dimension(:)       :: phi  
          complex(dpk), allocatable, dimension(:)    :: fplusz
          complex(dpk)                               :: alpha1,alpha2

          complex(dpk)                               :: s_GJ
          integer                                    :: GJ_num_int_pts 
          character(len=40)                          :: solution_mode

     end type input_params_t


     type :: contour_params_t

          complex(dpk), allocatable, dimension(:)    :: def_pts_IFT_cntr  !! the points defining the IFT contour
          complex(dpk), allocatable, dimension(:)    :: def_pts_ker_cntr  !! the points defining the Kernel contour
          real(dpk)                                  :: cont_cross_over_pt 
          integer                                    :: total_ker_points,tot_IFT_pts 
          complex(dpk), allocatable, dimension(:)    :: ker_int_points  !!location of the starting pts
          complex(dpk), allocatable, dimension(:)    :: iftpoints  
          complex(dpk)                               :: GJ_ref_pt
          complex(dpk)                               :: GJ_cntr_maxima
          integer(dpk)                               :: GJ_num_int_pts
     end type contour_params_t

end module input_params

