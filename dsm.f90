      !==============================================================!
      !              Dynamic Smagorinsky Model (DSM) Module          !
      !==============================================================!
      !            Balachandra R. Mettu <brmettu@ncsu.edu>           !
      !                    Last modified: 16 apr 2019                !
      !==============================================================!
      ! REFERENCE:
      !- Germano et al. (1991), Physics of Fluids - Original model
      !- Moin et al. (1991), Physics of Fluids - Compressible form
      !- Lily (1992) Physics of Fluids - Modification for stability 
      !- Park & Mahesh (2007), AIAA Forum - Filtering/averaging for unstructured grids
      !==============================================================!
      MODULE dsm
          ! PLACEHOLDER FOR IMPORTANT MODULES WHICH ARE NEEDED FOR DSM
          ! use MPI_module
          ! use connectivity_module
          ! use geometry_module
          ! use primitives_module
          ! use gradients_module 
          implicit none

          ! Upper bound for clipping Cs, Ci and Prt
          real(8), parameter :: dsm_clip_cs = 1.0d0 
          real(8), parameter :: dsm_clip_ci = 0.1d0
          real(8), parameter :: dsm_clip_pt = 1.0d1

          ! Arrays for Cs, Ci, Prt at each cell center
          real(8), allocatable :: cs_dsm(:),ci_dsm(:),pt_dsm(:)
          real(8), allocatable :: mcs(:),lcs(:)
          real(8), allocatable :: mci(:),lci(:)
          real(8), allocatable :: mpt(:),lpt(:)

          ! (test/grid filter)^2
          real(8), parameter :: ratf2 = 2.0d0**2 

      CONTAINS

          !-----------------------------------------------------------
          ! Handles the DSM routines
          !-----------------------------------------------------------
          subroutine handle_dsm(iact)
              implicit none
              integer, intent(in) :: iact

              integer :: i

              select case (iact)

              case(0)
                  call dsm_init() ! Initialize arrays
              case(1)
                  call dsm_fil()  ! Filtering for DSM
                  call dsm_avg()  ! Averaging for DSM
              case(2)

              case default
                  print*, "*** ERROR: INVALID 'IACT' SELECTED IN HANDLE_LES !!!"
                  call MPI_FINALIZE(ier)
                  stop
              end select

              return
          end subroutine handle_dsm


          !-----------------------------------------------------------
          ! Initialize DSM arrays
          !-----------------------------------------------------------
          subroutine dsm_init()
              implicit none

              !
              ! nel = # of interior cells
              ! net = # of interior + ghost + shared cells
              !

              allocate(cs_dsm(net))
              allocate(ci_dsm(net))
              allocate(pt_dsm(net))

              allocate(mcs(net),lcs(net))
              allocate(mci(net),lci(net))
              allocate(mpt(net),lpt(net))

              !
              !- Any other initalization routines can be placed here
              !

              return
          end subroutine dsm_init


          !-----------------------------------------------------------
          ! Filtering routine for DSM
          !-----------------------------------------------------------
          ! Outputs:- 
          !-----------
          ! lcs(1:net), mcs(1:net)
          ! lci(1:net), mci(1:net)
          ! lpt(1:net), mpt(1:net)
          !-----------------------------------------------------------
          subroutine dsm_fil()
              implicit none

              integer :: i,ii,iel,k,j
              real(8) :: rq,uq,vq,wq,tq,wgtf
              real(8) :: duij(3,3),dti(3),sdij(6),sskk,ssij
              real(8) :: rf,ruf(3),ruuf(6),sdf(6),ssf,rsssdf(6),rss2f
              real(8) :: dtf(3),rssdtf(3),rtf,rutf(3)
              real(8) :: lij(6),mij(6),lkk,lti(3),mti(3)

              !
              !- Loop over all interior cells
              !
              DO i = 1,nel

              rf     = 0.0d0  
              ruf    = 0.0d0   
              ruuf   = 0.0d0  
              sdf    = 0.0d0  
              ssf    = 0.0d0 
              rsssdf = 0.0d0
              rss2f  = 0.0d0  
              dtf    = 0.0d0 
              rssdtf = 0.0d0 
              rtf    = 0.0d0 
              rutf   = 0.0d0  

              !
              !- Start test filtering
              !- Data structures will vary based on code
              !- Loop over neighboring cells
              !
              do iel = 0,iee(0,i)     ! iee(0,i) = # of ngbrs

              if (iel==0) then
                  ii   = i            ! current cell
                  wgtf = 0.5d0 
              else
                  ii   = iee(iel,i)   ! ngbr cell
                  wgtf = 0.5d0/dble(iee(0,i))
              endif

              ! rho,u,v,w,t
              rq = primq(1,ii)
              uq = primq(2,ii)
              vq = primq(3,ii)
              wq = primq(4,ii)
              tq = primq(5,ii)

              ! dui/dxj
              duij(:,1) = grad(:,2,ii) ! gradients
              duij(:,2) = grad(:,3,ii)
              duij(:,3) = grad(:,4,ii)

              ! dt/dxi
              dti(:) = grad(:,5,ii)

              ! Sij
              sdij(1) = duij(1,1)
              sdij(2) = duij(2,2)
              sdij(3) = duij(3,3)

              sdij(4) = 0.5d0*( duij(1,2) + duij(2,1) )
              sdij(5) = 0.5d0*( duij(1,3) + duij(3,1) )
              sdij(6) = 0.5d0*( duij(2,3) + duij(3,2) )

              ! |S| = sqrt(2*Sij*Sij)
              ssij =  ( sdij(1)**2 + sdij(2)**2 + sdij(3)**2 )
     &        + 2.0d0*( sdij(4)**2 + sdij(5)**2 + sdij(6)**2 ) 

              ssij = sqrt(2.0d0*ssij)

              ! Sdij = Sij - delta_ij*Skk
              sskk = ( sdij(1) + sdij(1) + sdij(3) )*thd

              sdij(1) = sdij(1) - sskk 
              sdij(2) = sdij(2) - sskk 
              sdij(3) = sdij(3) - sskk 

              !
              !- Vars for Cs
              !
              ! rho
              rf = rf + rq*wgtf

              ! rho*ui
              ruf(1) = ruf(1) + rq*uq*wgtf
              ruf(2) = ruf(2) + rq*vq*wgtf
              ruf(3) = ruf(3) + rq*wq*wgtf

              ! rho*ui*uj
              ruuf(1) = ruuf(1) + rq*uq*uq*wgtf
              ruuf(2) = ruuf(2) + rq*vq*vq*wgtf
              ruuf(3) = ruuf(3) + rq*wq*wq*wgtf
              ruuf(4) = ruuf(4) + rq*uq*vq*wgtf
              ruuf(5) = ruuf(5) + rq*uq*wq*wgtf
              ruuf(6) = ruuf(6) + rq*vq*wq*wgtf

              ! Sdij = Sij - delta_ij*Skk
              sdf(:) = sdf(:) + sdij(:)*wgtf

              ! |S| = sqrt(2*Sij*Sij)
              ssf = ssf + ssij*wgtf

              ! rho*|s|*Sdij
              rsssdf(:) = rsssdf(:) + rq*ssij*sdij(:)*wgtf

              !
              !- Vars for Ci
              !
              ! rho*|S|*|S|
              rss2f = rss2f + rq*ssij*ssij*wgtf

              !
              !- Vars for Prt
              !
              ! dt/dxi 
              dtf(:) = dtf(:) + dti(:)*wgtf

              ! rho*|S|*dt/dxi
              rssdtf(:) = rssdtf(:) + rq*ssij*dti(:)*wgtf 

              ! rho*t
              rtf = rtf + rq*tq*wgtf

              ! rho*ui*t
              rutf(1) = rutf(1) + rq*uq*tq*wgtf
              rutf(2) = rutf(2) + rq*vq*tq*wgtf
              rutf(3) = rutf(3) + rq*wq*tq*wgtf

              enddo ! iel 

              !
              !- Vars for Cs
              !
              lij(1) = ruuf(1) - ruf(1)*ruf(1)/rf
              lij(2) = ruuf(2) - ruf(2)*ruf(2)/rf
              lij(3) = ruuf(3) - ruf(3)*ruf(3)/rf
              lij(4) = ruuf(4) - ruf(1)*ruf(2)/rf
              lij(5) = ruuf(5) - ruf(1)*ruf(3)/rf
              lij(6) = ruuf(6) - ruf(2)*ruf(3)/rf

              lkk = ( lij(1) + lij(2) + lij(3) )*thd

              lij(1) = lij(1) - lkk
              lij(2) = lij(2) - lkk
              lij(3) = lij(3) - lkk

              mij(:) = rsssdf(:) - ratf2*rf*ssf*sdf(:)

              lcs(i) =  lij(1)*mij(1) + lij(2)*mij(2) + lij(3)*mij(3)
     &        + 2.0d0*( lij(4)*mij(4) + lij(5)*mij(5) + lij(6)*mij(6) )

              mcs(i) =  mij(1)**2 + mij(2)**2 + mij(3)**2
     &        + 2.0d0*( mij(4)**2 + mij(5)**2 + mij(6)**2 )

              mcs(i) = 2.0d0*mcs(i)

              !
              !- Vars for Ci
              !
              lci(i) = 3.0d0*lkk

              mci(i) = 2.0d0*( rf*ratf2*ssf*ssf - rss2f )

              !
              !- Vars for Prt
              !
              lti(:) = rutf(:) - ruf(:)*rtf/rf

              mti(:) = rssdtf(:) - ratf2*rf*ssf*dtf(:)

              lpt(i) = mti(1)*lti(1) + mti(2)*lti(2) + mti(3)*lti(3)

              mpt(i) = mti(1)**2 + mti(2)**2 + mti(3)**2

              ENDDO ! nel

              !
              !- PLACEHOLDER FOR BC
              !- Data structures will vary based on code
              !- Set BCs for ghost cells: extrapolation
              !- for lcs(:),mcs(:),lci(:),mci(:),lpt(:),mpt(:)
              !
              do k = 1,nft-1
              do  j  = imaster(k,1),imaster(k,2)
              i  = ife(j,1) ! Face to element connectivity
              ii = ife(j,2)

              lcs(ii) = lcs(i)
              mcs(ii) = mcs(i)
              lci(ii) = lci(i)
              mci(ii) = mci(i)
              lpt(ii) = lpt(i)
              mpt(ii) = mpt(i)
              enddo
              enddo

              !
              !- PLACEHOLDER FOR MPI EXCHANGES FOR PERIODIC/SHARED BOUNDARIES
              !- Set shared cell values
              !- Set periodic ghost cell values
              !
              call dsm_exchange(1)

              return
          end subroutine dsm_fil

          !-----------------------------------------------------------
          ! Averaging routine for DSM                                
          !-----------------------------------------------------------
          ! Outputs:-                                               
          !-----------                                               
          ! cs_dsm(1:net), cs_dsm(1:net), pt_dsm(1:net)              
          !-----------------------------------------------------------
          subroutine dsm_avg()
              implicit none

              integer :: i,j,ii,iel,list(0:7),nli,k
              real(8) :: wgtf,dfw2

              ! Local averages
              real(8) :: lcs_l,mcs_l,lci_l,mci_l,lpt_l,mpt_l
              real(8) :: lcs_r,mcs_r,lci_r,mci_r,lpt_r,mpt_r
              real(8) :: lcs_loc,mcs_loc,lci_loc,mci_loc,lpt_loc,mpt_loc
              real(8) :: cs_loc,ci_loc,pt_loc

              !
              !- Local averaging 
              !- Loop over interior cells
              !
              DO i = 1,nel

              dfw2 = voli(i)**(-tthd) ! Filtering width
              list = ilst(0:7,i)
              nli  = list(0)

              lcs_loc = 0.0d0
              mcs_loc = 0.0d0
              lci_loc = 0.0d0
              mci_loc = 0.0d0
              lpt_loc = 0.0d0
              mpt_loc = 0.0d0

              lcs_l = lcs(i)
              mcs_l = mcs(i)
              lci_l = lci(i)
              mci_l = mci(i)
              lpt_l = lpt(i)
              mpt_l = mpt(i)

              !
              !- Loop over cell neighbors
              !- Data structures will vary based on code
              !
              do iel = 1,nli
              ii = list(iel)

              lcs_r = lcs(ii)
              mcs_r = mcs(ii)
              lci_r = lci(ii)
              mci_r = mci(ii)
              lpt_r = lpt(ii)
              mpt_r = mpt(ii)

              lcs_loc = lcs_loc + lcs_l + lcs_r
              mcs_loc = mcs_loc + mcs_l + mcs_r
              lci_loc = lci_loc + lci_l + lci_r
              mci_loc = mci_loc + mci_l + mci_r
              lpt_loc = lpt_loc + lpt_l + lpt_r
              mpt_loc = mpt_loc + mpt_l + mpt_r
              enddo 

              wgtf = 0.5d0/dble(nli)

              lcs_loc = lcs_loc*wgtf
              mcs_loc = mcs_loc*wgtf
              lci_loc = lci_loc*wgtf
              mci_loc = mci_loc*wgtf
              lpt_loc = lpt_loc*wgtf
              mpt_loc = mpt_loc*wgtf

              !- Cs clipping: 0.0 < Cs < 1.0
              if (lcs_loc<1d-9 .or. mcs_loc<1d-9) then
                  cs_dsm(i) = 0.0d0
              else
                  cs_dsm(i) = lcs_loc/mcs_loc
                  cs_loc    = sqrt(cs_dsm(i)/dfw2)
                  if (cs_loc>dsm_clip_cs) cs_dsm(i) = dfw2*dsm_clip_cs
              endif

              !- Ci clipping:  0.0 < Ci < 1.0
              if (lci_loc<1d-9 .or. mci_loc<1d-9) then
                  ci_dsm(i) = 0.0d0
              else
                  ci_dsm(i) = lci_loc/mci_loc
                  ci_loc    = ci_dsm(i)/dfw2
                  if (ci_loc>dsm_clip_ci) ci_dsm(i) = dfw2*dsm_clip_ci
              endif

              !- Prt clipping: 0.0 < Prt < 10.0
              if (lpt_loc<1d-9 .or. mpt_loc<1d-9) then
                  pt_dsm(i) = 0.0d0
              else
                  pt_dsm(i) = lpt_loc/mpt_loc
                  pt_loc    = cs_dsm(i)/pt_dsm(i)
                  if (pt_loc>dsm_clip_pt) pt_dsm(i) = cs_dsm(i)/dsm_clip_pt
              endif

              ENDDO ! nel

              !
              !- PLACEHOLDER FOR BC
              !- Data structures will vary based on code
              !- Set BCs for ghost cells
              !- for cs_dsm(:),ci_dsm(:),pt_dsm(:)
              !
              do k = 1,nft-1
              if (k>=11 .and. k<=16) then ! Skip periodic cells
                  cycle
              elseif (k==3) then ! BC for wall face 
                  do  j  = imaster(k,1),imaster(k,2)
                  i  = ife(j,1)
                  ii = ife(j,2)
                  cs_dsm(ii) = -cs_dsm(i)
                  ci_dsm(ii) = -ci_dsm(i)
                  pt_dsm(ii) = -pt_dsm(i)
                  enddo
              else ! Extrapolate for other faces
                  do  j  = imaster(k,1),imaster(k,2)
                  i  = ife(j,1)
                  ii = ife(j,2)
                  cs_dsm(ii) = cs_dsm(i)
                  ci_dsm(ii) = ci_dsm(i)
                  pt_dsm(ii) = pt_dsm(i)
                  enddo
              endif
              enddo

              !
              !- PLACEHOLDER FOR MPI EXCHANGES FOR PERIODIC/SHARED BOUNDARY
              !- Data structures will vary based on code
              !- Set shared cell values
              !- Set periodic cell values
              !- for cs_dsm(:),ci_dsm(:),pt_dsm(:)
              !
              call dsm_exchange(2)

              return
          end subroutine dsm_avg


          !
          !
          !- PLACEHOLDER FOR MPI EXCHANGES FOR PERIODIC/SHARED BOUNDARY
          !
          !
          subroutine dsm_exchange(iact)
              implicit none
              integer, intent(in) :: iact

          end subroutine dsm_exchange

          !==============================================================!
          !==============================================================!
