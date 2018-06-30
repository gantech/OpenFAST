!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of ExtLoads.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
!> ExtLoads is a time-domain loads module for horizontal-axis wind turbines.
module ExtLoads

   use NWTC_Library
   use ExtLoads_Types
   
   implicit none

   private

   ! ..... Public Subroutines ...................................................................................................

   public :: ExtLd_Init                           ! Initialization routine
   public :: ExtLd_End                            ! Ending routine (includes clean up)
   public :: ExtLd_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                     !   continuous states, and updating discrete states
   public :: ExtLd_CalcOutput                     ! Routine for computing outputs
  
contains    
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine sets the initialization output data structure, which contains data to be returned to the calling program (e.g.,
!! FAST)   
subroutine ExtLd_SetInitOut(p, InitOut, errStat, errMsg)

   type(ExtLd_InitOutputType),       intent(  out)  :: InitOut          ! output data
   type(ExtLd_ParameterType),        intent(in   )  :: p                ! Parameters
   integer(IntKi),                intent(  out)  :: errStat          ! Error status of the operation
   character(*),                  intent(  out)  :: errMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ExtLd_SetInitOut'
   
   
   
   integer(IntKi)                               :: i, j, k, f
   integer(IntKi)                               :: NumCoords
#ifdef DBG_OUTS
   integer(IntKi)                               :: m
   character(5)                                 ::chanPrefix
#endif   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
end subroutine ExtLd_SetInitOut

!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine ExtLd_Init( InitInp, u, p, y, m, interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   type(ExtLd_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   type(ExtLd_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(ExtLd_OutputType),          intent(  out) :: y             !< Initial system outputs (outputs are not calculated;
   type(ExtLd_MiscVarType),         intent(  out) :: m             !< Miscellaneous variables
   type(ExtLd_ParameterType),       intent(  out) :: p             !< Parameter variables
                                                                !!   only the output mesh is initialized)
   real(DbKi),                   intent(inout) :: interval      !< Coupling interval in seconds: the rate that
                                                                !!   (1) ExtLd_UpdateStates() is called in loose coupling &
                                                                !!   (2) ExtLd_UpdateDiscState() is called in tight coupling.
                                                                !!   Input is the suggested time from the glue code;
                                                                !!   Output is the actual coupling interval that will be used
                                                                !!   by the glue code.
   type(ExtLd_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),               intent(  out) :: errStat       !< Error status of the operation
   character(*),                 intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   

      ! Local variables
   integer(IntKi)                              :: i             ! loop counter
   
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
      
   character(*), parameter                     :: RoutineName = 'ExtLd_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

      ! Initialize the NWTC Subroutine Library

      ! Set parameters here
   p%NumBlds = InitInp%NumBlades
   call AllocAry(p%NumBldNds, p%NumBlds, 'NumBldNds', ErrStat2,ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
     if (ErrStat >= AbortErrLev) then
        call Cleanup()
        return
     end if
   p%NumBldNds(:) = InitInp%NumBldNodes(:)
   p%nTotBldNds = sum(p%NumBldNds(:))
   p%NumTwrNds = InitInp%NumTwrNds
   p%TwrAero = InitInp%TwrAero

      !............................................................................................
      ! Define and initialize inputs here 
      !............................................................................................
   
   call Init_u( u, p, InitInp, errStat2, errMsg2 ) 
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

      ! 
      !............................................................................................
      ! Define outputs here
      !............................................................................................
   call Init_y(y, u, p, errStat2, errMsg2) ! do this after input meshes have been initialized
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   
   
      !............................................................................................
      ! Define initialization output here
      !............................................................................................
   call ExtLd_SetInitOut(p, InitOut, errStat2, errMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   
   call Cleanup() 

 contains
   subroutine Cleanup()
     
   end subroutine Cleanup
   
end subroutine ExtLd_Init
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes ExtLoads meshes and output array variables for use during the simulation.
subroutine Init_y(y, u, p, errStat, errMsg)
   type(ExtLd_OutputType),           intent(  out)  :: y               !< Module outputs
   type(ExtLd_InputType),            intent(inout)  :: u               !< Module inputs -- intent(out) because of mesh sibling copy
   type(ExtLd_ParameterType),        intent(in   )  :: p               !< Parameters
   integer(IntKi),                intent(  out)  :: errStat         !< Error status of the operation
   character(*),                  intent(  out)  :: errMsg          !< Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: k                 ! loop counter for blades
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_y'

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

   if (p%TwrAero) then

      call MeshCopy ( SrcMesh  = u%TowerMotion    &
           , DestMesh = y%TowerLoad      &
           , CtrlCode = MESH_SIBLING     &
           , IOS      = COMPONENT_OUTPUT &
           , force    = .TRUE.           &
           , moment   = .TRUE.           &
           , ErrStat  = ErrStat2         &
           , ErrMess  = ErrMsg2          )

      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) RETURN         

      !y%TowerLoad%force = 0.0_ReKi  ! shouldn't have to initialize this
      !y%TowerLoad%moment= 0.0_ReKi  ! shouldn't have to initialize this
   else
      y%TowerLoad%nnodes = 0
   end if

   allocate( y%BladeLoad(p%NumBlds), stat=ErrStat2 )
   if (errStat2 /= 0) then
      call SetErrStat( ErrID_Fatal, 'Error allocating y%BladeLoad.', ErrStat, ErrMsg, RoutineName )      
      return
   end if

   do k = 1, p%NumBlds

      call MeshCopy ( SrcMesh  = u%BladeMotion(k) &
           , DestMesh = y%BladeLoad(k)   &
           , CtrlCode = MESH_SIBLING     &
           , IOS      = COMPONENT_OUTPUT &
           , force    = .TRUE.           &
           , moment   = .TRUE.           &
           , ErrStat  = ErrStat2         &
           , ErrMess  = ErrMsg2          )

      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 

   end do

   CALL AllocPAry( y%DX_y%twrLd, p%NumTwrNds*6, 'twrLd', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( y%DX_y%bldLd, p%nTotBldNds*6, 'bldLd', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! make sure the C versions are synced with these arrays
   y%DX_y%c_obj%twrLd_Len = p%NumTwrNds*6; y%DX_y%c_obj%twrLd = C_LOC( y%DX_y%twrLd(1) )
   y%DX_y%c_obj%bldLd_Len = p%nTotBldNds*6; y%DX_y%c_obj%bldLd = C_LOC( y%DX_y%bldLd(1) )

   call ConvertOpDataForExtProg(y, p, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   
end subroutine Init_y
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes ExtLoads meshes and input array variables for use during the simulation.
subroutine Init_u( u, p, InitInp, errStat, errMsg )
!..................................................................................................................................

  USE BeamDyn_IO, ONLY: BD_CrvExtractCrv
  
   type(ExtLd_InputType),           intent(  out)  :: u                 !< Input data
   type(ExtLd_ParameterType),       intent(in   )  :: p                 !< Parameters
   type(ExtLd_InitInputType),       intent(in   )  :: InitInp           !< Input data for ExtLd initialization routine
   integer(IntKi),               intent(  out)  :: errStat           !< Error status of the operation
   character(*),                 intent(  out)  :: errMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   real(reKi)                                   :: position(3)       ! node reference position
   real(reKi)                                   :: positionL(3)      ! node local position
   real(R8Ki)                                   :: theta(3)          ! Euler angles
   real(R8Ki)                                   :: orientation(3,3)  ! node reference orientation
   real(R8Ki)                                   :: orientationL(3,3) ! node local orientation
   
   real(R8Ki)                                   :: wm_crv(3)         ! Wiener-Milenkovic parameters
   integer(IntKi)                               :: j                 ! counter for nodes
   integer(IntKi)                               :: jTot              ! counter for blade nodes
   integer(IntKi)                               :: k                 ! counter for blades

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_u'

      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Meshes for motion inputs (ElastoDyn and/or BeamDyn)
         !................
         ! tower
         !................
   if (p%NumTwrNds > 0) then
      
      call MeshCreate ( BlankMesh = u%TowerMotion   &
                       ,IOS       = COMPONENT_INPUT &
                       ,Nnodes    = p%NumTwrNds     &
                       ,ErrStat   = ErrStat2        &
                       ,ErrMess   = ErrMsg2         &
                       ,Orientation     = .true.    &
                       ,TranslationDisp = .true.    &
                       ,TranslationVel  = .true.    &
                       ,RotationVel = .true.        &
                      )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      if (errStat >= AbortErrLev) return
            
         ! set node initial position/orientation
      position = 0.0_ReKi
      do j=1,p%NumTwrNds         
         position(:) = InitInp%TwrPos(:,j)
         
         call MeshPositionNode(u%TowerMotion, j, position, errStat2, errMsg2)  ! orientation is identity by default
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
         
         ! create point elements
      do j=1,p%NumTwrNds
         call MeshConstructElement( u%TowerMotion, ELEMENT_POINT, errStat2, errMsg2, p1=j )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
            
      call MeshCommit(u%TowerMotion, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      if (errStat >= AbortErrLev) return

      
      u%TowerMotion%Orientation     = u%TowerMotion%RefOrientation
      u%TowerMotion%TranslationDisp = 0.0_R8Ki
      u%TowerMotion%TranslationVel  = 0.0_ReKi
      u%TowerMotion%RotationVel = 0.0_ReKi
      
   end if ! we compute tower loads
   
         !................
         ! hub
         !................
   
      call MeshCreate ( BlankMesh = u%HubMotion     &
                       ,IOS       = COMPONENT_INPUT &
                       ,Nnodes    = 1               &
                       ,ErrStat   = ErrStat2        &
                       ,ErrMess   = ErrMsg2         &
                       ,Orientation     = .true.    &
                       ,TranslationDisp = .true.    &
                       ,TranslationVel  = .true.    &
                       ,RotationVel     = .true.    &
                      )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      if (errStat >= AbortErrLev) return
                     
      call MeshPositionNode(u%HubMotion, 1, InitInp%HubPos, errStat2, errMsg2, InitInp%HubOrient)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
      call MeshConstructElement( u%HubMotion, ELEMENT_POINT, errStat2, errMsg2, p1=1 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      call MeshCommit(u%HubMotion, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      if (errStat >= AbortErrLev) return

         
      u%HubMotion%Orientation     = u%HubMotion%RefOrientation
      u%HubMotion%TranslationDisp = 0.0_R8Ki
      u%HubMotion%TranslationVel = 0.0_R8Ki
      u%HubMotion%RotationVel     = 0.0_R8Ki   

         !................
         ! nacelle
         !................
   
      call MeshCreate ( BlankMesh = u%NacelleMotion     &
                       ,IOS       = COMPONENT_INPUT &
                       ,Nnodes    = 1               &
                       ,ErrStat   = ErrStat2        &
                       ,ErrMess   = ErrMsg2         &
                       ,Orientation     = .true.    &
                       ,TranslationDisp = .true.    &
                       ,TranslationVel  = .true.    &
                       ,RotationVel     = .true.    &
                      )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      if (errStat >= AbortErrLev) return
                     
      call MeshPositionNode(u%NacelleMotion, 1, InitInp%NacellePos, errStat2, errMsg2, InitInp%NacelleOrient)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
      call MeshConstructElement( u%NacelleMotion, ELEMENT_POINT, errStat2, errMsg2, p1=1 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      call MeshCommit(u%NacelleMotion, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      if (errStat >= AbortErrLev) return

         
      u%NacelleMotion%Orientation     = u%NacelleMotion%RefOrientation
      u%NacelleMotion%TranslationDisp = 0.0_R8Ki
      u%NacelleMotion%TranslationVel = 0.0_R8Ki
      u%NacelleMotion%RotationVel     = 0.0_R8Ki   
      
         !................
         ! blades
         !................

      allocate( u%BladeRootMotion(p%NumBlds), STAT = ErrStat2 )
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, 'Error allocating u%BladeRootMotion array.', ErrStat, ErrMsg, RoutineName )
         return
      end if
      
      allocate( u%BladeMotion(p%NumBlds), STAT = ErrStat2 )
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, 'Error allocating u%BladeMotion array.', ErrStat, ErrMsg, RoutineName )
         return
      end if
      
      do k=1,p%NumBlds


         call MeshCreate ( BlankMesh = u%BladeRootMotion(k)     &
              ,IOS       = COMPONENT_INPUT &
              ,Nnodes    = 1               &
              ,ErrStat   = ErrStat2        &
              ,ErrMess   = ErrMsg2         &
              ,Orientation     = .true.    &
              ,TranslationDisp = .true.    &
              ,TranslationVel  = .true.    &
              ,RotationVel     = .true.    &
              )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
         if (errStat >= AbortErrLev) return
         
         call MeshPositionNode(u%BladeRootMotion(k), 1, InitInp%BldRootPos(:,k), errStat2, errMsg2, InitInp%BldRootOrient(:,:,k))
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
         call MeshConstructElement( u%BladeRootMotion(k), ELEMENT_POINT, errStat2, errMsg2, p1=1 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
         call MeshCommit(u%BladeRootMotion(k), errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
         if (errStat >= AbortErrLev) return
         
         u%BladeRootMotion(k)%Orientation     = u%BladeRootMotion(k)%RefOrientation
         u%BladeRootMotion(k)%TranslationDisp = 0.0_R8Ki
         u%BladeRootMotion(k)%TranslationVel  = 0.0_R8Ki
         u%BladeRootMotion(k)%RotationVel     = 0.0_R8Ki   
         
         call MeshCreate ( BlankMesh = u%BladeMotion(k)                     &
                          ,IOS       = COMPONENT_INPUT                      &
                          ,Nnodes    = InitInp%NumBldNodes(k) &
                          ,ErrStat   = ErrStat2                             &
                          ,ErrMess   = ErrMsg2                              &
                          ,Orientation     = .true.                         &
                          ,TranslationDisp = .true.                         &
                          ,TranslationVel  = .true.                         &
                          ,RotationVel = .true.                             &
                         )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

         if (errStat >= AbortErrLev) return
            
                        
         do j=1,InitInp%NumBldNodes(k)

               ! reference position of the jth node in the kth blade:
            position(:) = InitInp%BldPos(:,j,k)
                                 
               ! reference orientation of the jth node in the kth blade
            orientation(:,:) = InitInp%BldOrient(:,:,j,k)

            
            call MeshPositionNode(u%BladeMotion(k), j, position, errStat2, errMsg2, orientation)
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
         end do ! j=blade nodes
         
            ! create point elements
         do j=1,InitInp%NumBldNodes(k)
            call MeshConstructElement( u%BladeMotion(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         end do !j
            
         call MeshCommit(u%BladeMotion(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
         if (errStat >= AbortErrLev) return
      
         u%BladeMotion(k)%Orientation     = u%BladeMotion(k)%RefOrientation
         u%BladeMotion(k)%TranslationDisp = 0.0_R8Ki
         u%BladeMotion(k)%TranslationVel  = 0.0_R8Ki
         u%BladeMotion(k)%RotationVel = 0.0_R8Ki
   
   end do !k=numBlades

   ! Set the parameters first
   u%DX_u%nTowerNodes = p%NumTwrNds
   u%DX_u%nBlades = p%NumBlds
   CALL AllocPAry( u%DX_u%nBladeNodes, p%NumBlds, 'nBladeNodes', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   u%DX_u%c_obj%nBladeNodes_Len = p%NumBlds; u%DX_u%c_obj%nBladeNodes = C_LOC( u%DX_u%nBladeNodes(1) )
   u%DX_u%nBladeNodes(:) = p%NumBldNds(:)

   ! Set the reference positions next
   CALL AllocPAry( u%DX_u%twrRefPos, p%NumTwrNds*6, 'twrRefPos', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( u%DX_u%bldRefPos, p%nTotBldNds*6, 'bldRefPos', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( u%DX_u%hubRefPos, 6, 'hubRefPos', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( u%DX_u%nacRefPos, 6, 'nacRefPos', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! make sure the C versions are synced with these arrays
   u%DX_u%c_obj%twrRefPos_Len = p%NumTwrNds*6; u%DX_u%c_obj%twrRefPos = C_LOC( u%DX_u%twrRefPos(1) )
   u%DX_u%c_obj%bldRefPos_Len = p%nTotBldNds*6; u%DX_u%c_obj%bldRefPos = C_LOC( u%DX_u%bldRefPos(1) )
   u%DX_u%c_obj%hubRefPos_Len = 6; u%DX_u%c_obj%hubRefPos = C_LOC( u%DX_u%hubRefPos(1) )
   u%DX_u%c_obj%nacRefPos_Len = 6; u%DX_u%c_obj%nacRefPos = C_LOC( u%DX_u%nacRefPos(1) )
   
   if (p%TwrAero) then
      do j=1,p%NumTwrNds
         call BD_CrvExtractCrv(u%TowerMotion%RefOrientation(:,:,j), wm_crv, ErrStat2, ErrMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         u%DX_u%twrRefPos((j-1)*6+1:(j-1)*6+3) = u%TowerMotion%Position(:,j)
         u%DX_u%twrRefPos((j-1)*6+4:(j-1)*6+6) = wm_crv
      end do
   end if

   jTot = 1
   do k=1,p%NumBlds
      do j=1,p%NumBldNds(k)
         call BD_CrvExtractCrv(u%BladeMotion(k)%RefOrientation(:,:,j), wm_crv, ErrStat2, ErrMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         u%DX_u%bldRefPos((jTot-1)*6+1:(jTot-1)*6+3) = u%BladeMotion(k)%Position(:,j)
         u%DX_u%bldRefPos((jTot-1)*6+4:(jTot-1)*6+6) = wm_crv
         jTot = jTot+1
      end do
   end do

   call BD_CrvExtractCrv(u%HubMotion%RefOrientation(:,:,1), wm_crv, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   u%DX_u%hubRefPos(1:3) = u%HubMotion%Position(:,1)
   u%DX_u%hubRefPos(4:6) = wm_crv

   call BD_CrvExtractCrv(u%NacelleMotion%RefOrientation(:,:,1), wm_crv, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   u%DX_u%nacRefPos(1:3) = u%NacelleMotion%Position(:,1)
   u%DX_u%nacRefPos(4:6) = wm_crv

   ! Now the displacements
   CALL AllocPAry( u%DX_u%twrDef, p%NumTwrNds*12, 'twrDef', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( u%DX_u%bldDef, p%nTotBldNds*12, 'bldDef', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( u%DX_u%hubDef, 12, 'hubDef', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( u%DX_u%nacDef, 12, 'nacDef', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   
   ! make sure the C versions are synced with these arrays
   u%DX_u%c_obj%twrDef_Len = p%NumTwrNds*12; u%DX_u%c_obj%twrDef = C_LOC( u%DX_u%twrDef(1) )
   u%DX_u%c_obj%bldDef_Len = p%nTotBldNds*12; u%DX_u%c_obj%bldDef = C_LOC( u%DX_u%bldDef(1) )
   u%DX_u%c_obj%hubDef_Len = 12; u%DX_u%c_obj%hubDef = C_LOC( u%DX_u%hubDef(1) )
   u%DX_u%c_obj%nacDef_Len = 12; u%DX_u%c_obj%nacDef = C_LOC( u%DX_u%nacDef(1) )      
   call ConvertInpDataForExtProg(u, p, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   
end subroutine Init_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine converts the displacement data in the meshes in the input into a simple array format that can be accessed by external programs
subroutine ConvertInpDataForExtProg(u, p, errStat, errMsg )
!..................................................................................................................................
  USE BeamDyn_IO, ONLY: BD_CrvExtractCrv
  
   type(ExtLd_InputType),           intent(inout)  :: u                 !< Input data
   type(ExtLd_ParameterType),       intent(in   )  :: p                 !< Parameters
   integer(IntKi),               intent(  out)  :: errStat           !< Error status of the operation
   character(*),                 intent(  out)  :: errMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   real(R8Ki)                                   :: wm_crv(3)         ! Wiener-Milenkovic parameters
   integer(intKi)                               :: j                 ! counter for nodes
   integer(intKi)                               :: jTot              ! counter for nodes
   integer(intKi)                               :: k                 ! counter for blades

   
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ConvertInpDataForExtProg'

      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (p%TwrAero) then
      do j=1,p%NumTwrNds
         call BD_CrvExtractCrv(u%TowerMotion%Orientation(:,:,j), wm_crv, ErrStat2, ErrMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         u%DX_u%twrDef((j-1)*12+1:(j-1)*12+3) = u%TowerMotion%TranslationDisp(:,j)
         u%DX_u%twrDef((j-1)*12+4:(j-1)*12+6) = u%TowerMotion%TranslationVel(:,j)
         u%DX_u%twrDef((j-1)*12+7:(j-1)*12+9) = wm_crv
         u%DX_u%twrDef((j-1)*12+10:(j-1)*12+12) = u%TowerMotion%RotationVel(:,j)
      end do
   end if

   jTot = 1
   do k=1,p%NumBlds
      do j=1,p%NumBldNds(k)
         call BD_CrvExtractCrv(u%BladeMotion(k)%Orientation(:,:,j), wm_crv, ErrStat2, ErrMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         u%DX_u%bldDef((jTot-1)*12+1:(jTot-1)*12+3) = u%BladeMotion(k)%TranslationDisp(:,j)
         u%DX_u%bldDef((jTot-1)*12+4:(jTot-1)*12+6) = u%BladeMotion(k)%TranslationVel(:,j)
         u%DX_u%bldDef((jTot-1)*12+7:(jTot-1)*12+9) = wm_crv
         u%DX_u%bldDef((jTot-1)*12+10:(jTot-1)*12+12) = u%BladeMotion(k)%RotationVel(:,j)
         jTot = jTot+1
      end do
   end do


   call BD_CrvExtractCrv(u%HubMotion%Orientation(:,:,1), wm_crv, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   u%DX_u%hubDef(1:3) = u%HubMotion%TranslationDisp(:,1)
   u%DX_u%hubDef(4:6) = u%HubMotion%TranslationVel(:,1)
   u%DX_u%hubDef(7:9) = wm_crv
   u%DX_u%hubDef(10:12) = u%HubMotion%RotationVel(:,1)

   call BD_CrvExtractCrv(u%NacelleMotion%Orientation(:,:,1), wm_crv, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   u%DX_u%nacDef(1:3) = u%NacelleMotion%TranslationDisp(:,1)
   u%DX_u%nacDef(4:6) = u%NacelleMotion%TranslationVel(:,1)
   u%DX_u%nacDef(7:9) = wm_crv
   u%DX_u%nacDef(10:12) = u%NacelleMotion%RotationVel(:,1)
   
   
   
end subroutine ConvertInpDataForExtProg
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine converts the data in the simple array format in the output data type into OpenFAST mesh format
subroutine ConvertOpDataForExtProg(y, p, errStat, errMsg )
!..................................................................................................................................
  
   type(ExtLd_OutputType),          intent(inout)  :: y                 !< Ouput data
   type(ExtLd_ParameterType),       intent(in   )  :: p                 !< Parameters
   integer(IntKi),               intent(  out)  :: errStat           !< Error status of the operation
   character(*),                 intent(  out)  :: errMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: j                 ! counter for nodes
   integer(intKi)                               :: jTot              ! counter for nodes
   integer(intKi)                               :: k                 ! counter for blades

   
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ConvertInpDataForExtProg'

      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (p%TwrAero) then
      do j=1,p%NumTwrNds
         y%TowerLoad%Force(:,j) = y%DX_y%twrLd((j-1)*6+1:(j-1)*6+3)
         y%TowerLoad%Moment(:,j) = y%DX_y%twrLd((j-1)*6+4:(j-1)*6+6)
      end do
   end if

   jTot = 1
   do k=1,p%NumBlds
      do j=1,p%NumBldNds(k)
         y%BladeLoad(k)%Force(:,j) = y%DX_y%bldLd((jTot-1)*6+1:(jTot-1)*6+3)
         y%BladeLoad(k)%Moment(:,j) = y%DX_y%bldLd((jTot-1)*6+4:(jTot-1)*6+6)
         jTot = jTot+1
      end do
   end do
   
   
end subroutine ConvertOpDataForExtProg
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine ExtLd_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(ExtLd_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(ExtLd_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
      TYPE(ExtLd_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(ExtLd_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(ExtLd_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(ExtLd_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(ExtLd_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(ExtLd_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:


         ! Destroy the input data:


END SUBROUTINE ExtLd_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine ExtLd_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!..................................................................................................................................

   real(DbKi),                     intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(ExtLd_InputType),             intent(inout) :: u(:)       !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                     intent(in   ) :: utimes(:)  !< Times associated with u(:), in seconds
   type(ExtLd_ParameterType),         intent(in   ) :: p          !< Parameters
   type(ExtLd_ContinuousStateType),   intent(inout) :: x          !< Input: Continuous states at t;
                                                               !!   Output: Continuous states at t + Interval
   type(ExtLd_DiscreteStateType),     intent(inout) :: xd         !< Input: Discrete states at t;
                                                               !!   Output: Discrete states at t  + Interval
   type(ExtLd_ConstraintStateType),   intent(inout) :: z          !< Input: Constraint states at t;
                                                               !!   Output: Constraint states at t+dt
   type(ExtLd_OtherStateType),        intent(inout) :: OtherState !< Input: Other states at t;
                                                               !!   Output: Other states at t+dt
   type(ExtLd_MiscVarType),           intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                 intent(  out) :: errStat    !< Error status of the operation
   character(*),                   intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(ExtLd_InputType)                           :: uInterp     ! Interpolated/Extrapolated input
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ExtLd_UpdateStates'
      
   ErrStat = ErrID_None
   ErrMsg  = ""
           
   
end subroutine ExtLd_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine ExtLd_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(ExtLd_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
   TYPE(ExtLd_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ExtLd_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(ExtLd_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(ExtLd_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(ExtLd_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(ExtLd_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(ExtLd_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   integer, parameter                           :: indx = 1  ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt
   integer(intKi)                               :: i
   integer(intKi)                               :: j

   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'ExtLd_CalcOutput'
   real(ReKi)                                   :: SigmaCavitCrit, SigmaCavit

   ErrStat = ErrID_None
   ErrMsg  = ""

 end subroutine ExtLd_CalcOutput
 
END MODULE ExtLoads
