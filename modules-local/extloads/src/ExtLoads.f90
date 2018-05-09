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

   public :: ExtLoads_Init                           ! Initialization routine
   public :: ExtLoads_End                            ! Ending routine (includes clean up)
   public :: ExtLoads_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                     !   continuous states, and updating discrete states
   public :: ExtLoads_CalcOutput                     ! Routine for computing outputs
  
contains    
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine sets the initialization output data structure, which contains data to be returned to the calling program (e.g.,
!! FAST)   
subroutine ExtLoads_SetInitOut(p, InitOut, errStat, errMsg)

   type(ExtLoads_InitOutputType),       intent(  out)  :: InitOut          ! output data
   type(ExtLoads_ParameterType),        intent(in   )  :: p                ! Parameters
   integer(IntKi),                intent(  out)  :: errStat          ! Error status of the operation
   character(*),                  intent(  out)  :: errMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ExtLoads_SetInitOut'
   
   
   
   integer(IntKi)                               :: i, j, k, f
   integer(IntKi)                               :: NumCoords
#ifdef DBG_OUTS
   integer(IntKi)                               :: m
   character(5)                                 ::chanPrefix
#endif   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
end subroutine ExtLoads_SetInitOut

!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine ExtLoads_Init( InitInp, u, p, y, m, interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   type(ExtLoads_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   type(ExtLoads_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(ExtLoads_OutputType),          intent(  out) :: y             !< Initial system outputs (outputs are not calculated;
   type(ExtLoads_MiscVarType),         intent(  out) :: m             !< Miscellaneous variables
   type(ExtLoads_ParameterType),       intent(  out) :: p             !< Parameter variables
                                                                !!   only the output mesh is initialized)
   real(DbKi),                   intent(inout) :: interval      !< Coupling interval in seconds: the rate that
                                                                !!   (1) ExtLoads_UpdateStates() is called in loose coupling &
                                                                !!   (2) ExtLoads_UpdateDiscState() is called in tight coupling.
                                                                !!   Input is the suggested time from the glue code;
                                                                !!   Output is the actual coupling interval that will be used
                                                                !!   by the glue code.
   type(ExtLoads_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),               intent(  out) :: errStat       !< Error status of the operation
   character(*),                 intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   

      ! Local variables
   integer(IntKi)                              :: i             ! loop counter
   
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
      
   character(*), parameter                     :: RoutineName = 'ExtLoads_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )


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
   call ExtLoads_SetInitOut(p, InitOut, errStat2, errMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   
   call Cleanup() 

 contains
   subroutine Cleanup()
     
   end subroutine Cleanup
   
end subroutine ExtLoads_Init
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes ExtLoads meshes and output array variables for use during the simulation.
subroutine Init_y(y, u, p, errStat, errMsg)
   type(ExtLoads_OutputType),           intent(  out)  :: y               !< Module outputs
   type(ExtLoads_InputType),            intent(inout)  :: u               !< Module inputs -- intent(out) because of mesh sibling copy
   type(ExtLoads_ParameterType),        intent(in   )  :: p               !< Parameters
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

   allocate( y%BladeLoad(p%numBlades), stat=ErrStat2 )
   if (errStat2 /= 0) then
      call SetErrStat( ErrID_Fatal, 'Error allocating y%BladeLoad.', ErrStat, ErrMsg, RoutineName )      
      return
   end if

   do k = 1, p%numBlades

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
   
end subroutine Init_y
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes ExtLoads meshes and input array variables for use during the simulation.
subroutine Init_u( u, p, InitInp, errStat, errMsg )
!..................................................................................................................................

   type(ExtLoads_InputType),           intent(  out)  :: u                 !< Input data
   type(ExtLoads_ParameterType),       intent(in   )  :: p                 !< Parameters
   type(ExtLoads_InitInputType),       intent(in   )  :: InitInp           !< Input data for ExtLoads initialization routine
   integer(IntKi),               intent(  out)  :: errStat           !< Error status of the operation
   character(*),                 intent(  out)  :: errMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   real(reKi)                                   :: position(3)       ! node reference position
   real(reKi)                                   :: positionL(3)      ! node local position
   real(R8Ki)                                   :: theta(3)          ! Euler angles
   real(R8Ki)                                   :: orientation(3,3)  ! node reference orientation
   real(R8Ki)                                   :: orientationL(3,3) ! node local orientation
   
   integer(intKi)                               :: j                 ! counter for nodes
   integer(intKi)                               :: k                 ! counter for blades
   
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
      u%HubMotion%RotationVel     = 0.0_ReKi   
      
         !................
         ! blades
         !................
   
      allocate( u%BladeMotion(p%NumBlades), STAT = ErrStat2 )
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, 'Error allocating u%BladeMotion array.', ErrStat, ErrMsg, RoutineName )
         return
      end if
      
      do k=1,p%NumBlades
         call MeshCreate ( BlankMesh = u%BladeMotion(k)                     &
                          ,IOS       = COMPONENT_INPUT                      &
                          ,Nnodes    = InitInp%NumBldNodes(k) &
                          ,ErrStat   = ErrStat2                             &
                          ,ErrMess   = ErrMsg2                              &
                          ,Orientation     = .true.                         &
                          ,TranslationDisp = .true.                         &
                          ,TranslationVel  = .true.                         &
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
         u%BladeMotion(k)%TranslationVel  = 0.0_ReKi
   
   end do !k=numBlades
   
   
end subroutine Init_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine ExtLoads_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(ExtLoads_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(ExtLoads_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
      TYPE(ExtLoads_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(ExtLoads_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(ExtLoads_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(ExtLoads_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(ExtLoads_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(ExtLoads_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:


         ! Destroy the input data:


END SUBROUTINE ExtLoads_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine ExtLoads_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!..................................................................................................................................

   real(DbKi),                     intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(ExtLoads_InputType),             intent(inout) :: u(:)       !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                     intent(in   ) :: utimes(:)  !< Times associated with u(:), in seconds
   type(ExtLoads_ParameterType),         intent(in   ) :: p          !< Parameters
   type(ExtLoads_ContinuousStateType),   intent(inout) :: x          !< Input: Continuous states at t;
                                                               !!   Output: Continuous states at t + Interval
   type(ExtLoads_DiscreteStateType),     intent(inout) :: xd         !< Input: Discrete states at t;
                                                               !!   Output: Discrete states at t  + Interval
   type(ExtLoads_ConstraintStateType),   intent(inout) :: z          !< Input: Constraint states at t;
                                                               !!   Output: Constraint states at t+dt
   type(ExtLoads_OtherStateType),        intent(inout) :: OtherState !< Input: Other states at t;
                                                               !!   Output: Other states at t+dt
   type(ExtLoads_MiscVarType),           intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                 intent(  out) :: errStat    !< Error status of the operation
   character(*),                   intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(ExtLoads_InputType)                           :: uInterp     ! Interpolated/Extrapolated input
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ExtLoads_UpdateStates'
      
   ErrStat = ErrID_None
   ErrMsg  = ""
           
   
end subroutine ExtLoads_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine ExtLoads_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(ExtLoads_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
   TYPE(ExtLoads_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(ExtLoads_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(ExtLoads_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(ExtLoads_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(ExtLoads_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(ExtLoads_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(ExtLoads_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   integer, parameter                           :: indx = 1  ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt
   integer(intKi)                               :: i
   integer(intKi)                               :: j

   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'ExtLoads_CalcOutput'
   real(ReKi)                                   :: SigmaCavitCrit, SigmaCavit

   ErrStat = ErrID_None
   ErrMsg  = ""

 end subroutine ExtLoads_CalcOutput
 
END MODULE ExtLoads
