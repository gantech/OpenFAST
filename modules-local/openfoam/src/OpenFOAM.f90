!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    Lidar module, a submodule of InflowWind
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
! File last committed: $Date: $
! (File) Revision #: $Rev: $
! URL: $HeadURL: $
!**********************************************************************************************************************************
MODULE OpenFOAM

! This is a pseudo module used to couple FAST v8 with OpenFOAM; it is considered part of the FAST glue code
   USE FAST_Types
!   USE OpenFOAM_IO

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: OpFM_Ver = ProgDesc( 'OpenFOAM Integration', '', '' )


! ==================================================================================================="


      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Init_OpFM                           ! Initialization routine
   PUBLIC :: OpFM_SetInputs                      ! Glue-code routine to update inputs for OpenFOAM
   PUBLIC :: OpFM_SetWriteOutput


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_OpFM( InitInp, p_FAST, AirDens, u_AD14, u_AD, initOut_AD, y_AD, y_ED, y_BD, OpFM, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(OpFM_InitInputType),        INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   REAL(ReKi),                      INTENT(IN   )  :: AirDens     ! Air Density kg/m^3
   TYPE(AD14_InputType),            INTENT(IN   )  :: u_AD14      ! AeroDyn14 input data
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! AeroDyn input data
   TYPE(AD_OutputType),             INTENT(IN   )  :: y_AD        ! AeroDyn output data (for mesh mapping)
   TYPE(AD_InitOutputType),         INTENT(IN   )  :: initOut_AD  ! AeroDyn InitOutput data (for BladeProps)
   TYPE(ED_OutputType),             INTENT(IN)     :: y_ED        ! The outputs of the structural dynamics module
   TYPE(BD_OutputType),             INTENT(IN)     :: y_BD(:)     ! The outputs of the structural dynamics module   
   TYPE(OpenFOAM_Data),             INTENT(INOUT)  :: OpFM        ! data for the OpenFOAM integration module
   TYPE(OpFM_InitOutputType),       INTENT(INOUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                   :: NMappings  ! number of blades
   INTEGER(IntKi)                                   :: k          ! blade loop counter
   INTEGER(IntKi)                                   :: j          ! node counter

   INTEGER(IntKi)                                   :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                        :: RoutineName = 'Init_OpFM'

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

      !............................................................................................
      ! Define parameters here:
      !............................................................................................

      ! number of velocity nodes in the interface:

   OpFM%p%NnodesVel = 1  ! always want the hub point
   IF ( p_FAST%CompAero  == Module_AD14 ) THEN ! AeroDyn 14 needs these velocities
      CALL SetErrStat(ErrID_Fatal, 'Error AeroDyn14 is not supported yet with different number of velocity and force actuator nodes', ErrStat, ErrMsg, RoutineName)
      RETURN
   ELSEIF ( p_FAST%CompAero  == Module_AD ) THEN ! AeroDyn 15 needs these velocities
      OpFM%p%NumBl = InitInp%NumBl

      OpFM%p%NnodesVel = OpFM%p%NnodesVel + y_AD%TowerLoad%NNodes                   ! tower nodes (if any)
      DO k=1,OpFM%p%NumBl
         OpFM%p%NnodesVel = OpFM%p%NnodesVel + u_AD%BladeMotion(k)%NNodes           ! blade nodes
      END DO
   END IF

      ! number of force nodes in the interface
   Opfm%p%NnodesForceBlade =  InitInp%NumActForcePtsBlade 
   OpFM%p%NnodesForceTower = InitInp%NumActForcePtsTower
   OpFM%p%NnodesForce = 1 + OpFM%p%NumBl * InitInp%NumActForcePtsBlade

   if ( y_AD%TowerLoad%NNodes > 0 ) then
      OpFM%p%NMappings = OpFM%p%NumBl + 1
      OpFM%p%NnodesForce = OpFM%p%NnodesForce + InitInp%NumActForcePtsTower
   else
      OpFM%p%NMappings = OpFM%p%NumBl
      OpFM%p%NnodesForceTower = 0
   end if

      ! air density, required for normalizing values sent to OpenFOAM:
   OpFM%p%AirDens = AirDens
   if ( EqualRealNos( AirDens, 0.0_ReKi ) ) &
      CALL SetErrStat( ErrID_Fatal, 'Air density cannot be zero for OpenFOAM integration. Check that AeroDyn is used and that air density is set properly', ErrStat,ErrMsg,RoutineName)

      !............................................................................................
      ! Allocate arrays and define initial guesses for the OpenFOAM inputs here:
      !............................................................................................
   CALL AllocPAry( OpFM%u%pxVel, OpFM%p%NnodesVel, 'pxVel', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pyVel, OpFM%p%NnodesVel, 'pyVel', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pzVel, OpFM%p%NnodesVel, 'pzVel', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pxDotVel, OpFM%p%NnodesVel, 'pxDotVel', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pyDotVel, OpFM%p%NnodesVel, 'pyDotVel', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pzDotVel, OpFM%p%NnodesVel, 'pzDotVel', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pxForce, OpFM%p%NnodesForce, 'pxForce', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pyForce, OpFM%p%NnodesForce, 'pyForce', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pzForce, OpFM%p%NnodesForce, 'pzForce', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pxDotForce, OpFM%p%NnodesForce, 'pxDotForce', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pyDotForce, OpFM%p%NnodesForce, 'pyDotForce', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pzDotForce, OpFM%p%NnodesForce, 'pzDotForce', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL AllocPAry( OpFM%u%pOrientation, 3*3*OpFM%p%NnodesForce, 'pOrientation', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%fx, OpFM%p%NnodesForce, 'fx', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%fy, OpFM%p%NnodesForce, 'fy', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%fz, OpFM%p%NnodesForce, 'fz', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%momentx, OpFM%p%NnodesForce, 'momentx', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%momenty, OpFM%p%NnodesForce, 'momenty', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%momentz, OpFM%p%NnodesForce, 'momentz', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%forceNodesChord, OpFM%p%NnodesForce, 'forceNodesChord', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   
   IF (InitInp%NumCtrl2SC > 0) THEN
      CALL AllocPAry( OpFM%u%SuperController, InitInp%NumCtrl2SC, 'u%SuperController', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF

   IF (ErrStat >= AbortErrLev) RETURN

      ! make sure the C versions are synced with these arrays
   OpFM%u%c_obj%pxVel_Len = OpFM%p%NnodesVel; OpFM%u%c_obj%pxVel = C_LOC( OpFM%u%pxVel(1) )
   OpFM%u%c_obj%pyVel_Len = OpFM%p%NnodesVel; OpFM%u%c_obj%pyVel = C_LOC( OpFM%u%pyVel(1) )
   OpFM%u%c_obj%pzVel_Len = OpFM%p%NnodesVel; OpFM%u%c_obj%pzVel = C_LOC( OpFM%u%pzVel(1) )
   OpFM%u%c_obj%pxDotVel_Len = OpFM%p%NnodesVel; OpFM%u%c_obj%pxDotVel = C_LOC( OpFM%u%pxDotVel(1) )
   OpFM%u%c_obj%pyDotVel_Len = OpFM%p%NnodesVel; OpFM%u%c_obj%pyDotVel = C_LOC( OpFM%u%pyDotVel(1) )
   OpFM%u%c_obj%pzDotVel_Len = OpFM%p%NnodesVel; OpFM%u%c_obj%pzDotVel = C_LOC( OpFM%u%pzDotVel(1) )
   OpFM%u%c_obj%pxForce_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%pxForce = C_LOC( OpFM%u%pxForce(1) )
   OpFM%u%c_obj%pyForce_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%pyForce = C_LOC( OpFM%u%pyForce(1) )
   OpFM%u%c_obj%pzForce_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%pzForce = C_LOC( OpFM%u%pzForce(1) )
   OpFM%u%c_obj%pxDotForce_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%pxDotForce = C_LOC( OpFM%u%pxDotForce(1) )
   OpFM%u%c_obj%pyDotForce_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%pyDotForce = C_LOC( OpFM%u%pyDotForce(1) )
   OpFM%u%c_obj%pzDotForce_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%pzDotForce = C_LOC( OpFM%u%pzDotForce(1) )
   OpFM%u%c_obj%pOrientation_Len = OpFM%p%NnodesForce*3*3; OpFM%u%c_obj%pOrientation = C_LOC( OpFM%u%pOrientation(1) )
   OpFM%u%c_obj%fx_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%fx = C_LOC( OpFM%u%fx(1) )
   OpFM%u%c_obj%fy_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%fy = C_LOC( OpFM%u%fy(1) )
   OpFM%u%c_obj%fz_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%fz = C_LOC( OpFM%u%fz(1) )
   OpFM%u%c_obj%momentx_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%momentx = C_LOC( OpFM%u%momentx(1) )
   OpFM%u%c_obj%momenty_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%momenty = C_LOC( OpFM%u%momenty(1) )
   OpFM%u%c_obj%momentz_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%momentz = C_LOC( OpFM%u%momentz(1) )
   OpFM%u%c_obj%forceNodesChord_Len = OpFM%p%NnodesForce; OpFM%u%c_obj%forceNodesChord = C_LOC( OpFM%u%forceNodesChord(1) )
   if (InitInp%NumCtrl2SC > 0) then
      OpFM%u%c_obj%SuperController_Len = InitInp%NumCtrl2SC
      OpFM%u%c_obj%SuperController     = C_LOC( OpFM%u%SuperController(1) )
      OpFM%u%SuperController = 0.0_ReKi
   end if

      ! initialize the arrays:
   call OpFM_CreateActForceBladeTowerNodes(OpFM%p, ErrStat2, ErrMsg2) !Creates the blade and tower nodes in radial and tower height co-ordinates
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call OpFM_InterpolateForceNodesChord(initOut_AD, OpFM%p, OpFM%u, ErrStat2, ErrMsg2) !Interpolates the chord distribution to the force nodes
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   call OpFM_CreateActForceMotionsMesh( p_FAST, y_ED, y_BD, InitInp, OpFM, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   

      !............................................................................................
   ! Allocate arrays and set up mappings to point loads (for AD15 only):
      ! (bjj: note that normally I'd put these things in the FAST_ModuleMapType, but I don't want
      ! to add OpenFOAM integrations in the rest fo the code).
      !............................................................................................
   ! Allocate space for mapping data structures
   ALLOCATE( OpFM%m%ActForceLoads(OpFM%p%NMappings), OpFM%m%Line2_to_Line2_Loads(OpFM%p%NMappings), OpFM%m%Line2_to_Line2_Motions(OpFM%p%NMappings),STAT=ErrStat2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ALLOCATE( OpFM%m%ActForceLoadsPoints(OpFM%p%NMappings), OpFM%m%Line2_to_Point_Loads(OpFM%p%NMappings), OpFM%m%Line2_to_Point_Motions(OpFM%p%NMappings),STAT=ErrStat2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   do k=1,OpFM%p%NMappings
      call MeshCopy (  SrcMesh  = OpFM%m%ActForceMotions(k)  &
           , DestMesh = OpFM%m%ActForceLoads(k) &
           , CtrlCode = MESH_SIBLING          &
           , IOS      = COMPONENT_OUTPUT      &
           , Force    = .true.                &
           , Moment   = .true.                &
           , Orientation = .true.             &
           , TranslationDisp = .true.         &
           , TranslationVel  = .true.         &
           , RotationVel     = .true.         &
           , ErrStat  = ErrStat2              &
           , ErrMess  = ErrMsg2               )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      OpFM%m%ActForceLoads(k)%RemapFlag = .true.
      call MeshCopy (  SrcMesh  = OpFM%m%ActForceMotionsPoints(k)  &
           , DestMesh = OpFM%m%ActForceLoadsPoints(k) &
           , CtrlCode = MESH_SIBLING          &
           , IOS      = COMPONENT_OUTPUT      &
           , Force    = .true.                &
           , Moment   = .true.                &
           , Orientation = .true.             &
           , TranslationDisp = .true.         &
           , TranslationVel  = .true.         &
           , RotationVel     = .true.         &
           , ErrStat  = ErrStat2              &
           , ErrMess  = ErrMsg2               )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      OpFM%m%ActForceLoadsPoints(k)%RemapFlag = .true.
   end do
   
   ! create the mapping data structures:
   DO k=1,OpFM%p%NumBl
      IF (p_FAST%CompElast == Module_ED ) THEN
         call MeshMapCreate( y_ED%BladeLn2Mesh(k), OpFM%m%ActForceMotions(k), OpFM%m%Line2_to_Line2_Motions(k),  ErrStat2, ErrMsg2 );
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ELSEIF (p_FAST%CompElast == Module_BD ) THEN
         call MeshMapCreate( y_BD(k)%BldMotion, OpFM%m%ActForceMotions(k), OpFM%m%Line2_to_Line2_Motions(k),  ErrStat2, ErrMsg2 );
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
      call MeshMapCreate( y_AD%BladeLoad(k), OpFM%m%ActForceLoads(k), OpFM%m%Line2_to_Line2_Loads(k),  ErrStat2, ErrMsg2 );
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call MeshMapCreate( OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k),  ErrStat2, ErrMsg2 );
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call MeshMapCreate( OpFM%m%ActForceLoads(k), OpFM%m%ActForceLoadsPoints(k), OpFM%m%Line2_to_Point_Loads(k),  ErrStat2, ErrMsg2 );
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!      OpFM%m%ActForceLoads(k)%RemapFlag = .false.
   END DO
   
   do k=OpFM%p%NumBl+1,OpFM%p%NMappings
      call MeshMapCreate( y_ED%TowerLn2Mesh, OpFM%m%ActForceMotions(k), OpFM%m%Line2_to_Line2_Motions(k),  ErrStat2, ErrMsg2 );
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call MeshMapCreate( OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k),  ErrStat2, ErrMsg2 );
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      if ( y_AD%TowerLoad%nnodes > 0 ) then ! we can have an input mesh on the tower without having an output mesh.
         call MeshMapCreate( y_AD%TowerLoad, OpFM%m%ActForceLoads(k), OpFM%m%Line2_to_Line2_Loads(k),  ErrStat2, ErrMsg2 );
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         call MeshMapCreate( OpFM%m%ActForceLoads(k), OpFM%m%ActForceLoadsPoints(k), OpFM%m%Line2_to_Point_Loads(k),  ErrStat2, ErrMsg2 );
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!         OpFM%m%ActForceLoads(k)%RemapFlag = .false.
      end if
      
   end do
   
   call SetOpFMPositions(p_FAST, u_AD14, u_AD, y_ED, y_BD, OpFM, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   OpFM%u%fx = 0.0_ReKi
   OpFM%u%fy = 0.0_ReKi
   OpFM%u%fz = 0.0_ReKi
   
      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................
   CALL AllocPAry( OpFM%y%u, OpFM%p%NnodesVel, 'u', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%y%v, OpFM%p%NnodesVel, 'v', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%y%w, OpFM%p%NnodesVel, 'w', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (InitInp%NumSC2Ctrl > 0) then
      CALL AllocPAry( OpFM%y%SuperController, InitInp%NumSC2Ctrl, 'y%SuperController', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if

   IF (ErrStat >= AbortErrLev) RETURN

      ! make sure the C versions are synced with these arrays
   OpFM%y%c_obj%u_Len = OpFM%p%NnodesVel; OpFM%y%c_obj%u = C_LOC( OpFM%y%u(1) )
   OpFM%y%c_obj%v_Len = OpFM%p%NnodesVel; OpFM%y%c_obj%v = C_LOC( OpFM%y%v(1) )
   OpFM%y%c_obj%w_Len = OpFM%p%NnodesVel; OpFM%y%c_obj%w = C_LOC( OpFM%y%w(1) )

   if (InitInp%NumSC2Ctrl > 0) then
      OpFM%y%c_obj%SuperController_Len = InitInp%NumSC2Ctrl
      OpFM%y%c_obj%SuperController     = C_LOC( OpFM%y%SuperController(1) )
   end if



      !............................................................................................
      ! Define initialization-routine output (including writeOutput array) here:
      !............................................................................................

   CALL AllocAry( InitOut%WriteOutputHdr, 3, 'WriteOutputHdr', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry( InitOut%WriteOutputUnt, 3, 'WriteOutputUnt', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry( OpFM%y%WriteOutput, 3, 'WriteOutput', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN

   InitOut%WriteOutputHdr(1) = 'Wind1VelX'; InitOut%WriteOutputUnt(1) = '(m/s)'
   InitOut%WriteOutputHdr(2) = 'Wind1VelY'; InitOut%WriteOutputUnt(2) = '(m/s)'
   InitOut%WriteOutputHdr(3) = 'Wind1VelZ'; InitOut%WriteOutputUnt(3) = '(m/s)'
   OpFM%y%WriteOutput = 0.0_ReKi

   InitOut%Ver = OpFM_Ver

   RETURN

END SUBROUTINE Init_OpFM
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_SetInputs( p_FAST, p_AD14, u_AD14, y_AD14, u_AD, y_AD, y_ED, y_BD, y_SrvD, OpFM, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD14_ParameterType),       INTENT(IN)      :: p_AD14      ! The parameters from AeroDyn14 (for mesh transfer with improperly set meshes)
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD14_OutputType),          INTENT(IN)      :: y_AD14      ! The output meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN)      :: y_AD        ! The output meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        ! The outputs of the structural dynamics module
   TYPE(BD_OutputType),            INTENT(IN)      :: y_BD(:)     ! The outputs of the structural dynamics module   
   TYPE(SrvD_OutputType),          INTENT(IN)      :: y_SrvD      ! The outputs of the ServoDyn module (control)
   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'OpFM_SetInputs'


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! set the positions
   call SetOpFMPositions(p_FAST, u_AD14, u_AD, y_ED, y_BD, OpFM, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
      ! set the forces
   call SetOpFMForces(p_FAST, p_AD14, u_AD14, y_AD14, u_AD, y_AD, y_ED, y_BD, OpFM, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
  
      ! set SuperController inputs
   if (p_FAST%CompServo == Module_SrvD) then
      if (allocated(y_SrvD%SuperController).and. associated(OpFM%u%SuperController)) OpFM%u%SuperController = y_SrvD%SuperController
   end if


END SUBROUTINE OpFM_SetInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOpFMPositions(p_FAST, u_AD14, u_AD, y_ED, y_BD, OpFM, ErrStat, ErrMsg)

  TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data
  TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        ! The outputs of the structural dynamics module
   TYPE(BD_OutputType),            INTENT(IN)      :: y_BD(:)     ! The outputs of the structural dynamics module   

   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables:

   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetOpFMPositions'


   ErrStat  = ErrID_None
   ErrMsg = ''

   ! Do the Velocity (AeroDyn) nodes first
   !-------------------------------------------------------------------------------------------------
   Node = 1   ! displaced hub position 
   OpFM%u%pxVel(Node) = y_ED%HubPtMotion%Position(1,1) + y_ED%HubPtMotion%TranslationDisp(1,1)
   OpFM%u%pyVel(Node) = y_ED%HubPtMotion%Position(2,1) + y_ED%HubPtMotion%TranslationDisp(2,1)
   OpFM%u%pzVel(Node) = y_ED%HubPtMotion%Position(3,1) + y_ED%HubPtMotion%TranslationDisp(3,1)
   OpFM%u%pxDotVel(Node) = 0.0 !y_ED%HubPtMotion%TranslationVel(1,1)
   OpFM%u%pyDotVel(Node) = 0.0 !y_ED%HubPtMotion%TranslationVel(2,1)
   OpFM%u%pzDotVel(Node) = 0.0! y_ED%HubPtMotion%TranslationVel(3,1)
   
   ! blade nodes
   DO K = 1,SIZE(u_AD%BladeMotion)
      DO J = 1,u_AD%BladeMotion(k)%Nnodes
         Node = Node + 1
         OpFM%u%pxVel(Node) = u_AD%BladeMotion(k)%TranslationDisp(1,j) + u_AD%BladeMotion(k)%Position(1,j)
         OpFM%u%pyVel(Node) = u_AD%BladeMotion(k)%TranslationDisp(2,j) + u_AD%BladeMotion(k)%Position(2,j)
         OpFM%u%pzVel(Node) = u_AD%BladeMotion(k)%TranslationDisp(3,j) + u_AD%BladeMotion(k)%Position(3,j)
         OpFM%u%pxDotVel(Node) = u_AD%BladeMotion(k)%TranslationVel(1,j)
         OpFM%u%pyDotVel(Node) = u_AD%BladeMotion(k)%TranslationVel(2,j)
         OpFM%u%pzDotVel(Node) = u_AD%BladeMotion(k)%TranslationVel(3,j)         
         
      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
   END DO !K = 1,p%NumBl
   
   if (OpFM%p%NMappings .gt. OpFM%p%NumBl) then
      ! tower nodes
      DO J=1,u_AD%TowerMotion%nnodes
         Node = Node + 1
         OpFM%u%pxVel(Node) = u_AD%TowerMotion%TranslationDisp(1,J) + u_AD%TowerMotion%Position(1,J)
         OpFM%u%pyVel(Node) = u_AD%TowerMotion%TranslationDisp(2,J) + u_AD%TowerMotion%Position(2,J)
         OpFM%u%pzVel(Node) = u_AD%TowerMotion%TranslationDisp(3,J) + u_AD%TowerMotion%Position(3,J)
         OpFM%u%pxDotVel(Node) = u_AD%TowerMotion%TranslationVel(1,J)
         OpFM%u%pyDotVel(Node) = u_AD%TowerMotion%TranslationVel(2,J)
         OpFM%u%pzDotVel(Node) = u_AD%TowerMotion%TranslationVel(3,J)
      END DO
   end if
   
   ! Do the Actuator Force nodes now
   Node = 1   ! displaced hub position 
   OpFM%u%pxForce(Node) = OpFM%u%pxVel(Node)
   OpFM%u%pyForce(Node) = OpFM%u%pyVel(Node) 
   OpFM%u%pzForce(Node) = OpFM%u%pzVel(Node)
   OpFM%u%pxDotForce(Node) = OpFM%u%pxDotVel(Node)
   OpFM%u%pyDotForce(Node) = OpFM%u%pyDotVel(Node) 
   OpFM%u%pzDotForce(Node) = OpFM%u%pzDotVel(Node)
   OpFM%u%pOrientation((Node-1)*9 + 1) = y_ED%HubPtMotion%Orientation(1,1,1)
   OpFM%u%pOrientation((Node-1)*9 + 2) = y_ED%HubPtMotion%Orientation(2,1,1)
   OpFM%u%pOrientation((Node-1)*9 + 3) = y_ED%HubPtMotion%Orientation(3,1,1)   
   OpFM%u%pOrientation((Node-1)*9 + 4) = y_ED%HubPtMotion%Orientation(1,2,1)
   OpFM%u%pOrientation((Node-1)*9 + 5) = y_ED%HubPtMotion%Orientation(2,2,1)
   OpFM%u%pOrientation((Node-1)*9 + 6) = y_ED%HubPtMotion%Orientation(3,2,1)   
   OpFM%u%pOrientation((Node-1)*9 + 7) = y_ED%HubPtMotion%Orientation(1,3,1)
   OpFM%u%pOrientation((Node-1)*9 + 8) = y_ED%HubPtMotion%Orientation(2,3,1)
   OpFM%u%pOrientation((Node-1)*9 + 9) = y_ED%HubPtMotion%Orientation(3,3,1)   

   
   DO K = 1,OpFM%p%NumBl
      ! mesh mapping from line2 mesh to point mesh
      IF (p_FAST%CompElast == Module_ED ) THEN
         call Transfer_Line2_to_Line2( y_ED%BladeLn2Mesh(k), OpFM%m%ActForceMotions(k), OpFM%m%Line2_to_Line2_Motions(k), ErrStat2, ErrMsg2 )
         if(ErrStat2 /= 0) then
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         end if
      ELSEIF (p_FAST%CompElast == Module_BD ) THEN
         call Transfer_Line2_to_Line2( y_BD(k)%BldMotion, OpFM%m%ActForceMotions(k), OpFM%m%Line2_to_Line2_Motions(k), ErrStat2, ErrMsg2 )
         if(ErrStat2 /= 0) then
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         end if
      END IF
      call Transfer_Line2_to_Point( OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 )
      if(ErrStat2 /= 0) then
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      
      DO J = 1, OpFM%p%NnodesForceBlade
         Node = Node + 1
         OpFM%u%pxForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(1,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(1,J)
         OpFM%u%pyForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(2,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(2,J)
         OpFM%u%pzForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(3,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(3,J)
         OpFM%u%pxDotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(1,J)
         OpFM%u%pyDotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(2,J)
         OpFM%u%pzDotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(3,J)         
         OpFM%u%pOrientation((Node-1)*9 + 1) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 2) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,1,J)         
         OpFM%u%pOrientation((Node-1)*9 + 3) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 4) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 5) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,2,J)         
         OpFM%u%pOrientation((Node-1)*9 + 6) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 7) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,3,J)
         OpFM%u%pOrientation((Node-1)*9 + 8) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,3,J)         
         OpFM%u%pOrientation((Node-1)*9 + 9) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,3,J)
      END DO
      
   END DO
   
   DO K = OpFM%p%NumBl+1,OpFM%p%NMappings

      call Transfer_Line2_to_Line2( y_ED%TowerLn2Mesh, OpFM%m%ActForceMotions(k), OpFM%m%Line2_to_Line2_Motions(k), ErrStat2, ErrMsg2 )
      if(ErrStat2 /= 0) then
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      call Transfer_Line2_to_Point( OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 )
      if(ErrStat2 /= 0) then
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      
      DO J=1,OpFM%p%NnodesForceTower
         Node = Node + 1
         OpFM%u%pxForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(1,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(1,J)
         OpFM%u%pyForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(2,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(2,J)
         OpFM%u%pzForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(3,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(3,J)
         OpFM%u%pxDotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(1,J)
         OpFM%u%pyDotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(2,J)
         OpFM%u%pzDotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(3,J)
         OpFM%u%pOrientation((Node-1)*9 + 1) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 2) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,1,J)         
         OpFM%u%pOrientation((Node-1)*9 + 3) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 4) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 5) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,2,J)         
         OpFM%u%pOrientation((Node-1)*9 + 6) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 7) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,3,J)
         OpFM%u%pOrientation((Node-1)*9 + 8) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,3,J)         
         OpFM%u%pOrientation((Node-1)*9 + 9) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,3,J)

      END DO
      
   END DO
   

END SUBROUTINE SetOpFMPositions
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOpFMForces(p_FAST, p_AD14, u_AD14, y_AD14, u_AD, y_AD, y_ED, y_BD, OpFM, ErrStat, ErrMsg)

   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   TYPE(AD14_ParameterType),       INTENT(IN)      :: p_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD14_OutputType),          INTENT(IN)      :: y_AD14      ! The output meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN)      :: y_AD        ! The output meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        ! The outputs of the structural dynamics module
   TYPE(BD_OutputType),            INTENT(IN)      :: y_BD(:)     ! The outputs of the structural dynamics module   
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data
   !TYPE(FAST_MiscVarType),         INTENT(IN   )   :: m_FAST      ! misc FAST data, including inputs from external codes like Simulink
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables:
   REAL(ReKi )                                     :: factor      ! scaling factor to get normalized forces for OpenFOAM
   REAL(ReKi)                                      :: dRforceNodes ! Uniform distance between two consecutive blade force nodes
   REAL(ReKi)                                      :: dHforceNodes ! Uniform distance between two consecutive tower force nodes

   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh
#ifdef DEBUG_OPENFOAM   
   INTEGER(IntKi)                                  :: actForcesFile, aerodynForcesFile ! Unit numbers for files containing actuator forces and aerodyn forces
#endif
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetOpFMForces'

   ErrStat = ErrID_None
   ErrMsg  = ''

   !-------------------------------------------------------------------------------------------------
   Node = 1   ! undisplaced hub position  (no aerodynamics computed here)
   OpFM%u%fx(Node) = 0.0_ReKi
   OpFM%u%fy(Node) = 0.0_ReKi
   OpFM%u%fz(Node) = 0.0_ReKi

   !.......................
   ! blade nodes
   !.......................

#ifdef DEBUG_OPENFOAM   
   CALL GetNewUnit( aerodynForcesFile )
   open(unit=aerodynForcesFile,file='fast_aerodyn_velocity_forces.csv')
   write(aerodynForcesFile,*) '#x, y, z, fx, fy, fz'

   CALL GetNewUnit( actForcesFile )
   open(unit=actForcesFile,file='fast_actuator_forces.csv')
   write(actForcesFile,*) '#x, y, z, fx, fy, fz'
#endif

   DO K = 1,OpFM%p%NumBl
      
#ifdef DEBUG_OPENFOAM   
      DO J = 1,u_AD%BladeMotion(k)%NNodes
        write(aerodynForcesFile,*) u_AD%BladeMotion(k)%TranslationDisp(1,j) + u_AD%BladeMotion(k)%Position(1,j), ', ', u_AD%BladeMotion(k)%TranslationDisp(2,j) + u_AD%BladeMotion(k)%Position(2,j), ', ', u_AD%BladeMotion(k)%TranslationDisp(3,j) + u_AD%BladeMotion(k)%Position(3,j), ', ', OpFM%y%u(1 + (k-1)*u_AD%BladeMotion(k)%NNodes + j), ', ', OpFM%y%v(1 + (k-1)*u_AD%BladeMotion(k)%NNodes + j), ', ', OpFM%y%w(1 + (k-1)*u_AD%BladeMotion(k)%NNodes + j), ', ', y_AD%BladeLoad(k)%Force(1,j), ', ', y_AD%BladeLoad(k)%Force(2,j), ', ', y_AD%BladeLoad(k)%Force(2,j)
      END DO
#endif

      call Transfer_Line2_to_Line2( y_AD%BladeLoad(k), OpFM%m%ActForceLoads(k), OpFM%m%Line2_to_Line2_Loads(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), OpFM%m%ActForceMotions(k) )
      if(ErrStat2 /= 0) then
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      call Transfer_Line2_to_Point( OpFM%m%ActForceLoads(k), OpFM%m%ActForceLoadsPoints(k), OpFM%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k) )
      if(ErrStat2 /= 0) then
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if

      DO J = 1, OpFM%p%NnodesForceBlade
         Node = Node + 1
         OpFM%u%fx(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(1,j)
         OpFM%u%fy(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(2,j)
         OpFM%u%fz(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(3,j)
         OpFM%u%momentx(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(1,j)
         OpFM%u%momenty(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(2,j)
         OpFM%u%momentz(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM   
         write(actForcesFile,*) OpFM%u%pxForce(Node), ', ', OpFM%u%pyForce(Node), ', ', OpFM%u%pzForce(Node), ', ', OpFM%u%fx(Node), ', ', OpFM%u%fy(Node), ', ', OpFM%u%fz(Node), ', '
#endif

      END DO 

   END DO !K = 1,OpFM%p%NumBl

   !.......................
   ! tower nodes
   !.......................

   ! mesh mapping from line2 mesh to point mesh
   DO K = OpFM%p%NumBl+1,OpFM%p%NMappings

#ifdef DEBUG_OPENFOAM
   DO J = 1,u_AD%TowerMotion%NNodes
      write(aerodynForcesFile,*) u_AD%TowerMotion%TranslationDisp(1,j) + u_AD%TowerMotion%Position(1,j), ', ', u_AD%TowerMotion%TranslationDisp(2,j) + u_AD%TowerMotion%Position(2,j), ', ', u_AD%TowerMotion%TranslationDisp(3,j) + u_AD%TowerMotion%Position(3,j), ', ', y_AD%TowerLoad%Force(1,j), ', ', y_AD%TowerLoad%Force(2,j), ', ', y_AD%TowerLoad%Force(2,j)
   END DO
#endif

   call Transfer_Line2_to_Line2( y_AD%TowerLoad, OpFM%m%ActForceLoads(k), OpFM%m%Line2_to_Line2_Loads(k), ErrStat2, ErrMsg2, u_AD%TowerMotion, OpFM%m%ActForceMotions(k) )
   if(ErrStat .ne. 0) then
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
   call Transfer_Line2_to_Point( OpFM%m%ActForceLoads(k), OpFM%m%ActForceLoadsPoints(k), OpFM%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k) )
   if(ErrStat .ne. 0) then
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
   
   DO J=1,OpFM%p%NnodesForceTower
      Node = Node + 1
      OpFM%u%fx(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(1,j)
      OpFM%u%fy(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(2,j)
      OpFM%u%fz(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(3,j)
      OpFM%u%momentx(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(1,j)
      OpFM%u%momenty(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(2,j)
      OpFM%u%momentz(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM   
      write(actForcesFile,*) OpFM%u%pxForce(Node), ', ', OpFM%u%pyForce(Node), ', ', OpFM%u%pzForce(Node), ', ', OpFM%u%fx(Node), ', ', OpFM%u%fy(Node), ', ', OpFM%u%fz(Node), ', '
#endif
   END DO

#ifdef DEBUG_OPENFOAM   
   close(aerodynForcesFile)
   close(actForcesFile)
#endif

   END DO

END SUBROUTINE SetOpFMForces
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_SetWriteOutput( OpFM )
!..................................................................................................................................

   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module

   ! set the hub-height wind speeds
   IF ( ALLOCATED( OpFM%y%WriteOutput ) ) THEN
      IF ( ASSOCIATED( OpFM%y%u ) ) then
         OpFM%y%WriteOutput(1) = OpFM%y%u(1)
         OpFM%y%WriteOutput(2) = OpFM%y%v(1)
         OpFM%y%WriteOutput(3) = OpFM%y%w(1)
      END IF
   END IF



END SUBROUTINE OpFM_SetWriteOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_CreateActForceMotionsMesh( p_FAST, y_ED, y_BD, InitIn_OpFM, OpFM, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(ED_OutputType),             INTENT(IN)     :: y_ED        ! The outputs of the structural dynamics module
   TYPE(BD_OutputType),             INTENT(IN)     :: y_BD(:)     ! The outputs of the structural dynamics module   
   TYPE(OpFM_InitInputType),        INTENT(IN   )  :: InitIn_OpFM ! InitInp data for the OpenFOAM integration module
   TYPE(OpenFOAM_Data),             INTENT(INOUT)  :: OpFM        ! data for the OpenFOAM integration module
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   TYPE(MeshType) , DIMENSION(:), ALLOCATABLE      :: tmpActForceMotionsMesh   !< temporary mesh for interpolating orientation to actuator force points [-]
   INTEGER(IntKi)                                  :: k          ! blade loop counter
   INTEGER(IntKi)                                  :: i,j          ! node counter

   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'OpFM_CreateActForceMotionsMesh'

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Allocate space for mapping data structures
      ALLOCATE(tmpActForceMotionsMesh(OpFM%p%NMappings) , STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal, 'Error allocating force nodes mesh mapping types', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
      CALL OpFM_CreateTmpActForceMotionsMesh( p_FAST, y_ED, y_BD, OpFM%p, InitIn_OpFM, tmpActForceMotionsMesh, ErrStat2, ErrMsg2 )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           if (errStat >= AbortErrLev) return

      ALLOCATE(OpFM%m%ActForceMotions(OpFM%p%NMappings), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal, 'Error allocating force nodes mesh', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
      ALLOCATE(OpFM%m%ActForceMotionsPoints(OpFM%p%NMappings), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal, 'Error allocating force nodes mesh', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
      DO k=1,OpFM%p%NumBl
         call MeshCreate ( BlankMesh = OpFM%m%ActForceMotions(k)         &
                          ,IOS       = COMPONENT_INPUT             &
                          ,Nnodes    = OpFM%p%NnodesForceBlade &
!                          ,Force           = .false.        &
!                          ,Moment          = .false.        &
                          ,Orientation     = .true.         &
                          ,TranslationDisp = .true.         &
                          ,TranslationVel  = .true.         &
                          ,RotationVel     = .true.         &
                          ,ErrStat   = ErrStat2                    &
                          ,ErrMess   = ErrMsg2                     &                          
                         )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
               OpFM%m%ActForceMotions(k)%RemapFlag = .false. 

               call MeshCreate ( BlankMesh = OpFM%m%ActForceMotionsPoints(k)         &
                          ,IOS       = COMPONENT_INPUT             &
                          ,Nnodes    = OpFM%p%NnodesForceBlade &
!                          ,Force           = .false.        &
!                          ,Moment          = .false.        &
                          ,Orientation     = .true.         &
                          ,TranslationDisp = .true.         &
                          ,TranslationVel  = .true.         &
                          ,RotationVel     = .true.         &
                          ,ErrStat   = ErrStat2                    &
                          ,ErrMess   = ErrMsg2                     &                          
                         )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
               OpFM%m%ActForceMotions(k)%RemapFlag = .false. 

         do j=1,OpFM%p%NnodesForceBlade
            call MeshPositionNode(OpFM%m%ActForceMotions(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                                  orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

            call MeshPositionNode(OpFM%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                                  orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            call MeshConstructElement(OpFM%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         end do !j

        ! create elements:      
        DO J = 2,OpFM%p%NnodesForceBlade
           call MeshConstructElement ( Mesh      = OpFM%m%ActForceMotions(k)  &
                                                 , Xelement = ELEMENT_LINE2      &
                                                 , P1       = J-1                &   ! node1 number
                                                 , P2       = J                  &   ! node2 number
                                                 , ErrStat  = ErrStat2           &
                                                 , ErrMess  = ErrMsg2            )
           call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
        END DO ! J (blade nodes)
        call MeshCommit(OpFM%m%ActForceMotions(k), errStat2, errMsg2 )
        call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
        if (errStat >= AbortErrLev) return
        call MeshCommit(OpFM%m%ActForceMotionsPoints(k), errStat2, errMsg2 )
        call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
        if (errStat >= AbortErrLev) return
      END DO

      DO k=OpFM%p%NumBl+1,OpFM%p%NMappings !Tower if present
         call MeshCreate ( BlankMesh = OpFM%m%ActForceMotions(k)         &
                          ,IOS       = COMPONENT_INPUT             &
                          ,Nnodes    = OpFM%p%NnodesForceTower &
                          ,Orientation     = .true.         &
                          ,TranslationDisp = .true.         &
                          ,TranslationVel  = .true.         &
                          ,RotationVel     = .true.         &
                          ,ErrStat   = ErrStat2                    &
                          ,ErrMess   = ErrMsg2                     &
                         )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
               OpFM%m%ActForceMotions(k)%RemapFlag = .false. 

         call MeshCreate ( BlankMesh = OpFM%m%ActForceMotionsPoints(k)         &
                          ,IOS       = COMPONENT_INPUT             &
                          ,Nnodes    = OpFM%p%NnodesForceTower &
                          ,Orientation     = .true.         &
                          ,TranslationDisp = .true.         &
                          ,TranslationVel  = .true.         &
                          ,RotationVel     = .true.         &
                          ,ErrStat   = ErrStat2                    &
                          ,ErrMess   = ErrMsg2                     &
                         )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
               OpFM%m%ActForceMotionsPoints(k)%RemapFlag = .false. 

         do j=1,OpFM%p%NnodesForceTower
            call MeshPositionNode(OpFM%m%ActForceMotions(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                                  orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

            call MeshPositionNode(OpFM%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                                  orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            call MeshConstructElement(OpFM%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         end do !j
        ! create elements:      
        DO J = 2,OpFM%p%NnodesForceTower
           call MeshConstructElement ( Mesh      = OpFM%m%ActForceMotions(k)  &
                                                 , Xelement = ELEMENT_LINE2      &
                                                 , P1       = J-1                &   ! node1 number
                                                 , P2       = J                  &   ! node2 number
                                                 , ErrStat  = ErrStat2           &
                                                 , ErrMess  = ErrMsg2            )
           call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
        END DO ! J (tower nodes)
         
         call MeshCommit(OpFM%m%ActForceMotions(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            if (errStat >= AbortErrLev) return
         call MeshCommit(OpFM%m%ActForceMotionsPoints(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            if (errStat >= AbortErrLev) return
      END DO

      DO k=1,OpFM%p%NMappings
         call MeshDestroy ( tmpActForceMotionsMesh(k), ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      DEALLOCATE(tmpActForceMotionsMesh)

END SUBROUTINE OpFM_CreateActForceMotionsMesh
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_CreateTmpActForceMotionsMesh( p_FAST, y_ED, y_BD, p_OpFM, InitIn_OpFM, tmpActForceMotions, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST        ! Parameters for the glue code
   TYPE(ED_OutputType),             INTENT(IN   )  :: y_ED          ! The outputs of the structural dynamics module
   TYPE(BD_OutputType),             INTENT(IN   )  :: y_BD(:)       ! The outputs of the structural dynamics module   
   TYPE(OpFM_ParameterType),        INTENT(IN   )  :: p_OpFM        ! data for the OpenFOAM integration module
   TYPE(OpFM_InitInputType),        INTENT(IN   )  :: InitIn_OpFM   ! InitInp data for the OpenFOAM integration module
   TYPE(MeshType),                  INTENT(INOUT)  :: tmpActForceMotions(:) ! temporary mesh to create the actuator force nodes
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   TYPE(MeshMapType) , DIMENSION(:), ALLOCATABLE   :: tmp_line2_to_point_Motions    !< mapping data structure to convert orientation of structural nodes to actuator force nodes [-]
   TYPE(MeshType) , DIMENSION(:), ALLOCATABLE      :: tmp_StructModelMesh   !< temporary mesh copying Structural model mesh
   REAL(ReKi), DIMENSION(:,:), ALLOCATABLE         :: forceNodePositions  ! new positions for the force actuator nodes
   INTEGER(IntKi)                                  :: NumBl      ! number of blades
   INTEGER(IntKi)                                  :: k          ! blade loop counter
   INTEGER(IntKi)                                  :: i,j          ! node counter
   INTEGER(IntKi)                                  :: sBldNodes    ! Temporary variable to keep track of number of struct blade nodes

   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'OpFM_CreateTmpActForceMotionsMesh'

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Make a copy of the Structural model mesh with the reference orientation set to zero
   ALLOCATE(tmp_StructModelMesh(p_OpFM%NMappings) , STAT=ErrStat2)
   IF (ErrStat2 /= 0) THEN
      CALL SetErrStat(ErrID_Fatal, 'Error allocating temporary copy of ElastoDyn mesh type', ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF
   CALL CreateTmpStructModelMesh(p_FAST, y_ED, y_BD, p_OpFM, tmp_StructModelMesh, ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF (ErrStat >= AbortErrLev) RETURN
 
   ! Allocate space for mapping data structures
   ALLOCATE( tmp_line2_to_point_Motions(p_OpFM%NMappings),STAT=ErrStat2)
   IF (ErrStat2 /= 0) THEN
      CALL SetErrStat(ErrID_Fatal, 'Error allocating temporary actuator force mesh mapping types', ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF

   ! create meshes to map:
   sBldNodes = 0
   ALLOCATE(forceNodePositions(3,p_OpFM%NnodesForceBlade)) ! Allocate space to create new positions
   DO k=1,p_OpFM%NumBl
      call MeshCreate ( BlankMesh = tmpActForceMotions(k)         &
           ,IOS       = COMPONENT_INPUT             &
           ,Nnodes    = p_OpFM%NnodesForceBlade &
           ,ErrStat   = ErrStat2                    &
           ,ErrMess   = ErrMsg2                     &
           ,force     = .false.                     &
           ,moment    = .false.                     &
           ,orientation = .true.                    &
           ,translationvel  = .true.                &
           ,rotationvel     = .true.                &
           )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
      tmpActForceMotions(k)%RemapFlag = .false.
      call CalcForceActuatorPositionsBlade(p_OpFM, InitIn_OpFM%nStructBldEtaNodes(k), InitIn_OpFM%StructBldEtaNodes(sBldNodes+1:sBldNodes+InitIn_OpFM%nStructBldEtaNodes(k)), tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      do j=1,p_OpFM%NnodesForceBlade
         call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
         call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
      
      call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return
      sBldNodes = sBldNodes + InitIn_OpFM%nStructBldEtaNodes(k)
   end do
   DEALLOCATE(forceNodePositions) ! Free space

   ALLOCATE(forceNodePositions(3,p_OpFM%NnodesForceTower)) ! Allocate space to create new positions
   DO k=p_OpFM%NumBl+1,p_OpFM%NMappings   
      call CalcForceActuatorPositionsTower(p_OpFM, SIZE(InitIn_OpFM%StructTwrEtaNodes), InitIn_OpFM%StructTwrEtaNodes, tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      call MeshCreate ( BlankMesh = tmpActForceMotions(k)        &
           ,IOS       = COMPONENT_INPUT             &
           ,Nnodes    = p_OpFM%NnodesForceTower &
           ,ErrStat   = ErrStat2                    &
           ,ErrMess   = ErrMsg2                     &
           ,force     = .false.                     &
           ,moment    = .false.                     &
           ,orientation = .true.                    &
           ,translationvel  = .true.                &
           ,rotationvel     = .true.                &
           )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      tmpActForceMotions(k)%RemapFlag = .false.
      do j=1,p_OpFM%NnodesForceTower
         call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
         call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
      
      call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if (errStat >= AbortErrLev) return
   END DO
   DEALLOCATE(forceNodePositions) ! Free space
   
   ! create the mapping data structures:
   DO k=1,p_OpFM%NumBl
      call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 );
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO
   
   DO k=p_OpFM%NumBl+1,p_OpFM%NMappings
      call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 );
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO
   
   ! Map the orientation
   DO K = 1,p_OpFM%NMappings
      ! mesh mapping from line2 mesh to point mesh
      call Transfer_Line2_to_Point( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k), ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   END DO

   DO k=1,p_OpFM%NMappings
      call MeshDestroy ( tmp_StructModelMesh(k), ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call MeshMapDestroy ( tmp_line2_to_point_Motions(k), ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO
   DEALLOCATE(tmp_StructModelMesh)
   DEALLOCATE(tmp_line2_to_point_Motions)   

   RETURN

END SUBROUTINE OpFM_CreateTmpActForceMotionsMesh
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateTmpStructModelMesh(p_FAST, y_ED, y_BD, p_OpFM, tmpStructModelMesh, ErrStat, ErrMsg )

  TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
  TYPE(ED_OutputType),             INTENT(IN   )  :: y_ED        ! The outputs of the structural dynamics module
  TYPE(BD_OutputType),             INTENT(IN   )  :: y_BD(:)     ! The outputs of the structural dynamics module  
  TYPE(OpFM_ParameterType),        INTENT(IN   )  :: p_OpFM      ! Parameters of the OpenFOAM integration module
  TYPE(MeshType),                  INTENT(INOUT)  :: tmpStructModelMesh(:) ! temporary copy of structural model mesh
  INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
  CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


  !Local variables
  INTEGER(IntKi)                                  :: nNodesStructModel ! Number of nodes (tower/blade) in the structural model mesh

  INTEGER(IntKi)                                  :: i,j          ! node counter
  INTEGER(IntKi)                                  :: k            ! blade counter
  INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

  CHARACTER(*),   PARAMETER                       :: RoutineName = 'CreateTmpStructModelMesh'
  

  IF (p_FAST%CompElast == Module_ED ) THEN

     DO K = 1,p_OpFM%NumBl

        nNodesStructModel = SIZE(y_ED%BladeLn2Mesh(K)%position(1,:))
         
        CALL MeshCreate( BlankMesh       = tmpStructModelMesh(K)  &
                       , NNodes          = nNodesStructModel      &
                       , IOS             = COMPONENT_OUTPUT       &
                       , Orientation     = .TRUE.                 &
                       , ErrStat         = ErrStat2               &
                       , ErrMess         = ErrMsg2                )
        IF (ErrStat >= AbortErrLev) RETURN

         tmpStructModelMesh(K)%RemapFlag = .false.        
         !For some reason, ElastoDyn keeps the last point as the blade/tower root
         CALL MeshPositionNode ( tmpStructModelMesh(K), 1, y_ED%BladeLn2Mesh(K)%Position(:,nNodesStructModel), ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         DO J = 1,nNodesStructModel-1
           CALL MeshPositionNode ( tmpStructModelMesh(K), J+1, y_ED%BladeLn2Mesh(K)%Position(:,J), ErrStat2, ErrMsg2 )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
        
        ! create elements:      
        DO J = 2,nNodesStructModel
           
           CALL MeshConstructElement ( Mesh      = tmpStructModelMesh(K)  &
                                                 , Xelement = ELEMENT_LINE2      &
                                                 , P1       = J-1                &   ! node1 number
                                                 , P2       = J                  &   ! node2 number
                                                 , ErrStat  = ErrStat2           &
                                                 , ErrMess  = ErrMsg2            )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        END DO ! J (blade nodes)
        
        ! that's our entire mesh:
        CALL MeshCommit ( tmpStructModelMesh(K), ErrStat2, ErrMsg2 )
        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        ! Copy the orientation
        tmpStructModelMesh(K)%Orientation(:,:,1) = y_ED%BladeLn2Mesh(k)%RefOrientation(:,:,nNodesStructModel)
        DO J=1,nNodesStructModel-1
           tmpStructModelMesh(K)%Orientation(:,:,J+1) = y_ED%BladeLn2Mesh(K)%RefOrientation(:,:,J)           
        END DO
        
     END DO

  ELSEIF (p_FAST%CompElast == Module_BD ) THEN

     DO K = 1,p_OpFM%NumBl

        nNodesStructModel = SIZE(y_BD(K)%BldMotion%position(1,:))
         
        CALL MeshCreate( BlankMesh       = tmpStructModelMesh(K)  &
                       , NNodes          = nNodesStructModel      &
                       , IOS             = COMPONENT_OUTPUT       &
                       , Orientation     = .TRUE.                 &
                       , ErrStat         = ErrStat2               &
                       , ErrMess         = ErrMsg2                )
        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        IF (ErrStat >= AbortErrLev) RETURN

         tmpStructModelMesh(K)%RemapFlag = .false.        
         !For some reason, BeamDyn keeps the last point as the blade root
         DO J = 1,nNodesStructModel
           CALL MeshPositionNode ( tmpStructModelMesh(K), J, y_BD(K)%BldMotion%Position(:,J), ErrStat2, ErrMsg2 )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        END DO
        
        ! create elements:      
        DO J = 2,nNodesStructModel
           
           CALL MeshConstructElement ( Mesh      = tmpStructModelMesh(K)  &
                                                 , Xelement = ELEMENT_LINE2      &
                                                 , P1       = J-1                &   ! node1 number
                                                 , P2       = J                  &   ! node2 number
                                                 , ErrStat  = ErrStat2           &
                                                 , ErrMess  = ErrMsg2            )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        END DO ! J (blade nodes)
        
        ! that's our entire mesh:
        CALL MeshCommit ( tmpStructModelMesh(K), ErrStat2, ErrMsg2 )
        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
        ! Copy the orientation
        tmpStructModelMesh(K)%Orientation(:,:,1) = y_BD(K)%BldMotion%RefOrientation(:,:,nNodesStructModel)
        DO J=1,nNodesStructModel-1
           tmpStructModelMesh(K)%Orientation(:,:,J+1) = y_BD(K)%BldMotion%RefOrientation(:,:,J)           
        END DO
        
     END DO     

  END IF
     
     
  DO K = p_OpFM%NumBl+1, p_OpFM%NMappings
     
     nNodesStructModel = SIZE(y_ED%TowerLn2Mesh%position(1,:))
     
     CALL MeshCreate( BlankMesh       = tmpStructModelMesh(K)      &
          , NNodes          = nNodesStructModel          &
          , IOS             = COMPONENT_OUTPUT           &
          , Orientation     = .TRUE.                 &
          , ErrStat         = ErrStat2               &
          , ErrMess         = ErrMsg2                )
     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
     
     tmpStructModelMesh(K)%RemapFlag = .false.        
     !For some reason, ElastoDyn keeps the last point as the blade/tower root        
     CALL MeshPositionNode ( tmpStructModelMesh(K), 1, y_ED%TowerLn2Mesh%Position(:,nNodesStructModel), ErrStat2, ErrMsg2 )
     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
     DO J = 1,nNodesStructModel-1
        CALL MeshPositionNode ( tmpStructModelMesh(K), J+1, y_ED%TowerLn2Mesh%Position(:,J), ErrStat2, ErrMsg2 )
        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
     END DO
     
     ! create elements:      
     DO J = 2,nNodesStructModel
        
        CALL MeshConstructElement ( Mesh      = tmpStructModelMesh(K)  &
             , Xelement = ELEMENT_LINE2      &
             , P1       = J-1                &   ! node1 number
             , P2       = J                  &   ! node2 number
             , ErrStat  = ErrStat2           &
             , ErrMess  = ErrMsg2            )
        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
     END DO ! J (blade nodes)
     
     ! that's our entire mesh:
     CALL MeshCommit ( tmpStructModelMesh(K), ErrStat2, ErrMsg2 )
     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
     
     ! Copy the orientation
     tmpStructModelMesh(K)%Orientation(:,:,1) = y_ED%TowerLn2Mesh%RefOrientation(:,:,nNodesStructModel)
     DO J=1,nNodesStructModel-1
        tmpStructModelMesh(K)%Orientation(:,:,J+1) = y_ED%TowerLn2Mesh%RefOrientation(:,:,J)           
     END DO
     
  END DO     

  RETURN 
END SUBROUTINE CreateTmpStructModelMesh
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcForceActuatorPositionsBlade(p_OpFM, nBldEtaNodes, sBldEtaNodes, structPositions, forceNodePositions, ErrStat2, ErrMsg2)

  TYPE(OpFM_ParameterType), INTENT(IN )  :: p_OpFM        ! data for the OpenFOAM integration module
  INTEGER(IntKi),        INTENT(IN)      :: nBldEtaNodes  ! Number of structural model nodes
  REAL(ReKi),            INTENT(IN)      :: sBldEtaNodes(:) ! The non-dimensional co-ordinates at which the structural model positions are defined.
  REAL(ReKi),   POINTER, INTENT(IN)      :: structPositions(:,:)     ! structural model positions
  REAL(ReKi),            INTENT(INOUT)   :: forceNodePositions(:,:)  ! Array to store the newly created positions
  INTEGER(IntKi)                         :: ErrStat2    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                   :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None


  !Local variables
  INTEGER(IntKi)                         :: nStructNodes    ! Number of structural model nodes

  INTEGER(IntKI)                         :: i,j,k        ! Loop variables
  INTEGER(IntKI)                         :: jLower    ! Index of the struct node just smaller than the force node
  REAL(ReKi)                             :: rInterp      ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes

  CHARACTER(*),   PARAMETER                       :: RoutineName = 'CalcForceActuatorPositionsBlade'
  
  ! Now calculate the positions of the force nodes based on interpolation
  forceNodePositions(:,1) = structPositions(:,1)
  DO I=2,p_OpFM%NnodesForceBlade-1 ! Calculate the position of the force nodes
     jLower=1
     do while ( ( (sBldEtaNodes(jLower) - p_OpFM%forceBldEtaNodes(I))*(sBldEtaNodes(jLower+1) - p_OpFM%forceBldEtaNodes(I)) .gt. 0) .and. (jLower .lt. nBldEtaNodes) )
        jLower = jLower + 1
     end do
     rInterp =  (p_OpFM%forceBldEtaNodes(I) - sBldEtaNodes(jLower))/(sBldEtaNodes(jLower+1)-sBldEtaNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
     forceNodePositions(:,I) = structPositions(:,jLower) + rInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
  END DO
  forceNodePositions(:,p_OpFM%NnodesForceBlade) = structPositions(:,nBldEtaNodes)

  RETURN

END SUBROUTINE CalcForceActuatorPositionsBlade
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcForceActuatorPositionsTower(p_OpFM, nTwrEtaNodes, sTwrEtaNodes, structPositions, forceNodePositions, ErrStat2, ErrMsg2)

  TYPE(OpFM_ParameterType), INTENT(IN )  :: p_OpFM        ! data for the OpenFOAM integration module 
  INTEGER(IntKi),           INTENT(IN)   :: nTwrEtaNodes  ! Number of tower structural model nodes 
  REAL(ReKi), POINTER,      INTENT(IN)   :: sTwrEtaNodes(:) ! Location of tower structural model nodes in [0-1] coordinates
  REAL(ReKi),   POINTER                  :: structPositions(:,:)     ! structural model positions
  REAL(ReKi),             INTENT(INOUT)  :: forceNodePositions(:,:)  ! Array to store the newly created positions
  INTEGER(IntKi)                         :: ErrStat2    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                   :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None


  !Local variables
  INTEGER(IntKI)                         :: i,j,k        ! Loop variables
  INTEGER(IntKI)                         :: jLower    ! Index of the struct node just smaller than the force node
  REAL(ReKi)                             :: hInterp      ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes

  CHARACTER(*),   PARAMETER                       :: RoutineName = 'CalcForceActuatorPositionsTower'
  
  ! Now calculate the positions of the force nodes based on interpolation
  forceNodePositions(:,1) = structPositions(:,1)
  DO I=2,p_OpFM%NnodesForceTower-1 ! Calculate the position of the force nodes
     jLower=1
     do while ( ((sTwrEtaNodes(jLower) - p_OpFM%forceTwrEtaNodes(I))*(sTwrEtaNodes(jLower+1) - p_OpFM%forceTwrEtaNodes(I)) .gt. 0) .and. (jLower .lt. nTwrEtaNodes))
        jLower = jLower + 1
     end do
     hInterp =  (p_OpFM%forceTwrEtaNodes(I) - sTwrEtaNodes(jLower))/(sTwrEtaNodes(jLower+1)-sTwrEtaNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
     forceNodePositions(:,I) = structPositions(:,jLower) + hInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
  END DO
  forceNodePositions(:,p_OpFM%NnodesForceTower) = structPositions(:,nTwrEtaNodes)

  RETURN

END SUBROUTINE CalcForceActuatorPositionsTower

SUBROUTINE OpFM_CreateActForceBladeTowerNodes(p_OpFM, ErrStat, ErrMsg) 
!Creates the blade and tower nodes in blade local and tower local co-ordinates

  TYPE(OpFM_ParameterType), INTENT(INOUT)  :: p_OpFM        ! data for the OpenFOAM integration module
  INTEGER(IntKi)                         :: ErrStat    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                   :: ErrMsg     ! temporary Error message if ErrStat /= ErrID_None

  !Local variables
  REAL(ReKi)                             :: dEtaForceNodes ! Uniform distance between two consecutive force nodes for the blade
  REAL(ReKi)                             :: dHForceNodes ! Uniform distance between two consecutive force nodes for the tower
  INTEGER(IntKI)                         :: i            ! Loop variables
  INTEGER(IntKi)                         :: ErrStat2    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                   :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

  CHARACTER(*),   PARAMETER                       :: RoutineName = 'OpFM_CreateActForceBladeTowerNodes'
  

  ! Line2 to Line2 mapping expects the destination mesh to be smaller than the source mesh for deformation mapping and larger than the source mesh for load mapping. This forces me to create nodes at the very ends of the blade and tower.

  !Uniform distribution for now. Change to desired distribution here
  !Do the blade first
  allocate(p_OpFM%forceBldEtaNodes(p_OpFM%NnodesForceBlade), stat=errStat2)
  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
  
  dEtaForceNodes = 1.0/(p_OpFM%NnodesForceBlade-1)
  do i=1,p_OpFM%NnodesForceBlade-1
     p_OpFM%forceBldEtaNodes(i) = (i-1)*dEtaForceNodes 
  end do
  p_OpFM%forceBldEtaNodes(p_OpFM%NnodesForceBlade) = 1.0


  if (p_OpFM%NMappings .gt. p_OpFM%NumBl) then
     !Do the tower now
     allocate(p_OpFM%forceTwrEtaNodes(p_OpFM%NnodesForceTower), stat=errStat2)
     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
     dHforceNodes = 1.0/(p_OpFM%NnodesForceTower-1)
     do i=1,p_OpFM%NnodesForceTower-1
        p_OpFM%forceTwrEtaNodes(i) = (i-1)*dHForceNodes
     end do
     p_OpFM%forceTwrEtaNodes(p_OpFM%NnodesForceTower) = 1.0
  end if

  return

END SUBROUTINE OpFM_CreateActForceBladeTowerNodes

SUBROUTINE OpFM_InterpolateForceNodesChord(InitOut_AD, p_OpFM, u_OpFM, ErrStat, ErrMsg) 

  !Interpolates the chord distribution to the force nodes
  ! This routine assumes that the definition of chord along the blade and diameter along the tower are defined at the end points of the blade and tower as well.
  ! This routine approximates the interpolation of chord when BeamDyn is used. BeamDyn uses the curved blade length parametrization for the blade nodes, while AeroDyn uses the radius along the pitch axis. The two definitions will lead to differences in the interpolation of chord when the blade is curved. If the exact chord is needed at the actuator points, please do a formal mesh mapping with blade chord as a scalar on the mesh nodes.
  
  TYPE(AD_InitOutputType),  INTENT(IN   ) :: InitOut_AD ! InitOut  data for the OpenFOAM integration module
  TYPE(OpFM_ParameterType), INTENT(IN   ) :: p_OpFM     ! Input data for the OpenFOAM integration module
  TYPE(OpFM_InputType),     INTENT(INOUT) :: u_OpFM     ! Parameter data for the OpenFOAM integration module
  INTEGER(IntKi)                          :: ErrStat    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                    :: ErrMsg     ! temporary Error message if ErrStat /= ErrID_None

  !Local variables
  INTEGER(IntKI)                         :: i,j,k,node  ! Loop variables
  INTEGER(IntKI)                         :: nNodesBladeProps ! Number of nodes in the blade properties for a given blade
  INTEGER(IntKI)                         :: nNodesTowerProps ! Number of nodes in the tower properties 
  INTEGER(IntKI)                         :: jLower      ! Index of the blade properties node just smaller than the force node
  INTEGER(IntKi)                         :: ErrStat2    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                   :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
  REAL(ReKi)                             :: rInterp     ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
  REAL(ReKi), ALLOCATABLE                :: AD_etaNodes(:)  ! Non-dimensional co-ordinates eta at which the blade and tower chord are defined
  
  ! Set the chord for the hub node to be zero. Ideally, I'd like this to be the hub radius. Will figure this out later.
  Node = 1
  u_OpFM%forceNodesChord(Node) = 0.0_ReKi

  ! The blades first
  do k = 1, p_OpFM%NumBl
     ! Calculate the chord at the force nodes based on interpolation
     nNodesBladeProps = SIZE(InitOut_AD%BladeProps(k)%BlChord)
     allocate(AD_etaNodes(nNodesBladeProps))
     AD_etaNodes = InitOut_AD%BladeProps(k)%BlSpn(:)/InitOut_AD%BladeProps(k)%BlSpn(nNodesBladeProps)
     DO I=1,p_OpFM%NnodesForceBlade
        Node = Node + 1
        jLower=1
        do while ( ( (AD_etaNodes(jLower) - p_OpFM%forceBldEtaNodes(I))*(AD_etaNodes(jLower+1) - p_OpFM%forceBldEtaNodes(I)) .gt. 0 ) .and. (jLower .lt. nNodesBladeProps) )!Determine the closest two nodes at which the blade properties are specified
           jLower = jLower + 1
        end do
        if (jLower .lt. nNodesBladeProps) then
           rInterp =  (p_OpFM%forceBldEtaNodes(I) - AD_etaNodes(jLower))/(AD_etaNodes(jLower+1)-AD_etaNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           u_OpFM%forceNodesChord(Node) = InitOut_AD%BladeProps(k)%BlChord(jLower) + rInterp * (InitOut_AD%BladeProps(k)%BlChord(jLower+1) - InitOut_AD%BladeProps(k)%BlChord(jLower))
        else
           u_OpFM%forceNodesChord(Node) = InitOut_AD%BladeProps(k)%BlChord(nNodesBladeProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn blade properties. Surprisingly this is not an issue with the tower.
        end if
     END DO
     deallocate(AD_etaNodes)     
  end do
     

  ! The tower now
  do k = p_OpFM%NumBl+1,p_OpFM%NMappings
     nNodesTowerProps = SIZE(InitOut_AD%TwrElev)
     allocate(AD_etaNodes(nNodesTowerProps))
     ! Calculate the chord at the force nodes based on interpolation
     AD_etaNodes = InitOut_AD%TwrElev(:)/InitOut_AD%TwrElev(nNodesTowerProps) ! Non-dimensionalize the tower elevation array
     DO I=1,p_OpFM%NnodesForceTower
        Node = Node + 1
        jLower=1
        do while ( ( (AD_etaNodes(jLower) - p_OpFM%forceTwrEtaNodes(I))*(AD_etaNodes(jLower+1) - p_OpFM%forceTwrEtaNodes(I)) .gt. 0) .and. (jLower .lt. nNodesTowerProps) ) !Determine the closest two nodes at which the blade properties are specified
           jLower = jLower + 1
        end do
        if (jLower .lt. nNodesTowerProps) then
           rInterp =  (p_OpFM%forceTwrEtaNodes(I) - AD_etaNodes(jLower))/(AD_etaNodes(jLower+1)-AD_etaNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           u_OpFM%forceNodesChord(Node) = InitOut_AD%TwrDiam(jLower) + rInterp * (InitOut_AD%TwrDiam(jLower+1) - InitOut_AD%TwrDiam(jLower))
        else
           u_OpFM%forceNodesChord(Node) = InitOut_AD%TwrDiam(nNodesTowerProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn tower properties.
        end if
     END DO
  end do

END SUBROUTINE OpFM_InterpolateForceNodesChord

END MODULE OpenFOAM
!**********************************************************************************************************************************
