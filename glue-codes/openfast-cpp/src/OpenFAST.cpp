#include "OpenFAST.H"
#include "hdf5.h"
#include <iostream>
#include <fstream>
#include <cmath>

int fast::OpenFAST::AbortErrLev = ErrID_Fatal; // abort error level; compare with NWTC Library

//Constructor 
fast::fastInputs::fastInputs():
    nTurbinesGlob(0),
    dryRun(false),
    debug(false),
    tStart(-1.0),
    nEveryCheckPoint(-1),
    tMax(0.0),
    dtFAST(0.0),
    scStatus(false),
    scLibFile(""),
    numScInputs(0),
    numScOutputs(0)
{
    //Nothing to do here
}


//Constructor
fast::OpenFAST::OpenFAST():
    nTurbinesGlob(0),
    nTurbinesProc(0),
    scStatus(false),
    simStart(fast::init),
    timeZero(false),
    firstPass_(true)
{
}

inline bool fast::OpenFAST::checkFileExists(const std::string& name) {
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
}

void fast::OpenFAST::init() {
    
    allocateMemory();
    
    if (!dryRun) {
        switch (simStart) {
            
        case fast::trueRestart:
            
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                
                if (turbineData[iTurb].sType == EXTINFLOW) {
                    /* note that this will set nt_global inside the FAST library */
                    FAST_AL_CFD_Restart(&iTurb, turbineData[iTurb].FASTRestartFileName.data(), &AbortErrLev, &dtFAST, &turbineData[iTurb].inflowType, &turbineData[iTurb].numBlades, &turbineData[iTurb].numVelPtsBlade, &turbineData[iTurb].numVelPtsTwr, &ntStart, &extinfw_i_f_FAST[iTurb], &extinfw_o_t_FAST[iTurb], &sc_i_f_FAST[iTurb], &sc_o_t_FAST[iTurb], &ErrStat, ErrMsg);
                    checkError(ErrStat, ErrMsg);
                    nt_global = ntStart;

                    if (nTurbinesProc > 0) velNodeDataFile = openVelocityDataFile(false);                                
                    
                } else if(turbineData[iTurb].sType == EXTLOADS) {
                    FAST_BR_CFD_Restart(&iTurb, turbineData[iTurb].FASTRestartFileName.data(), &AbortErrLev, &dtFAST, &turbineData[iTurb].numBlades, &ntStart, &extld_i_f_FAST[iTurb], &extld_o_t_FAST[iTurb], &sc_i_f_FAST[iTurb], &sc_o_t_FAST[iTurb], &ErrStat, ErrMsg);
                    turbineData[iTurb].inflowType = 0;
                    checkError(ErrStat, ErrMsg);
                    nt_global = ntStart;
                    
                }
                    
                allocateMemory2(iTurb);
                
                readRestartFile(iTurb, nt_global);
                
            }
            
            if(scStatus) {
                sc->readRestartFile(nt_global);
            }
            
            break ;
            
        case fast::init:
            
            // this calls the Init() routines of each module
            
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                if (turbineData[iTurb].sType == EXTINFLOW) {
                
                    FAST_AL_CFD_Init(&iTurb, &tMax, turbineData[iTurb].FASTInputFileName.data(), &turbineData[iTurb].TurbID, &numScOutputs, &numScInputs, &turbineData[iTurb].numForcePtsBlade, &turbineData[iTurb].numForcePtsTwr, turbineData[iTurb].TurbineBasePos.data(), &AbortErrLev, &dtFAST, &turbineData[iTurb].inflowType, &turbineData[iTurb].numBlades, &turbineData[iTurb].numVelPtsBlade, &turbineData[iTurb].numVelPtsTwr, &extinfw_i_f_FAST[iTurb], &extinfw_o_t_FAST[iTurb], &sc_i_f_FAST[iTurb], &sc_o_t_FAST[iTurb], &ErrStat, ErrMsg);
                    checkError(ErrStat, ErrMsg);

                    timeZero = true;

                    turbineData[iTurb].numVelPtsTwr = extinfw_o_t_FAST[iTurb].u_Len - turbineData[iTurb].numBlades*turbineData[iTurb].numVelPtsBlade - 1;
                    if(turbineData[iTurb].numVelPtsTwr == 0) {
                        turbineData[iTurb].numForcePtsTwr = 0;
                        std::cout << "Aerodyn doesn't want to calculate forces on the tower. All actuator points on the tower are turned off for turbine " << turbineMapProcToGlob[iTurb] << "." << std::endl ;
                    }

                    if (nTurbinesProc > 0) velNodeDataFile = openVelocityDataFile(true);                    
                    
                }  else if(turbineData[iTurb].sType == EXTLOADS) {

                    FAST_BR_CFD_Init(&iTurb, &tMax, turbineData[iTurb].FASTInputFileName.data(), &turbineData[iTurb].TurbID, turbineData[iTurb].TurbineBasePos.data(), &AbortErrLev, &dtFAST, &turbineData[iTurb].numBlades, &extld_i_f_FAST[iTurb], &extld_o_t_FAST[iTurb], &sc_i_f_FAST[iTurb], &sc_o_t_FAST[iTurb], &ErrStat, ErrMsg);
                    checkError(ErrStat, ErrMsg);

                    turbineData[iTurb].inflowType = 0;
                    timeZero = true;
                    
                }
                    
                allocateMemory2(iTurb);

                get_refPositions(iTurb);
                
                get_data_from_openfast(fast::nm2);
                get_data_from_openfast(fast::nm1);
                get_data_from_openfast(fast::n);
                get_data_from_openfast(fast::np1);
                
            }           
            
            break ;
            
        case fast::restartDriverInitFAST:
            
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                if (turbineData[iTurb].sType == EXTINFLOW) {
                
                    FAST_AL_CFD_Init(&iTurb, &tMax, turbineData[iTurb].FASTInputFileName.data(), &turbineData[iTurb].TurbID, &numScOutputs, &numScInputs, &turbineData[iTurb].numForcePtsBlade, &turbineData[iTurb].numForcePtsTwr, turbineData[iTurb].TurbineBasePos.data(), &AbortErrLev, &dtFAST, &turbineData[iTurb].inflowType, &turbineData[iTurb].numBlades, &turbineData[iTurb].numVelPtsBlade, &turbineData[iTurb].numVelPtsTwr, &extinfw_i_f_FAST[iTurb], &extinfw_o_t_FAST[iTurb], &sc_i_f_FAST[iTurb], &sc_o_t_FAST[iTurb], &ErrStat, ErrMsg);
                    checkError(ErrStat, ErrMsg);
                
                    timeZero = true;
                    
                    turbineData[iTurb].numVelPtsTwr = extinfw_o_t_FAST[iTurb].u_Len - turbineData[iTurb].numBlades*turbineData[iTurb].numVelPtsBlade - 1;
                    if(turbineData[iTurb].numVelPtsTwr == 0) {
                        turbineData[iTurb].numForcePtsTwr = 0;
                        std::cout << "Aerodyn doesn't want to calculate forces on the tower. All actuator points on the tower are turned off for turbine " << turbineMapProcToGlob[iTurb] << "." << std::endl ;
                    }

                    allocateMemory2(iTurb);
                
                    get_data_from_openfast(fast::nm2);
                    get_data_from_openfast(fast::nm1);
                    get_data_from_openfast(fast::n);
                    get_data_from_openfast(fast::np1);

                    int nTimesteps;
            
                    if (nTurbinesProc > 0) {
                        readVelocityData(ntStart);
                    }
                    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                        applyVelocityData(0, iTurb, extinfw_o_t_FAST[iTurb], velNodeData[iTurb]);
                    }
                    solution0() ;
            
                    for (int iPrestart=0 ; iPrestart < ntStart; iPrestart++) {
                        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                            applyVelocityData(iPrestart, iTurb, extinfw_o_t_FAST[iTurb], velNodeData[iTurb]);
                        }
                        stepNoWrite();
                    }
            
                    if (nTurbinesProc > 0) velNodeDataFile = openVelocityDataFile(false);
                    
                } else if (turbineData[iTurb].sType == EXTLOADS) {

                    throw std::runtime_error("Option restartDriverInitFAST will not work with simulation type EXTLOADS.");
                }

            }
            
            break;
            
        case fast::simStartType_END:
            
            break;
            
        }
        
    }
}

void fast::OpenFAST::solution0() {
    
    if (!dryRun) {
        
        // set wind speeds at initial locations
        // for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        //     setExpLawWindSpeed(iTurb);
        // }

        if(scStatus) {
            
            sc->init(nTurbinesGlob, numScInputs, numScOutputs);
            
            sc->calcOutputs(scOutputsGlob);
            fillScOutputsLoc();
        }
        
        // Unfortunately setVelocity only sets the velocity at 'n+1'. Need to copy 'n+1' to 'n'
        init_velForceNodeData() ;
        send_data_to_openfast(fast::np1);
        
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            
            FAST_CFD_Solution0(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
            
            FAST_CFD_InitIOarrays_SS(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);       
        }
        
        get_data_from_openfast(fast::n);
        get_data_from_openfast(fast::nm1);
        get_data_from_openfast(fast::nm2);        
        
        timeZero = false;

        if (scStatus) {
            fillScInputsGlob(); // Update inputs to super controller
        }
    }
    
}

void fast::OpenFAST::get_refPositions(int iTurb) {
    
    if(turbineData[iTurb].sType == EXTLOADS) {
        
        int nBlades = turbineData[iTurb].numBlades;
        int iRunTot = 0;
        for (int i=0; i < nBlades; i++) {
            int nPtsBlade = turbineData[iTurb].nBRfsiPtsBlade[i];
            for (int j=0; j < nPtsBlade; j++) {
                for (int k=0; k < 3; k++) {
                    brFSIData[iTurb][fast::np1].bld_ref_pos[iRunTot*6+k] = extld_i_f_FAST[iTurb].bldRefPos[iRunTot*6+k] + turbineData[iTurb].TurbineBasePos[k];
                    brFSIData[iTurb][fast::np1].bld_ref_pos[iRunTot*6+k+3] = extld_i_f_FAST[iTurb].bldRefPos[iRunTot*6+k+3];                    
                }
                iRunTot++;
            }
        }
        
        int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
        for (int i=0; i < nPtsTwr; i++) {
            for (int j = 0; j < 3; j++) {
                brFSIData[iTurb][fast::np1].twr_ref_pos[i*6+j] = extld_i_f_FAST[iTurb].twrRefPos[i*6+j] + turbineData[iTurb].TurbineBasePos[j];
                brFSIData[iTurb][fast::np1].twr_ref_pos[i*6+j+3] = extld_i_f_FAST[iTurb].twrRefPos[i*6+j+3];
            }
        }
        
    }
    
}

void fast::OpenFAST::init_velForceNodeData() {

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        int nvelpts = get_numVelPtsLoc(iTurb);
        int nfpts = get_numForcePtsLoc(iTurb);
        for (int i=0; i<nvelpts; i++) {
            for (int j=0 ; j < 3; j++) {
                velForceNodeData[iTurb][fast::n].x_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].x_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm1].x_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].x_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm2].x_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].x_vel[i*3+j];
                velForceNodeData[iTurb][fast::n].xdot_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm1].xdot_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm2].xdot_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_vel[i*3+j];
                velForceNodeData[iTurb][fast::n].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm1].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm2].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_vel[i*3+j];
            }
        }
        for (int i=0; i<nfpts; i++) {
            for (int j=0 ; j < 3; j++) {
                velForceNodeData[iTurb][fast::n].x_force[i*3+j] = velForceNodeData[iTurb][fast::np1].x_force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].x_force[i*3+j] = velForceNodeData[iTurb][fast::np1].x_force[i*3+j];
                velForceNodeData[iTurb][fast::nm2].x_force[i*3+j] = velForceNodeData[iTurb][fast::np1].x_force[i*3+j];
                velForceNodeData[iTurb][fast::n].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_force[i*3+j];
                velForceNodeData[iTurb][fast::nm2].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_force[i*3+j];
                velForceNodeData[iTurb][fast::n].vel_force[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].vel_force[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_force[i*3+j];
                velForceNodeData[iTurb][fast::nm2].vel_force[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_force[i*3+j];
                velForceNodeData[iTurb][fast::n].force[i*3+j] = velForceNodeData[iTurb][fast::np1].force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].force[i*3+j] = velForceNodeData[iTurb][fast::np1].force[i*3+j];
                velForceNodeData[iTurb][fast::nm2].force[i*3+j] = velForceNodeData[iTurb][fast::np1].force[i*3+j];
            }
            for (int j=0;j<9;j++) {
                velForceNodeData[iTurb][fast::n].orient_force[i*9+j] = velForceNodeData[iTurb][fast::np1].orient_force[i*9+j];
                velForceNodeData[iTurb][fast::nm1].orient_force[i*9+j] = velForceNodeData[iTurb][fast::np1].orient_force[i*9+j];
                velForceNodeData[iTurb][fast::nm2].orient_force[i*9+j] = velForceNodeData[iTurb][fast::np1].orient_force[i*9+j];                 
            }
        }
    }
    
}

void fast::OpenFAST::predict_states() {
    
    if (firstPass_) {
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            int nfpts = get_numForcePtsLoc(iTurb);
            for (int i=0; i<nvelpts; i++) {
                for (int j=0 ; j < 3; j++) {
                    velForceNodeData[iTurb][fast::np1].x_vel[i*3+j] = velForceNodeData[iTurb][fast::nm2].x_vel[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].x_vel[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].x_vel[i*3+j];
                    velForceNodeData[iTurb][fast::np1].xdot_vel[i*3+j] = velForceNodeData[iTurb][fast::nm2].xdot_vel[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].xdot_vel[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].xdot_vel[i*3+j];
                    velForceNodeData[iTurb][fast::np1].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::nm2].vel_vel[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].vel_vel[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].vel_vel[i*3+j];
                }
            }
            velForceNodeData[iTurb][fast::np1].x_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].xdot_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].vel_vel_resid = 0.0;        
            for (int i=0; i<nfpts; i++) {
                for (int j=0 ; j < 3; j++) {
                    velForceNodeData[iTurb][fast::np1].x_force[i*3+j] = velForceNodeData[iTurb][fast::nm2].x_force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].x_force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].x_force[i*3+j];
                    velForceNodeData[iTurb][fast::np1].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::nm2].xdot_force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].xdot_force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].xdot_force[i*3+j];
                    velForceNodeData[iTurb][fast::np1].vel_force[i*3+j] = velForceNodeData[iTurb][fast::nm2].vel_force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].vel_force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].vel_force[i*3+j];
                    velForceNodeData[iTurb][fast::np1].force[i*3+j] = velForceNodeData[iTurb][fast::nm2].force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].force[i*3+j];
                }
                for (int j=0;j<9;j++)
                    velForceNodeData[iTurb][fast::np1].orient_force[i*9+j] = velForceNodeData[iTurb][fast::nm2].orient_force[i*3+j] + 3.0*velForceNodeData[iTurb][fast::n].orient_force[i*3+j] - 3.0*velForceNodeData[iTurb][fast::nm1].orient_force[i*3+j];
            }
            velForceNodeData[iTurb][fast::np1].x_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].xdot_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].orient_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].vel_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].force_resid = 0.0;
        }
    }
}

void fast::OpenFAST::prework() {
    
    if (nSubsteps_ > 1) {
        
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CFD_Store_SS(&iTurb, &nt_global, &ErrStat, ErrMsg) ;
            checkError(ErrStat, ErrMsg);
        }
        
    } else {
        
        if(scStatus) {
            sc->calcOutputs(scOutputsGlob);
            fillScOutputsLoc();
        }
        
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CFD_Prework(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
        }
    }
    
}

void fast::OpenFAST::update_states_driver_time_step() {
    
    if (firstPass_)
        prework();
    
    if (nSubsteps_ > 1) {

        if (!firstPass_) {
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                FAST_CFD_Reset_SS(&iTurb, &nSubsteps_, &ErrStat, ErrMsg);
                checkError(ErrStat, ErrMsg);
            }
            nt_global = nt_global - nSubsteps_;
        }
        
        for (int iSubstep=1; iSubstep < nSubsteps_+1; iSubstep++) {
            double ss_time = double(iSubstep)/double(nSubsteps_);
            step(ss_time);
        }
        
        get_data_from_openfast(fast::np1);

        if ( isDebug() ) {
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                std::ofstream fastcpp_velocity_file;
                fastcpp_velocity_file.open("fastcpp_residual." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv", std::ios_base::app) ;
                fastcpp_velocity_file << "Time step " << nt_global << " Velocity residual at the force nodes = " << velForceNodeData[iTurb][fast::np1].vel_force_resid << std::endl ;
                fastcpp_velocity_file << "          " << nt_global << " Position residual at the force nodes = " << velForceNodeData[iTurb][fast::np1].x_force_resid << std::endl ;
                fastcpp_velocity_file << "          " << nt_global << " Force residual at the force nodes = " << velForceNodeData[iTurb][fast::np1].force_resid << std::endl ;
                fastcpp_velocity_file.close() ;
            }
        }
        
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            velForceNodeData[iTurb][fast::np1].x_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].xdot_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].vel_vel_resid = 0.0;        
            velForceNodeData[iTurb][fast::np1].x_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].xdot_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].orient_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].vel_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].force_resid = 0.0;
        }        
    } else {
        
        send_data_to_openfast(fast::np1);
        
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CFD_UpdateStates(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);

            // Compute the force from the nacelle only if the drag coefficient is
            //   greater than zero
            if (get_nacelleCdLoc(iTurb) > 0.) {

                calc_nacelle_force (
                             
                                    extinfw_o_t_FAST[iTurb].u[0], 
                                    extinfw_o_t_FAST[iTurb].v[0], 
                                    extinfw_o_t_FAST[iTurb].w[0], 
                                    get_nacelleCdLoc(iTurb), 
                                    get_nacelleAreaLoc(iTurb), 
                                    get_airDensityLoc(iTurb), 
                                    extinfw_i_f_FAST[iTurb].fx[0], 
                                    extinfw_i_f_FAST[iTurb].fy[0], 
                                    extinfw_i_f_FAST[iTurb].fz[0]

                                    );

            }
            
        }
        
        get_data_from_openfast(fast::np1);

        if ( isDebug() ) {
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                std::ofstream fastcpp_velocity_file;
                fastcpp_velocity_file.open("fastcpp_residual." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv", std::ios_base::app) ;
                fastcpp_velocity_file << "Time step " << nt_global << " Velocity residual at the force nodes = " << velForceNodeData[iTurb][fast::np1].vel_force_resid << std::endl ;
                fastcpp_velocity_file << "          " << nt_global << " Position residual at the force nodes = " << velForceNodeData[iTurb][fast::np1].x_force_resid << std::endl ;
                fastcpp_velocity_file << "          " << nt_global << " Force residual at the force nodes = " << velForceNodeData[iTurb][fast::np1].force_resid << std::endl ;
                fastcpp_velocity_file.close() ;
            }
        }

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            velForceNodeData[iTurb][fast::np1].x_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].xdot_vel_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].vel_vel_resid = 0.0;        
            velForceNodeData[iTurb][fast::np1].x_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].xdot_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].orient_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].vel_force_resid = 0.0;
            velForceNodeData[iTurb][fast::np1].force_resid = 0.0;
        }        

    }
    firstPass_ = false;    
}

void fast::OpenFAST::advance_to_next_driver_time_step() {
    
    if (nSubsteps_ > 1) {
        //Nothing to do here
        
    } else {
        
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            
            if (turbineData[iTurb].inflowType == 2)
                writeVelocityData(velNodeDataFile, iTurb, nt_global, extinfw_i_f_FAST[iTurb], extinfw_o_t_FAST[iTurb]);
            
            if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {
                
                std::ofstream fastcpp_velocity_file;
                fastcpp_velocity_file.open("fastcpp_velocity." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
                fastcpp_velocity_file << "# x, y, z, Vx, Vy, Vz" << std::endl ;
                for (int iNode=0; iNode < get_numVelPtsLoc(iTurb); iNode++) {
                    fastcpp_velocity_file << extinfw_i_f_FAST[iTurb].pxVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzVel[iNode] << ", " << extinfw_o_t_FAST[iTurb].u[iNode] << ", " << extinfw_o_t_FAST[iTurb].v[iNode] << ", " << extinfw_o_t_FAST[iTurb].w[iNode] << " " << std::endl ;           
                }
                fastcpp_velocity_file.close() ;
                
            }
            
            FAST_CFD_AdvanceToNextTimeStep(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
        
            if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {
                std::ofstream actuatorForcesFile;
                actuatorForcesFile.open("actuator_forces." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
                actuatorForcesFile << "# x, y, z, fx, fy, fz" << std::endl ;
                for (int iNode=0; iNode < get_numForcePtsLoc(iTurb); iNode++) {
                    actuatorForcesFile << extinfw_i_f_FAST[iTurb].pxForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].fx[iNode] << ", " << extinfw_i_f_FAST[iTurb].fy[iNode] << ", " << extinfw_i_f_FAST[iTurb].fz[iNode] << " " << std::endl ;           
                }
                actuatorForcesFile.close() ;
            }
            
        }
        
        if(scStatus) {
            sc->updateStates(scInputsGlob); // Go from 'n' to 'n+1' based on input at previous time step
            fillScInputsGlob(); // Update inputs to super controller for 'n+1'
        }
        
        nt_global = nt_global + 1;
    }

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        FAST_CFD_WriteOutput(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
    }
        
    if ( (((nt_global - ntStart) % nEveryCheckPoint) == 0 )  && (nt_global != ntStart) ) {
        //sprintf(CheckpointFileRoot, "../../CertTest/Test18.%d", nt_global);
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            turbineData[iTurb].FASTRestartFileName = " "; // if blank, it will use FAST convention <RootName>.nt_global
            FAST_CreateCheckpoint(&iTurb, turbineData[iTurb].FASTRestartFileName.data(), &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
            writeRestartFile(iTurb, nt_global);
        }
        if(scStatus) {
            if (fastMPIRank == 0) {
                sc->writeRestartFile(nt_global);
            }
        }
    }
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        int nvelpts = get_numVelPtsLoc(iTurb);
        int nfpts = get_numForcePtsLoc(iTurb);
        for (int i=0; i<nvelpts; i++) {
            for (int j=0 ; j < 3; j++) {
                velForceNodeData[iTurb][fast::nm2].x_vel[i*3+j] = velForceNodeData[iTurb][fast::nm1].x_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm1].x_vel[i*3+j] = velForceNodeData[iTurb][fast::n].x_vel[i*3+j];
                velForceNodeData[iTurb][fast::n].x_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].x_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm2].xdot_vel[i*3+j] = velForceNodeData[iTurb][fast::nm1].xdot_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm1].xdot_vel[i*3+j] = velForceNodeData[iTurb][fast::n].xdot_vel[i*3+j];
                velForceNodeData[iTurb][fast::n].xdot_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm2].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::nm1].vel_vel[i*3+j];
                velForceNodeData[iTurb][fast::nm1].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::n].vel_vel[i*3+j];
                velForceNodeData[iTurb][fast::n].vel_vel[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_vel[i*3+j];
            }
        }
        
        for (int i=0; i<nfpts; i++) {
            for (int j=0 ; j < 3; j++) {
                velForceNodeData[iTurb][fast::nm2].x_force[i*3+j] = velForceNodeData[iTurb][fast::nm1].x_force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].x_force[i*3+j] = velForceNodeData[iTurb][fast::n].x_force[i*3+j];
                velForceNodeData[iTurb][fast::n].x_force[i*3+j] = velForceNodeData[iTurb][fast::np1].x_force[i*3+j];
                velForceNodeData[iTurb][fast::nm2].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::nm1].xdot_force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::n].xdot_force[i*3+j];
                velForceNodeData[iTurb][fast::n].xdot_force[i*3+j] = velForceNodeData[iTurb][fast::np1].xdot_force[i*3+j];
                velForceNodeData[iTurb][fast::nm2].vel_force[i*3+j] = velForceNodeData[iTurb][fast::nm1].vel_force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].vel_force[i*3+j] = velForceNodeData[iTurb][fast::n].vel_force[i*3+j];
                velForceNodeData[iTurb][fast::n].vel_force[i*3+j] = velForceNodeData[iTurb][fast::np1].vel_force[i*3+j];
                velForceNodeData[iTurb][fast::nm2].force[i*3+j] = velForceNodeData[iTurb][fast::nm1].force[i*3+j];
                velForceNodeData[iTurb][fast::nm1].force[i*3+j] = velForceNodeData[iTurb][fast::n].force[i*3+j];
                velForceNodeData[iTurb][fast::n].force[i*3+j] = velForceNodeData[iTurb][fast::np1].force[i*3+j];
            }
            for (int j=0;j<9;j++) {
                velForceNodeData[iTurb][fast::nm2].orient_force[i*9+j] = velForceNodeData[iTurb][fast::nm1].orient_force[i*9+j];
                velForceNodeData[iTurb][fast::nm1].orient_force[i*9+j] = velForceNodeData[iTurb][fast::n].orient_force[i*9+j];
                velForceNodeData[iTurb][fast::n].orient_force[i*9+j] = velForceNodeData[iTurb][fast::np1].orient_force[i*9+j];
            }
        }
    }

    firstPass_ = true ; // Set firstPass_ to true for the next time step
    
}

/* A version of step allowing for sub-timesteps when the driver program has a larger time step than OpenFAST */
void fast::OpenFAST::step(double ss_time) {
    
    /* ******************************
       set inputs from this code and call FAST:
       ********************************* */
    
    if(scStatus) {
        sc->calcOutputs(scOutputsGlob);
        fillScOutputsLoc();
    }
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        
        //  set wind speeds at original locations 
        // setExpLawWindSpeed(iTurb);
        
        // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
        // (note CFD could do subcycling around this step)
        
        if (turbineData[iTurb].inflowType == 2)
            writeVelocityData(velNodeDataFile, iTurb, nt_global, extinfw_i_f_FAST[iTurb], extinfw_o_t_FAST[iTurb]);
        
        if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {
            
            std::ofstream fastcpp_velocity_file;
            fastcpp_velocity_file.open("fastcpp_velocity." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
            fastcpp_velocity_file << "# x, y, z, Vx, Vy, Vz" << std::endl ;
            for (int iNode=0; iNode < get_numVelPtsLoc(iTurb); iNode++) {
                fastcpp_velocity_file << extinfw_i_f_FAST[iTurb].pxVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzVel[iNode] << ", " << extinfw_o_t_FAST[iTurb].u[iNode] << ", " << extinfw_o_t_FAST[iTurb].v[iNode] << ", " << extinfw_o_t_FAST[iTurb].w[iNode] << " " << std::endl ;           
            }
            fastcpp_velocity_file.close() ;
        }
        
        FAST_CFD_Prework(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        send_data_to_openfast(ss_time);
        FAST_CFD_UpdateStates(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        FAST_CFD_AdvanceToNextTimeStep(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        
        if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {
            std::ofstream actuatorForcesFile;
            actuatorForcesFile.open("actuator_forces." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
            actuatorForcesFile << "# x, y, z, fx, fy, fz" << std::endl ;
            for (int iNode=0; iNode < get_numForcePtsLoc(iTurb); iNode++) {
                actuatorForcesFile << extinfw_i_f_FAST[iTurb].pxForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].fx[iNode] << ", " << extinfw_i_f_FAST[iTurb].fy[iNode] << ", " << extinfw_i_f_FAST[iTurb].fz[iNode] << " " << std::endl ;           
            }
            actuatorForcesFile.close() ;
        }
        
    }
    
    if(scStatus) {
        sc->updateStates(scInputsGlob); // Go from 'n' to 'n+1' based on input at previous time step
        fillScInputsGlob(); // Update inputs to super controller for 'n+1'
    }
    
    nt_global = nt_global + 1;
    
}

void fast::OpenFAST::step() {
    
    /* ******************************
       set inputs from this code and call FAST:
       ********************************* */
    
    if(scStatus) {
        sc->calcOutputs(scOutputsGlob);
        fillScOutputsLoc();
    }
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        
        //  set wind speeds at original locations 
        // setExpLawWindSpeed(iTurb);
        
        // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
        // (note CFD could do subcycling around this step)
        
        if (turbineData[iTurb].inflowType == 2)
            writeVelocityData(velNodeDataFile, iTurb, nt_global, extinfw_i_f_FAST[iTurb], extinfw_o_t_FAST[iTurb]);
        
        if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {
            
            std::ofstream fastcpp_velocity_file;
            fastcpp_velocity_file.open("fastcpp_velocity." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
            fastcpp_velocity_file << "# x, y, z, Vx, Vy, Vz" << std::endl ;
            for (int iNode=0; iNode < get_numVelPtsLoc(iTurb); iNode++) {
                fastcpp_velocity_file << extinfw_i_f_FAST[iTurb].pxVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyVel[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzVel[iNode] << ", " << extinfw_o_t_FAST[iTurb].u[iNode] << ", " << extinfw_o_t_FAST[iTurb].v[iNode] << ", " << extinfw_o_t_FAST[iTurb].w[iNode] << " " << std::endl ;           
            }
            fastcpp_velocity_file.close() ;
        }
        
        FAST_CFD_Prework(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        send_data_to_openfast(fast::np1);
        FAST_CFD_UpdateStates(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        get_data_from_openfast(fast::np1);
        FAST_CFD_AdvanceToNextTimeStep(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        FAST_CFD_WriteOutput(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        

        // Compute the force from the nacelle only if the drag coefficient is
        //   greater than zero
        if (get_nacelleCdLoc(iTurb) > 0.) {

            calc_nacelle_force (
                             
                extinfw_o_t_FAST[iTurb].u[0], 
                extinfw_o_t_FAST[iTurb].v[0], 
                extinfw_o_t_FAST[iTurb].w[0], 
                get_nacelleCdLoc(iTurb), 
                get_nacelleAreaLoc(iTurb), 
                get_airDensityLoc(iTurb), 
                extinfw_i_f_FAST[iTurb].fx[0], 
                extinfw_i_f_FAST[iTurb].fy[0], 
                extinfw_i_f_FAST[iTurb].fz[0]

                );

        }
        
        if ( isDebug() && (turbineData[iTurb].inflowType == 2) ) {
            std::ofstream actuatorForcesFile;
            actuatorForcesFile.open("actuator_forces." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
            actuatorForcesFile << "# x, y, z, fx, fy, fz" << std::endl ;
            for (int iNode=0; iNode < get_numForcePtsLoc(iTurb); iNode++) {
                actuatorForcesFile << extinfw_i_f_FAST[iTurb].pxForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pyForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].pzForce[iNode] << ", " << extinfw_i_f_FAST[iTurb].fx[iNode] << ", " << extinfw_i_f_FAST[iTurb].fy[iNode] << ", " << extinfw_i_f_FAST[iTurb].fz[iNode] << " " << std::endl ;           
            }
            actuatorForcesFile.close() ;
        }
        
    }
    
    if(scStatus) {
        sc->updateStates(scInputsGlob); // Go from 'n' to 'n+1' based on input at previous time step
        fillScInputsGlob(); // Update inputs to super controller for 'n+1'
    }
    
    nt_global = nt_global + 1;
    
    if ( (((nt_global - ntStart) % nEveryCheckPoint) == 0 )  && (nt_global != ntStart) ) {
        //sprintf(FASTRestartFileName, "../../CertTest/Test18.%d", nt_global);
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            turbineData[iTurb].FASTRestartFileName = " "; // if blank, it will use FAST convention <RootName>.nt_global
            FAST_CreateCheckpoint(&iTurb, turbineData[iTurb].FASTRestartFileName.data(), &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
            writeRestartFile(iTurb, nt_global);
        }
        if(scStatus) {
            if (fastMPIRank == 0) {
                sc->writeRestartFile(nt_global);
            }
        }
    }
    
}

void fast::OpenFAST::stepNoWrite() {
    
    /* ******************************
       set inputs from this code and call FAST:
       ********************************* */
    
    if(scStatus) {
        sc->calcOutputs(scOutputsGlob);
        fillScOutputsLoc();
    }
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        
        //  set wind speeds at original locations 
        // setExpLawWindSpeed(iTurb);
        
        // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
        // (note CFD could do subcycling around this step)
        FAST_CFD_Prework(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        send_data_to_openfast(fast::np1);
        FAST_CFD_UpdateStates(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        get_data_from_openfast(fast::np1);
        FAST_CFD_AdvanceToNextTimeStep(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);
        FAST_CFD_WriteOutput(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);

        // Compute the force from the nacelle only if the drag coefficient is
        //   greater than zero
        if (get_nacelleCdLoc(iTurb) > 0.) {

            calc_nacelle_force (
                             
                extinfw_o_t_FAST[iTurb].u[0], 
                extinfw_o_t_FAST[iTurb].v[0], 
                extinfw_o_t_FAST[iTurb].w[0], 
                get_nacelleCdLoc(iTurb), 
                get_nacelleAreaLoc(iTurb), 
                get_airDensityLoc(iTurb), 
                extinfw_i_f_FAST[iTurb].fx[0], 
                extinfw_i_f_FAST[iTurb].fy[0], 
                extinfw_i_f_FAST[iTurb].fz[0]

                );

        }
        
    }
    
    if(scStatus) {
        sc->updateStates(scInputsGlob); // Go from 'n' to 'n+1' based on input at previous time step
        fillScInputsGlob(); // Update inputs to super controller for 'n+1'
    }
    
    nt_global = nt_global + 1;
    
}

void fast::OpenFAST::calc_nacelle_force(
        const float & u, 
        const float & v, 
        const float & w, 
        const float & cd, 
        const float & area, 
        const float & rho,
        float & fx, 
        float & fy, 
        float & fz) {
            // Calculate the force on the nacelle (fx,fy,fz) given the 
            //   velocity sampled at the nacelle point (u,v,w), 
            //   drag coefficient 'cd' and nacelle area 'area'
    
            // The velocity magnitude
            float Vmag = std::sqrt(u * u + v * v + w * w);
    
            // Velocity correction based on Martinez-Tossas PhD Thesis 2017
            // The correction samples the velocity at the center of the
            // Gaussian kernel and scales it to obtain the inflow velocity 
            float epsilon_d = std::sqrt(2.0 / M_PI * cd * area);
            float correction = 1. / (1.0 - cd * area /
                                        (4.0 * M_PI * epsilon_d * epsilon_d));
    
            // Compute the force for each velocity component
            fx = rho * 1./2. * cd * area * Vmag * u * correction * correction;
            fy = rho * 1./2. * cd * area * Vmag * v * correction * correction;
            fz = rho * 1./2. * cd * area * Vmag * w * correction * correction;
        }

void fast::OpenFAST::setInputs(const fast::fastInputs & fi ) {
    
    
    mpiComm = fi.comm;
    
    MPI_Comm_rank(mpiComm, &worldMPIRank);
    MPI_Comm_group(mpiComm, &worldMPIGroup);
    
    nTurbinesGlob = fi.nTurbinesGlob;
    
    if (nTurbinesGlob > 0) {

        dryRun = fi.dryRun;
        
        debug = fi.debug;
        
        tStart = fi.tStart;
        simStart = fi.simStart;
        nEveryCheckPoint = fi.nEveryCheckPoint;
        tMax = fi.tMax;
        loadSuperController(fi);
        dtFAST = fi.dtFAST;
        nSubsteps_ = fi.nSubsteps;
        
        ntStart = int(tStart/dtFAST);
        
        if (simStart == fast::restartDriverInitFAST) {
            nt_global = 0;
        } else {
            nt_global = ntStart;
        }
        
        globTurbineData.resize(nTurbinesGlob);
        globTurbineData = fi.globTurbineData;
        
    } else {
        throw std::runtime_error("Number of turbines < 0 ");
    }
    
}

void fast::OpenFAST::get_turbineParams(int iTurbGlob, turbineDataType & turbData) {

    //TODO: Figure out a better copy operator for the turbineDataType struct
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    turbData.TurbID = turbineData[iTurbLoc].TurbID;
    turbData.FASTInputFileName = turbineData[iTurbLoc].FASTInputFileName;
    turbData.FASTRestartFileName = turbineData[iTurbLoc].FASTRestartFileName;
    if(turbineData[iTurbLoc].TurbineBasePos.size() > 0) {
        turbData.TurbineBasePos.resize(turbineData[iTurbLoc].TurbineBasePos.size());
        for(int i=0; i < turbineData[iTurbLoc].TurbineBasePos.size(); i++)
            turbData.TurbineBasePos[i] = turbineData[iTurbLoc].TurbineBasePos[i];
    }
    if(turbineData[iTurbLoc].TurbineHubPos.size() > 0) {
        turbData.TurbineHubPos.resize(turbineData[iTurbLoc].TurbineHubPos.size());
        for(int i=0; i < turbineData[iTurbLoc].TurbineHubPos.size(); i++)
            turbData.TurbineHubPos[i] = turbineData[iTurbLoc].TurbineHubPos[i];
    }
    turbData.sType = turbineData[iTurbLoc].sType;
    turbData.numBlades = turbineData[iTurbLoc].numBlades;
    turbData.numVelPtsBlade = turbineData[iTurbLoc].numVelPtsBlade;
    turbData.numVelPtsTwr = turbineData[iTurbLoc].numVelPtsTwr;
    turbData.numVelPts = turbineData[iTurbLoc].numVelPts;
    turbData.numForcePtsBlade = turbineData[iTurbLoc].numForcePtsBlade;
    turbData.numForcePtsTwr = turbineData[iTurbLoc].numForcePtsTwr;
    turbData.numForcePts = turbineData[iTurbLoc].numForcePts;
    turbData.inflowType = turbineData[iTurbLoc].inflowType;
    turbData.nacelle_cd = turbineData[iTurbLoc].nacelle_cd;
    turbData.nacelle_area = turbineData[iTurbLoc].nacelle_area;
    turbData.air_density = turbineData[iTurbLoc].air_density;
    turbData.nBRfsiPtsBlade.resize(turbData.numBlades);
    for (int i=0; i < turbData.numBlades; i++)
        turbData.nBRfsiPtsBlade[i] = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
    turbData.nBRfsiPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    
}


void fast::OpenFAST::checkError(const int ErrStat, const char * ErrMsg){
    
    if (ErrStat != ErrID_None){
        
        if (ErrStat >= AbortErrLev){
            throw std::runtime_error(std::string(ErrMsg));
        }
        
    }
    
}

void fast::OpenFAST::setExpLawWindSpeed(){
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        // routine sets the u-v-w wind speeds used in FAST 
        int nVelPts = get_numVelPts(iTurb);
        int iTurbGlob = turbineMapProcToGlob[iTurb];
        for (int j = 0; j < nVelPts; j++){
            std::vector<double> coords(3,0.0);
            std::vector<double> tmpVel(3,0.0);
            getVelNodeCoordinates(coords, j, iTurbGlob, fast::np1);
            tmpVel[0] = (float) 10.0*pow((coords[2] / 90.0), 0.2); // 0.2 power law wind profile using reference 10 m/s at 90 meters
            setVelocity(tmpVel, j, iTurbGlob);
        }
    }
}

void fast::OpenFAST::getApproxHubPos(std::vector<double> & currentCoords, int iTurbGlob) {
    
    // Get hub position of Turbine 'iTurbGlob'
    currentCoords[0] = globTurbineData[iTurbGlob].TurbineHubPos[0];
    currentCoords[1] = globTurbineData[iTurbGlob].TurbineHubPos[1];
    currentCoords[2] = globTurbineData[iTurbGlob].TurbineHubPos[2];
    
}

void fast::OpenFAST::getHubPos(std::vector<double> & currentCoords, int iTurbGlob, fast::timeStep t) {
    
    // Get hub position of Turbine 'iTurbGlob'
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    currentCoords[0] = velForceNodeData[iTurbLoc][t].x_force[0] + turbineData[iTurbLoc].TurbineBasePos[0] ;
    currentCoords[1] = velForceNodeData[iTurbLoc][t].x_force[1] + turbineData[iTurbLoc].TurbineBasePos[1] ;
    currentCoords[2] = velForceNodeData[iTurbLoc][t].x_force[2] + turbineData[iTurbLoc].TurbineBasePos[2] ;
    
}

void fast::OpenFAST::getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob, fast::timeStep t) {
    
    // Get hub shaft direction of current turbine - pointing downwind
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    hubShftVec[0] = velForceNodeData[iTurbLoc][t].orient_force[0] ;
    hubShftVec[1] = velForceNodeData[iTurbLoc][t].orient_force[3] ;
    hubShftVec[2] = velForceNodeData[iTurbLoc][t].orient_force[6] ;
    
}


void fast::OpenFAST::getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t) {
    
    // Set coordinates at current node of current turbine 
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
    currentCoords[0] = velForceNodeData[iTurbLoc][t].x_vel[iNode*3+0] + turbineData[iTurbLoc].TurbineBasePos[0] ;
    currentCoords[1] = velForceNodeData[iTurbLoc][t].x_vel[iNode*3+1] + turbineData[iTurbLoc].TurbineBasePos[1] ;
    currentCoords[2] = velForceNodeData[iTurbLoc][t].x_vel[iNode*3+2] + turbineData[iTurbLoc].TurbineBasePos[2] ;
    
}

void fast::OpenFAST::getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t) {
    
    // Set coordinates at current node of current turbine 
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    currentCoords[0] = velForceNodeData[iTurbLoc][t].x_force[iNode*3+0] + turbineData[iTurbLoc].TurbineBasePos[0] ;
    currentCoords[1] = velForceNodeData[iTurbLoc][t].x_force[iNode*3+1] + turbineData[iTurbLoc].TurbineBasePos[1] ;
    currentCoords[2] = velForceNodeData[iTurbLoc][t].x_force[iNode*3+2] + turbineData[iTurbLoc].TurbineBasePos[2] ;
    
}

void fast::OpenFAST::getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob, fast::timeStep t) {
    
    // Set orientation at current node of current turbine 
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    for(int i=0;i<9;i++) {
        currentOrientation[i] = velForceNodeData[iTurbLoc][t].orient_force[iNode*9+i] ;
    }
    
}

void fast::OpenFAST::getForce(std::vector<double> & currentForce, int iNode, int iTurbGlob, fast::timeStep t) {
    
    // Set forces at current node of current turbine 
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    currentForce[0] = -velForceNodeData[iTurbLoc][t].force[iNode*3+0] ;
    currentForce[1] = -velForceNodeData[iTurbLoc][t].force[iNode*3+1] ;
    currentForce[2] = -velForceNodeData[iTurbLoc][t].force[iNode*3+2] ;
    
}

double fast::OpenFAST::getChord(int iNode, int iTurbGlob) {
    
    // Return blade chord/tower diameter at current node of current turbine 
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    return extinfw_i_f_FAST[iTurbLoc].forceNodesChord[iNode] ;
    
}

void fast::OpenFAST::setVelocity(std::vector<double> & currentVelocity, int iNode, int iTurbGlob) {
    
    // Set velocity at current node of current turbine - 
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
    for(int k=0; k < 3; k++) {
        velForceNodeData[iTurbLoc][fast::np1].vel_vel_resid += (velForceNodeData[iTurbLoc][fast::np1].vel_vel[iNode*3+k] - currentVelocity[k])*(velForceNodeData[iTurbLoc][fast::np1].vel_vel[iNode*3+k] - currentVelocity[k]);
        velForceNodeData[iTurbLoc][fast::np1].vel_vel[iNode*3+k] = currentVelocity[k];
    }
    
    // Put this in send_data_to_openfast
    extinfw_o_t_FAST[iTurbLoc].u[iNode] = currentVelocity[0];
    extinfw_o_t_FAST[iTurbLoc].v[iNode] = currentVelocity[1];
    extinfw_o_t_FAST[iTurbLoc].w[iNode] = currentVelocity[2];
}

void fast::OpenFAST::setVelocityForceNode(std::vector<double> & currentVelocity, int iNode, int iTurbGlob) {
    
    // Set velocity at current node of current turbine - 
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    for(int k=0; k < 3; k++) {
        velForceNodeData[iTurbLoc][fast::np1].vel_force_resid += (velForceNodeData[iTurbLoc][fast::np1].vel_force[iNode*3+k] - currentVelocity[k])*(velForceNodeData[iTurbLoc][fast::np1].vel_force[iNode*3+k] - currentVelocity[k]);
        velForceNodeData[iTurbLoc][fast::np1].vel_force[iNode*3+k] = currentVelocity[k];
    }
}

void fast::OpenFAST::interpolateVel_ForceToVelNodes() {
    
    // Interpolates the velocity from the force nodes to the velocity nodes
    
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        // Hub location
        
        for (int k=0; k < 3; k++) {
            double tmp = velForceNodeData[iTurb][fast::np1].vel_force[k];
            velForceNodeData[iTurb][fast::np1].vel_vel_resid += (velForceNodeData[iTurb][fast::np1].vel_vel[k] - tmp)*(velForceNodeData[iTurb][fast::np1].vel_vel[k] - tmp);
            velForceNodeData[iTurb][fast::np1].vel_vel[k] = tmp;
        }
        
        if ( isDebug() ) {
            std::ofstream actuatorVelFile;
            actuatorVelFile.open("actuator_velocity." + std::to_string(turbineMapProcToGlob[iTurb]) + ".csv") ;
            actuatorVelFile << "# x, y, z, Vx, Vy, Vz" << std::endl ;
            for (int iNode=0; iNode < get_numForcePtsLoc(iTurb); iNode++) {
                actuatorVelFile << velForceNodeData[iTurb][fast::np1].force[iNode*3+0] << ", " << velForceNodeData[iTurb][fast::np1].force[iNode*3+0] << ", " << velForceNodeData[iTurb][fast::np1].force[iNode*3+0] << ", " << velForceNodeData[iTurb][fast::np1].vel_force[iNode*3+0] << ", " << velForceNodeData[iTurb][fast::np1].vel_force[iNode*3+0] << ", " << velForceNodeData[iTurb][fast::np1].vel_force[iNode*3+0] << " " << std::endl ;
            }
            actuatorVelFile.close() ;
        }
        
        // Do the blades first
        int nBlades = get_numBladesLoc(iTurb);
        for(int iBlade=0; iBlade < nBlades; iBlade++) {
            
            // Create interpolating parameter - Distance from hub
            int nForcePtsBlade = get_numForcePtsBladeLoc(iTurb);
            std::vector<double> rDistForce(nForcePtsBlade) ;
            for(int j=0; j < nForcePtsBlade; j++) {
                int iNodeForce = 1 + iBlade * nForcePtsBlade + j ; //The number of actuator force points is always the same for all blades
                rDistForce[j] = sqrt( 
                    (extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[0])*(extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[0])  
                    + (extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[0])*(extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[0])  
                    + (extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[0])*(extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[0])  			
                    );
            }
            
            // Interpolate to the velocity nodes
            int nVelPtsBlade = get_numVelPtsBladeLoc(iTurb);
            for(int j=0; j < nVelPtsBlade; j++) {
                int iNodeVel = 1 + iBlade * nVelPtsBlade + j ; //Assumes the same number of velocity (Aerodyn) nodes for all blades
                double rDistVel = sqrt( 
                    (extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[0])*(extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[0])  
                    + (extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[0])*(extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[0])  
                    + (extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[0])*(extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[0])  			
                    );
                //Find nearest two force nodes
                int jForceLower = 0;
                while ( (rDistForce[jForceLower+1] < rDistVel) && ( jForceLower < (nForcePtsBlade-2)) )   {
                    jForceLower = jForceLower + 1;
                }
                int iNodeForceLower = 1 + iBlade * nForcePtsBlade + jForceLower ; 
                double rInterp = (rDistVel - rDistForce[jForceLower])/(rDistForce[jForceLower+1]-rDistForce[jForceLower]);
                for (int k=0; k < 3; k++) {
                    double tmp = velForceNodeData[iTurb][fast::np1].vel_force[iNodeForceLower*3+k] + rInterp * (velForceNodeData[iTurb][fast::np1].vel_force[(iNodeForceLower+1)*3+k] - velForceNodeData[iTurb][fast::np1].vel_force[iNodeForceLower*3+k]);
                    velForceNodeData[iTurb][fast::np1].vel_vel_resid += (velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+k] - tmp)*(velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+k] - tmp);
                    velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+k] = tmp;
                }
            }
        }
        
        // Now the tower if present and used
        int nVelPtsTower = get_numVelPtsTwrLoc(iTurb);
        if ( nVelPtsTower > 0 ) {
            
            // Create interpolating parameter - Distance from first node from ground
            int nForcePtsTower = get_numForcePtsTwrLoc(iTurb);
            std::vector<double> hDistForce(nForcePtsTower) ;
            int iNodeBotTowerForce = 1 + nBlades * get_numForcePtsBladeLoc(iTurb); // The number of actuator force points is always the same for all blades
            for(int j=0; j < nForcePtsTower; j++) {
                int iNodeForce = iNodeBotTowerForce + j ; 
                hDistForce[j] = sqrt( 
                    (extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[iNodeBotTowerForce])*(extinfw_i_f_FAST[iTurb].pxForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pxForce[iNodeBotTowerForce])  
                    + (extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[iNodeBotTowerForce])*(extinfw_i_f_FAST[iTurb].pyForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pyForce[iNodeBotTowerForce])  
                    + (extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[iNodeBotTowerForce])*(extinfw_i_f_FAST[iTurb].pzForce[iNodeForce] - extinfw_i_f_FAST[iTurb].pzForce[iNodeBotTowerForce])	
                    );
            }
            
            
            int iNodeBotTowerVel = 1 + nBlades * get_numVelPtsBladeLoc(iTurb); // Assumes the same number of velocity (Aerodyn) nodes for all blades
            for(int j=0; j < nVelPtsTower; j++) {
                int iNodeVel = iNodeBotTowerVel + j ; 
                double hDistVel = sqrt( 
                    (extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[iNodeBotTowerVel])*(extinfw_i_f_FAST[iTurb].pxVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pxVel[iNodeBotTowerVel])  
                    + (extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[iNodeBotTowerVel])*(extinfw_i_f_FAST[iTurb].pyVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pyVel[iNodeBotTowerVel])  
                    + (extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[iNodeBotTowerVel])*(extinfw_i_f_FAST[iTurb].pzVel[iNodeVel] - extinfw_i_f_FAST[iTurb].pzVel[iNodeBotTowerVel])  			
                    );
                //Find nearest two force nodes
                int jForceLower = 0;
                while ( (hDistForce[jForceLower+1] < hDistVel) && ( jForceLower < (nForcePtsTower-2)) )   {
                    jForceLower = jForceLower + 1;
                }
                int iNodeForceLower = iNodeBotTowerForce + jForceLower ; 
                double rInterp = (hDistVel - hDistForce[jForceLower])/(hDistForce[jForceLower+1]-hDistForce[jForceLower]);
                for (int k=0; k < 3; k++) {
                    double tmp = velForceNodeData[iTurb][fast::np1].vel_force[iNodeForceLower*3+k] + rInterp * (velForceNodeData[iTurb][fast::np1].vel_force[(iNodeForceLower+1)*3+k] - velForceNodeData[iTurb][fast::np1].vel_force[iNodeForceLower*3+k]);
                    velForceNodeData[iTurb][fast::np1].vel_vel_resid += (velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+k] - tmp)*(velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+k] - tmp);
                    velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+k] = tmp;
                }
            }
        }
        
    } // End loop over turbines
    
}

void fast::OpenFAST::computeTorqueThrust(int iTurbGlob, std::vector<double> & torque, std::vector<double> & thrust) {
    
    //Compute the torque and thrust based on the forces at the actuator nodes
    std::vector<double> relLoc(3,0.0);
    std::vector<double> rPerpShft(3);
    thrust[0] = 0.0; thrust[1] = 0.0; thrust[2] = 0.0;
    torque[0] = 0.0; torque[1] = 0.0; torque[2] = 0.0;    
    
    std::vector<double> hubShftVec(3);
    getHubShftDir(hubShftVec, iTurbGlob, fast::np1);
    
    int iTurbLoc = get_localTurbNo(iTurbGlob) ;
    int nfpts = get_numForcePtsBlade(iTurbLoc);
    for (int k=0; k < get_numBladesLoc(iTurbLoc); k++) {
        for (int j=0; j < nfpts; j++) {
            int iNode = 1 + nfpts*k + j ;
            
            thrust[0] = thrust[0] + extinfw_i_f_FAST[iTurbLoc].fx[iNode] ;
            thrust[1] = thrust[1] + extinfw_i_f_FAST[iTurbLoc].fy[iNode] ;
            thrust[2] = thrust[2] + extinfw_i_f_FAST[iTurbLoc].fz[iNode] ;
            
            relLoc[0] = extinfw_i_f_FAST[iTurbLoc].pxForce[iNode] - extinfw_i_f_FAST[iTurbLoc].pxForce[0] ;
            relLoc[1] = extinfw_i_f_FAST[iTurbLoc].pyForce[iNode] - extinfw_i_f_FAST[iTurbLoc].pyForce[0];
            relLoc[2] = extinfw_i_f_FAST[iTurbLoc].pzForce[iNode] - extinfw_i_f_FAST[iTurbLoc].pzForce[0];            
            
	    double rDotHubShftVec = relLoc[0]*hubShftVec[0] + relLoc[1]*hubShftVec[1] + relLoc[2]*hubShftVec[2]; 
	    for (int j=0; j < 3; j++)  rPerpShft[j] = relLoc[j] - rDotHubShftVec * hubShftVec[j];
            
            torque[0] = torque[0] + rPerpShft[1] * extinfw_i_f_FAST[iTurbLoc].fz[iNode] - rPerpShft[2] * extinfw_i_f_FAST[iTurbLoc].fy[iNode] + extinfw_i_f_FAST[iTurbLoc].momentx[iNode] ;
            torque[1] = torque[1] + rPerpShft[2] * extinfw_i_f_FAST[iTurbLoc].fx[iNode] - rPerpShft[0] * extinfw_i_f_FAST[iTurbLoc].fz[iNode] + extinfw_i_f_FAST[iTurbLoc].momenty[iNode] ;
            torque[2] = torque[2] + rPerpShft[0] * extinfw_i_f_FAST[iTurbLoc].fy[iNode] - rPerpShft[1] * extinfw_i_f_FAST[iTurbLoc].fx[iNode] + extinfw_i_f_FAST[iTurbLoc].momentz[iNode] ;
            
        }
    }
}

fast::ActuatorNodeType fast::OpenFAST::getVelNodeType(int iTurbGlob, int iNode) {
    // Return the type of velocity node for the given node number. The node ordering (from FAST) is 
    // Node 0 - Hub node
    // Blade 1 nodes
    // Blade 2 nodes
    // Blade 3 nodes
    // Tower nodes
    
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbGlob);
    if (iNode) {
        if ( (iNode + 1 - (get_numVelPts(iTurbLoc) - get_numVelPtsTwr(iTurbLoc)) ) > 0) {
            return TOWER; 
        }
        else {
            return BLADE;
        }
    }
    else {
        return HUB; 
    }
    
}

fast::ActuatorNodeType fast::OpenFAST::getForceNodeType(int iTurbGlob, int iNode) {
    // Return the type of actuator force node for the given node number. The node ordering (from FAST) is 
    // Node 0 - Hub node
    // Blade 1 nodes
    // Blade 2 nodes
    // Blade 3 nodes
    // Tower nodes
    
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbGlob);
    if (iNode) {
        if ( (iNode + 1 - (get_numForcePts(iTurbLoc) - get_numForcePtsTwr(iTurbLoc)) ) > 0) {
            return TOWER; 
        }
        else {
            return BLADE;
        }
    }
    else {
        return HUB; 
    }
    
}

void fast::OpenFAST::getBladeRefPositions(std::vector<double> & bldRefPos, int iTurbGlob) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j<nPtsBlade; j++) {
            for (int k=0; k < 6; k++) {
                bldRefPos[iRunTot*6+k] = brFSIData[iTurbLoc][fast::np1].bld_ref_pos[iRunTot*6+k];
            }
            iRunTot++;
        }
    }
    
}

void fast::OpenFAST::getBladeDisplacements(std::vector<double> & bldDefl, std::vector<double> & bldVel, int iTurbGlob, fast::timeStep t) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j<nPtsBlade; j++) {
            for (int k=0; k < 6; k++) {
                bldDefl[iRunTot*6+k] = brFSIData[iTurbLoc][t].bld_def[iRunTot*6+k];
                bldVel[iRunTot*6+k] = brFSIData[iTurbLoc][t].bld_vel[iRunTot*6+k];
            }
            iRunTot++;
        }
    }
    
}

void fast::OpenFAST::getTowerRefPositions(std::vector<double> & twrRefPos, int iTurbGlob) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    for (int i=0; i < nPtsTwr; i++) {
        for (int j=0; j < 6; j++) {
            twrRefPos[i*6+j] = brFSIData[iTurbLoc][fast::np1].twr_ref_pos[i*6+j];
        }
    }
    
}

void fast::OpenFAST::getTowerDisplacements(std::vector<double> & twrDefl, std::vector<double> & twrVel, int iTurbGlob, fast::timeStep t) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    for (int i=0; i < nPtsTwr; i++) {
        for (int j=0; j < 6; j++) {
            twrDefl[i*6+j] = brFSIData[iTurbLoc][t].twr_def[i*6+j];
            twrVel[i*6+j] = brFSIData[iTurbLoc][t].twr_vel[i*6+j];            
        }
    }
    
}

void fast::OpenFAST::getHubRefPosition(std::vector<double> & hubRefPos, int iTurbGlob) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < 6; j++) 
        hubRefPos[j] = brFSIData[iTurbLoc][fast::np1].hub_ref_pos[j];
    
}

void fast::OpenFAST::getHubDisplacement(std::vector<double> & hubDefl, std::vector<double> & hubVel, int iTurbGlob, fast::timeStep t) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < 6; j++) {
        hubDefl[j] = brFSIData[iTurbLoc][t].hub_def[j];
        hubVel[j] = brFSIData[iTurbLoc][t].hub_vel[j];            
    }
    
}

void fast::OpenFAST::getNacelleRefPosition(std::vector<double> & nacRefPos, int iTurbGlob) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < 6; j++) 
        nacRefPos[j] = brFSIData[iTurbLoc][fast::np1].nac_ref_pos[j];
    
}
    

void fast::OpenFAST::getNacelleDisplacement(std::vector<double> & nacDefl, std::vector<double> & nacVel, int iTurbGlob, fast::timeStep t) {
    
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for (int j=0; j < 6; j++) {
        nacDefl[j] = brFSIData[iTurbLoc][t].nac_def[j];
        nacVel[j] = brFSIData[iTurbLoc][t].nac_vel[j];            
    }
    
}

void fast::OpenFAST::setBladeForces(std::vector<double> & bldForces, int iTurbGlob, fast::timeStep t) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nBlades = get_numBladesLoc(iTurbLoc);
    int iRunTot = 0;
    for (int i=0; i < nBlades; i++) {
        int nPtsBlade = turbineData[iTurbLoc].nBRfsiPtsBlade[i];
        for(int j=0; j < nPtsBlade; j++) {
            for(int k=0; k < 6; k++) {
                brFSIData[iTurbLoc][t].bld_ld[6*iRunTot+k] = bldForces[6*iRunTot+k];
            }
            iRunTot++;
        }
    }
    
    //TODO: May be calculate the residual as well. 
}

void fast::OpenFAST::setTowerForces(std::vector<double> & twrForces, int iTurbGlob, fast::timeStep t) {

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    int nPtsTwr = turbineData[iTurbLoc].nBRfsiPtsTwr;
    for (int i=0; i < nPtsTwr; i++)
        for (int j=0; j < 6; j++) 
            brFSIData[iTurbLoc][t].twr_ld[i*6+j] = twrForces[i*6+j];
    //TODO: May be calculate the residual as well.     
    
}

//! Sets a uniform X force at all blade nodes 
void fast::OpenFAST::setUniformXBladeForces() {
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        int iTurbGlob = turbineMapProcToGlob[iTurb];
        int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
        std::vector<double> fsiForceTower(6*nPtsTwr,0.0);
        setTowerForces(fsiForceTower, iTurbGlob, fast::np1);
            
        int nBlades = get_numBladesLoc(iTurb);
        int nTotPtsBlade = 0;
        for(int iBlade=0; iBlade < nBlades; iBlade++)
            nTotPtsBlade += turbineData[iTurb].nBRfsiPtsBlade[iBlade];
        std::vector<double> fsiForceBlade(6*nTotPtsBlade, 0.0);
        for(int i=0; i < nTotPtsBlade; i++)
            fsiForceBlade[i*6] = 5e-4; // X component of force
        setBladeForces(fsiForceBlade, iTurbGlob, fast::np1);
        
    }
}

void fast::OpenFAST::allocateMemory() {
    
    for (int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
        if (dryRun) {
            if(worldMPIRank == 0) {
                std::cout << "iTurb = " << iTurb << " turbineMapGlobToProc[iTurb] = " << turbineMapGlobToProc[iTurb] << std::endl ;
            }
        }
        if(worldMPIRank == turbineMapGlobToProc[iTurb]) {
            turbineMapProcToGlob[nTurbinesProc] = iTurb;
            reverseTurbineMapProcToGlob[iTurb] = nTurbinesProc;
            nTurbinesProc++ ;
        }
        turbineSetProcs.insert(turbineMapGlobToProc[iTurb]);
    }
    
    int nProcsWithTurbines=0;
    turbineProcs.resize(turbineSetProcs.size());
    
    for (std::set<int>::const_iterator p = turbineSetProcs.begin(); p != turbineSetProcs.end(); p++) {
        turbineProcs[nProcsWithTurbines] = *p;
        nProcsWithTurbines++ ;
    }

    if (dryRun) {
        if (nTurbinesProc > 0) {
            std::ofstream turbineAllocFile;
            turbineAllocFile.open("turbineAlloc." + std::to_string(worldMPIRank) + ".txt") ;
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                turbineAllocFile << "Proc " << worldMPIRank << " loc iTurb " << iTurb << " glob iTurb " << turbineMapProcToGlob[iTurb] << std::endl ;
            }
            turbineAllocFile.flush();
            turbineAllocFile.close() ;
        }
    }
    
    
    // Construct a group containing all procs running atleast 1 turbine in FAST
    MPI_Group_incl(worldMPIGroup, nProcsWithTurbines, &turbineProcs[0], &fastMPIGroup) ;
    int fastMPIcommTag = MPI_Comm_create(mpiComm, fastMPIGroup, &fastMPIComm);
    if (MPI_COMM_NULL != fastMPIComm) {
        MPI_Comm_rank(fastMPIComm, &fastMPIRank);
    }

    turbineData.resize(nTurbinesProc);
    velForceNodeData.resize(nTurbinesProc);
    brFSIData.resize(nTurbinesProc);
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        turbineData[iTurb].TurbineBasePos.resize(3);
        
        int iTurbGlob = turbineMapProcToGlob[iTurb];
        turbineData[iTurb].sType = globTurbineData[iTurbGlob].sType;
        turbineData[iTurb].TurbID = globTurbineData[iTurbGlob].TurbID;
        turbineData[iTurb].FASTInputFileName = globTurbineData[iTurbGlob].FASTInputFileName ;
        turbineData[iTurb].FASTRestartFileName = globTurbineData[iTurbGlob].FASTRestartFileName ;
        for(int i=0;i<3;i++)
            turbineData[iTurb].TurbineBasePos[i] = globTurbineData[iTurbGlob].TurbineBasePos[i];
        turbineData[iTurb].numForcePtsBlade = globTurbineData[iTurbGlob].numForcePtsBlade;
        turbineData[iTurb].numForcePtsTwr = globTurbineData[iTurbGlob].numForcePtsTwr;

        // To hold data for 4 time steps
        velForceNodeData[iTurb].resize(4); 
        brFSIData[iTurb].resize(4); 
        
    }

    // Allocate memory for Turbine datastructure for all turbines
    FAST_AllocateTurbines(&nTurbinesProc, &ErrStat, ErrMsg);

    // Allocate memory for ExtInfw Input types in FAST
    extinfw_i_f_FAST.resize(nTurbinesProc) ;
    extinfw_o_t_FAST.resize(nTurbinesProc) ;

    // Allocate memory for ExtLd Input types in FAST
    extld_i_f_FAST.resize(nTurbinesProc) ;
    extld_o_t_FAST.resize(nTurbinesProc) ;

    sc_i_f_FAST.resize(nTurbinesProc) ;
    sc_o_t_FAST.resize(nTurbinesProc) ;
    
}

void fast::OpenFAST::allocateMemory2(int iTurbLoc) {

    if (turbineData[iTurbLoc].sType == EXTINFLOW) {
        turbineData[iTurbLoc].nBRfsiPtsBlade = std::vector<int>(turbineData[iTurbLoc].numBlades,0);
        turbineData[iTurbLoc].nBRfsiPtsTwr = 0;
    } else if (turbineData[iTurbLoc].sType == EXTLOADS) {
        turbineData[iTurbLoc].nBRfsiPtsBlade.resize(turbineData[iTurbLoc].numBlades);
        int nTotBldNds = 0;
        for(int i=0; i < turbineData[iTurbLoc].numBlades; i++) {
            nTotBldNds += extld_i_f_FAST[iTurbLoc].nBladeNodes[i];
            turbineData[iTurbLoc].nBRfsiPtsBlade[i] = extld_i_f_FAST[iTurbLoc].nBladeNodes[i];
        }
        turbineData[iTurbLoc].nBRfsiPtsTwr = extld_i_f_FAST[iTurbLoc].nTowerNodes[0];
        std::cout << "Number of tower nodes = " << turbineData[iTurbLoc].nBRfsiPtsTwr;

        // Allocate memory for reference position only for 1 time step - np1
        brFSIData[iTurbLoc][3].twr_ref_pos.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
        brFSIData[iTurbLoc][3].bld_ref_pos.resize(6*nTotBldNds);
        brFSIData[iTurbLoc][3].hub_ref_pos.resize(6);
        brFSIData[iTurbLoc][3].nac_ref_pos.resize(6);
        for(int k=0; k<4; k++) {
            brFSIData[iTurbLoc][k].twr_def.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
            brFSIData[iTurbLoc][k].twr_vel.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
            brFSIData[iTurbLoc][k].bld_def.resize(6*nTotBldNds);
            brFSIData[iTurbLoc][k].bld_vel.resize(6*nTotBldNds);
            brFSIData[iTurbLoc][k].twr_ld.resize(6*turbineData[iTurbLoc].nBRfsiPtsTwr);
            brFSIData[iTurbLoc][k].bld_ld.resize(6*nTotBldNds);
            brFSIData[iTurbLoc][k].hub_def.resize(6);
            brFSIData[iTurbLoc][k].hub_vel.resize(6);
            brFSIData[iTurbLoc][k].nac_def.resize(6);
            brFSIData[iTurbLoc][k].nac_vel.resize(6);
        }
    }

    if ( (turbineData[iTurbLoc].sType == EXTINFLOW) && (turbineData[iTurbLoc].inflowType == 2) ) {
        //Inflow data is coming from external program like a CFD solver
        turbineData[iTurbLoc].numForcePts = 1 + turbineData[iTurbLoc].numForcePtsTwr + turbineData[iTurbLoc].numBlades * turbineData[iTurbLoc].numForcePtsBlade ;
        turbineData[iTurbLoc].numVelPts = 1 + turbineData[iTurbLoc].numVelPtsTwr + turbineData[iTurbLoc].numBlades * turbineData[iTurbLoc].numVelPtsBlade ;
        
        int nfpts = get_numForcePtsLoc(iTurbLoc);
        int nvelpts = get_numVelPtsLoc(iTurbLoc);
        
        for(int k=0; k<4; k++) {
            velForceNodeData[iTurbLoc][k].x_vel.resize(3*nvelpts) ;
            velForceNodeData[iTurbLoc][k].xdot_vel.resize(3*nvelpts) ;
            velForceNodeData[iTurbLoc][k].vel_vel.resize(3*nvelpts) ;
            velForceNodeData[iTurbLoc][k].x_force.resize(3*nfpts) ;
            velForceNodeData[iTurbLoc][k].xdot_force.resize(3*nfpts) ;
            velForceNodeData[iTurbLoc][k].orient_force.resize(9*nfpts) ;
            velForceNodeData[iTurbLoc][k].vel_force.resize(3*nfpts) ;
            velForceNodeData[iTurbLoc][k].force.resize(3*nfpts) ;
        }
        
        if ( isDebug() ) {
            for (int iNode=0; iNode < get_numVelPtsLoc(iTurbLoc); iNode++) {
                std::cout << "Node " << iNode << " Position = " << extinfw_i_f_FAST[iTurbLoc].pxVel[iNode] << " " << extinfw_i_f_FAST[iTurbLoc].pyVel[iNode] << " " << extinfw_i_f_FAST[iTurbLoc].pzVel[iNode] << " " << std::endl ;
            }
        }
        
    } else {
        // Inflow data is coming from inflow module
        turbineData[iTurbLoc].numForcePtsTwr = 0;
        turbineData[iTurbLoc].numForcePtsBlade = 0;
        turbineData[iTurbLoc].numForcePts = 0;
        turbineData[iTurbLoc].numVelPtsTwr = 0;
        turbineData[iTurbLoc].numVelPtsBlade = 0;
        turbineData[iTurbLoc].numVelPts = 0;
    }

}

void fast::OpenFAST::allocateTurbinesToProcsSimple() {
    
    // Allocate turbines to each processor - round robin fashion
    int nProcs ;
    MPI_Comm_size(mpiComm, &nProcs);
    for(int j = 0; j < nTurbinesGlob; j++)  turbineMapGlobToProc[j] = j % nProcs ;
    
}

void fast::OpenFAST::end() {
    
    // Deallocate types we allocated earlier
    
    if (nTurbinesProc > 0) closeVelocityDataFile(nt_global, velNodeDataFile);
    
    if ( !dryRun) {
        bool stopTheProgram = false;
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_End(&iTurb, &stopTheProgram);
        }
    }
    
    MPI_Group_free(&fastMPIGroup);
    if (MPI_COMM_NULL != fastMPIComm) {
        MPI_Comm_free(&fastMPIComm);
    }
    MPI_Group_free(&worldMPIGroup);
    
    if(scStatus) {
        
        destroy_SuperController(sc) ;
        
        if(scLibHandle != NULL) {
            // close the library
            std::cout << "Closing library...\n";
            dlclose(scLibHandle);
        }
        
    }
    
}

void fast::OpenFAST::readVelocityData(int nTimesteps) {
    
    int nTurbines;
    
    hid_t velDataFile = H5Fopen(("velDatafile." + std::to_string(worldMPIRank) + ".h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    
    {
        hid_t attr = H5Aopen(velDataFile, "nTurbines", H5P_DEFAULT);
        herr_t ret = H5Aread(attr, H5T_NATIVE_INT, &nTurbines) ;
        H5Aclose(attr);
        
    }
    
    // Allocate memory and read the velocity data. 
    velNodeData.resize(nTurbines);
    for (int iTurb=0; iTurb < nTurbines; iTurb++) {
        int nVelPts = get_numVelPtsLoc(iTurb) ;
        velNodeData[iTurb].resize(nTimesteps*nVelPts*6) ;
        hid_t dset_id = H5Dopen2(velDataFile, ("/turbine" + std::to_string(turbineMapProcToGlob[iTurb])).c_str(), H5P_DEFAULT);
        hid_t dspace_id = H5Dget_space(dset_id);
        
        hsize_t start[3]; start[1] = 0; start[2] = 0;
        hsize_t count[3]; count[0] = 1; count[1] = nVelPts; count[2] = 6;
        hid_t mspace_id = H5Screate_simple(3, count, NULL); 
        
        for (int iStep=0; iStep < nTimesteps; iStep++) {
            start[0] = iStep;
            H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
            herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mspace_id, dspace_id, H5P_DEFAULT, &velNodeData[iTurb][iStep*nVelPts*6] );
        }
        herr_t status = H5Dclose(dset_id);
        
        
    }
    
}

hid_t fast::OpenFAST::openVelocityDataFile(bool createFile) {
    
    hid_t velDataFile;
    if (createFile) {
        // Open the file in create mode
        velDataFile = H5Fcreate(("velDatafile." + std::to_string(worldMPIRank) + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        
        {
            hsize_t dims[1];
            dims[0] = 1;
            hid_t dataSpace = H5Screate_simple(1, dims, NULL);
            hid_t attr = H5Acreate2(velDataFile, "nTurbines", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            herr_t status = H5Awrite(attr, H5T_NATIVE_INT, &nTurbinesProc);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);
            
            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(velDataFile, "nTimesteps", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);
        }
        
        int ntMax = tMax/dtFAST ;
        
        for (int iTurb = 0; iTurb < nTurbinesProc; iTurb++) {
            int nVelPts = get_numVelPtsLoc(iTurb);
            hsize_t dims[3];
            dims[0] = ntMax; dims[1] = nVelPts; dims[2] = 6 ;
            
            hsize_t chunk_dims[3];
            chunk_dims[0] = 1; chunk_dims[1] = nVelPts; chunk_dims[2] = 6;
            hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(dcpl_id, 3, chunk_dims);
            
            hid_t dataSpace = H5Screate_simple(3, dims, NULL);
            hid_t dataSet = H5Dcreate(velDataFile, ("/turbine" + std::to_string(turbineMapProcToGlob[iTurb])).c_str(), H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);    
            
            herr_t status = H5Pclose(dcpl_id);
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
        }
        
    } else {
        // Open the file in append mode
        velDataFile = H5Fopen(("velDatafile." + std::to_string(worldMPIRank) + ".h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    
    return velDataFile;
    
}

herr_t fast::OpenFAST::closeVelocityDataFile(int nt_global, hid_t velDataFile) {
    
    herr_t status = H5Fclose(velDataFile) ;
    return status;
}


void fast::OpenFAST::backupVelocityDataFile(int curTimeStep, hid_t & velDataFile) {

    closeVelocityDataFile(curTimeStep, velDataFile);
        
    std::ifstream source("velDatafile." + std::to_string(worldMPIRank) + ".h5", std::ios::binary);
    std::ofstream dest("velDatafile." + std::to_string(worldMPIRank) + ".h5." + std::to_string(curTimeStep) + ".bak", std::ios::binary);

    dest << source.rdbuf();
    source.close();
    dest.close();

    velDataFile = openVelocityDataFile(false);
}

void fast::OpenFAST::writeVelocityData(hid_t h5File, int iTurb, int iTimestep, ExtInfw_InputType_t iData, ExtInfw_OutputType_t oData) {
    
    hsize_t start[3]; start[0] = iTimestep; start[1] = 0; start[2] = 0;
    int nVelPts = get_numVelPtsLoc(iTurb) ;
    hsize_t count[3]; count[0] = 1; count[1] = nVelPts; count[2] = 6;
    
    std::vector<double> tmpVelData;
    tmpVelData.resize(nVelPts * 6);
    
    for (int iNode=0 ; iNode < nVelPts; iNode++) {
        tmpVelData[iNode*6 + 0] = iData.pxVel[iNode];
        tmpVelData[iNode*6 + 1] = iData.pyVel[iNode];
        tmpVelData[iNode*6 + 2] = iData.pzVel[iNode];
        tmpVelData[iNode*6 + 3] = oData.u[iNode];
        tmpVelData[iNode*6 + 4] = oData.v[iNode];
        tmpVelData[iNode*6 + 5] = oData.w[iNode];
    }
    
    hid_t dset_id = H5Dopen2(h5File, ("/turbine" + std::to_string(turbineMapProcToGlob[iTurb])).c_str(), H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t mspace_id = H5Screate_simple(3, count, NULL);  
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mspace_id, dspace_id, H5P_DEFAULT, tmpVelData.data());
    
    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Sclose(mspace_id);
    
    hid_t attr_id = H5Aopen_by_name(h5File, ".", "nTimesteps", H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr_id, H5T_NATIVE_INT, &iTimestep);
    status = H5Aclose(attr_id);
    
}

void fast::OpenFAST::applyVelocityData(int iPrestart, int iTurb, ExtInfw_OutputType_t extinfw_o_t_FAST, std::vector<double> & velData) {
    
    int nVelPts = get_numVelPtsLoc(iTurb);
    for (int j = 0; j < nVelPts; j++){
        extinfw_o_t_FAST.u[j] = velData[(iPrestart*nVelPts+j)*6 + 3]; 
        extinfw_o_t_FAST.v[j] = velData[(iPrestart*nVelPts+j)*6 + 4];
        extinfw_o_t_FAST.w[j] = velData[(iPrestart*nVelPts+j)*6 + 5];
    }
    
}

void fast::OpenFAST::send_data_to_openfast(fast::timeStep t) {
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        if (turbineData[iTurb].inflowType == 2) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            for (int iNodeVel=0; iNodeVel < nvelpts; iNodeVel++) {
                extinfw_o_t_FAST[iTurb].u[iNodeVel] = velForceNodeData[iTurb][t].vel_vel[iNodeVel*3+0];
                extinfw_o_t_FAST[iTurb].v[iNodeVel] = velForceNodeData[iTurb][t].vel_vel[iNodeVel*3+1];
                extinfw_o_t_FAST[iTurb].w[iNodeVel] = velForceNodeData[iTurb][t].vel_vel[iNodeVel*3+2];
            }
        }

        if(turbineData[iTurb].sType == EXTLOADS) {

            int nBlades = turbineData[iTurb].numBlades;
            int iRunTot = 0;
            for(int i=0; i < nBlades; i++) {
                int nPtsBlade = turbineData[iTurb].nBRfsiPtsBlade[i];
                for (int j=0; j < nPtsBlade; j++) {
                    for (int k=0; k<6; k++) {
                        extld_o_t_FAST[iTurb].bldLd[iRunTot*6+k] = brFSIData[iTurb][t].bld_ld[i];
                    }
                    iRunTot++;
                }
            }
            
            int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
            for (int i=0; i < nPtsTwr*6; i++)
                extld_o_t_FAST[iTurb].twrLd[i] = brFSIData[iTurb][t].twr_ld[i];

        }
        
    }
}


void fast::OpenFAST::send_data_to_openfast(double ss_time) {
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        if (turbineData[iTurb].inflowType == 2) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            for (int iNodeVel=0; iNodeVel < nvelpts; iNodeVel++) {
                extinfw_o_t_FAST[iTurb].u[iNodeVel] = velForceNodeData[iTurb][fast::n].vel_vel[iNodeVel*3+0] + ss_time * (velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+0] - velForceNodeData[iTurb][fast::n].vel_vel[iNodeVel*3+0]);
                extinfw_o_t_FAST[iTurb].v[iNodeVel] = velForceNodeData[iTurb][fast::n].vel_vel[iNodeVel*3+1] + ss_time * (velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+1] - velForceNodeData[iTurb][fast::n].vel_vel[iNodeVel*3+1]);
                extinfw_o_t_FAST[iTurb].w[iNodeVel] = velForceNodeData[iTurb][fast::n].vel_vel[iNodeVel*3+2] + ss_time * (velForceNodeData[iTurb][fast::np1].vel_vel[iNodeVel*3+2] - velForceNodeData[iTurb][fast::n].vel_vel[iNodeVel*3+2]);
            }
        }
    }
}

void fast::OpenFAST::get_data_from_openfast(timeStep t) {
    
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        if (turbineData[iTurb].inflowType == 2) {
            int nvelpts = get_numVelPtsLoc(iTurb);
            int nfpts = get_numForcePtsLoc(iTurb);
            for (int i=0; i<nvelpts; i++) {
                velForceNodeData[iTurb][t].x_vel_resid += (velForceNodeData[iTurb][t].x_vel[i*3+0] - extinfw_i_f_FAST[iTurb].pxVel[i])*(velForceNodeData[iTurb][t].x_vel[i*3+0] - extinfw_i_f_FAST[iTurb].pxVel[i]);
                velForceNodeData[iTurb][t].x_vel[i*3+0] = extinfw_i_f_FAST[iTurb].pxVel[i];
                velForceNodeData[iTurb][t].x_vel_resid += (velForceNodeData[iTurb][t].x_vel[i*3+1] - extinfw_i_f_FAST[iTurb].pyVel[i])*(velForceNodeData[iTurb][t].x_vel[i*3+1] - extinfw_i_f_FAST[iTurb].pyVel[i]);
                velForceNodeData[iTurb][t].x_vel[i*3+1] = extinfw_i_f_FAST[iTurb].pyVel[i];
                velForceNodeData[iTurb][t].x_vel_resid += (velForceNodeData[iTurb][t].x_vel[i*3+2] - extinfw_i_f_FAST[iTurb].pzVel[i])*(velForceNodeData[iTurb][t].x_vel[i*3+2] - extinfw_i_f_FAST[iTurb].pzVel[i]);
                velForceNodeData[iTurb][t].x_vel[i*3+2] = extinfw_i_f_FAST[iTurb].pzVel[i];
                velForceNodeData[iTurb][t].xdot_vel_resid += (velForceNodeData[iTurb][t].xdot_vel[i*3+0] - extinfw_i_f_FAST[iTurb].pxdotVel[i])*(velForceNodeData[iTurb][t].xdot_vel[i*3+0] - extinfw_i_f_FAST[iTurb].pxdotVel[i]);
                velForceNodeData[iTurb][t].xdot_vel[i*3+0] = extinfw_i_f_FAST[iTurb].pxdotVel[i];
                velForceNodeData[iTurb][t].xdot_vel_resid += (velForceNodeData[iTurb][t].xdot_vel[i*3+1] - extinfw_i_f_FAST[iTurb].pydotVel[i])*(velForceNodeData[iTurb][t].xdot_vel[i*3+1] - extinfw_i_f_FAST[iTurb].pydotVel[i]);
                velForceNodeData[iTurb][t].xdot_vel[i*3+1] = extinfw_i_f_FAST[iTurb].pydotVel[i];
                velForceNodeData[iTurb][t].xdot_vel_resid += (velForceNodeData[iTurb][t].xdot_vel[i*3+2] - extinfw_i_f_FAST[iTurb].pzdotVel[i])*(velForceNodeData[iTurb][t].xdot_vel[i*3+2] - extinfw_i_f_FAST[iTurb].pzdotVel[i]);
                velForceNodeData[iTurb][t].xdot_vel[i*3+2] = extinfw_i_f_FAST[iTurb].pzdotVel[i];
            }
            
            for (int i=0; i<nfpts; i++) {
                velForceNodeData[iTurb][t].x_force_resid += (velForceNodeData[iTurb][t].x_force[i*3+0] - extinfw_i_f_FAST[iTurb].pxForce[i])*(velForceNodeData[iTurb][t].x_force[i*3+0] - extinfw_i_f_FAST[iTurb].pxForce[i]);
                velForceNodeData[iTurb][t].x_force[i*3+0] = extinfw_i_f_FAST[iTurb].pxForce[i];
                velForceNodeData[iTurb][t].x_force_resid += (velForceNodeData[iTurb][t].x_force[i*3+1] - extinfw_i_f_FAST[iTurb].pyForce[i])*(velForceNodeData[iTurb][t].x_force[i*3+1] - extinfw_i_f_FAST[iTurb].pyForce[i]);
                velForceNodeData[iTurb][t].x_force[i*3+1] = extinfw_i_f_FAST[iTurb].pyForce[i];
                velForceNodeData[iTurb][t].x_force_resid += (velForceNodeData[iTurb][t].x_force[i*3+2] - extinfw_i_f_FAST[iTurb].pzForce[i])*(velForceNodeData[iTurb][t].x_force[i*3+2] - extinfw_i_f_FAST[iTurb].pzForce[i]);
                velForceNodeData[iTurb][t].x_force[i*3+2] = extinfw_i_f_FAST[iTurb].pzForce[i];
                velForceNodeData[iTurb][t].xdot_force_resid += (velForceNodeData[iTurb][t].xdot_force[i*3+0] - extinfw_i_f_FAST[iTurb].pxdotForce[i])*(velForceNodeData[iTurb][t].xdot_force[i*3+0] - extinfw_i_f_FAST[iTurb].pxdotForce[i]);
                velForceNodeData[iTurb][t].xdot_force[i*3+0] = extinfw_i_f_FAST[iTurb].pxdotForce[i];
                velForceNodeData[iTurb][t].xdot_force_resid += (velForceNodeData[iTurb][t].xdot_force[i*3+1] - extinfw_i_f_FAST[iTurb].pydotForce[i])*(velForceNodeData[iTurb][t].xdot_force[i*3+1] - extinfw_i_f_FAST[iTurb].pydotForce[i]);
                velForceNodeData[iTurb][t].xdot_force[i*3+1] = extinfw_i_f_FAST[iTurb].pydotForce[i];
                velForceNodeData[iTurb][t].xdot_force_resid += (velForceNodeData[iTurb][t].xdot_force[i*3+2] - extinfw_i_f_FAST[iTurb].pzdotForce[i])*(velForceNodeData[iTurb][t].xdot_force[i*3+2] - extinfw_i_f_FAST[iTurb].pzdotForce[i]);
                velForceNodeData[iTurb][t].xdot_force[i*3+2] = extinfw_i_f_FAST[iTurb].pzdotForce[i];
                for (int j=0;j<9;j++) {
                    velForceNodeData[iTurb][t].orient_force_resid += (velForceNodeData[iTurb][t].orient_force[i*9+j] - extinfw_i_f_FAST[iTurb].pOrientation[i*9+j])*(velForceNodeData[iTurb][t].orient_force[i*9+j] - extinfw_i_f_FAST[iTurb].pOrientation[i*9+j]);
                    velForceNodeData[iTurb][t].orient_force[i*9+j] = extinfw_i_f_FAST[iTurb].pOrientation[i*9+j];                    
                }
                velForceNodeData[iTurb][t].force_resid += (velForceNodeData[iTurb][t].force[i*3+0] - extinfw_i_f_FAST[iTurb].fx[i])*(velForceNodeData[iTurb][t].force[i*3+0] - extinfw_i_f_FAST[iTurb].fx[i]);
                velForceNodeData[iTurb][t].force[i*3+0] = extinfw_i_f_FAST[iTurb].fx[i];
                velForceNodeData[iTurb][t].force_resid += (velForceNodeData[iTurb][t].force[i*3+1] - extinfw_i_f_FAST[iTurb].fy[i])*(velForceNodeData[iTurb][t].force[i*3+1] - extinfw_i_f_FAST[iTurb].fy[i]);
                velForceNodeData[iTurb][t].force[i*3+1] = extinfw_i_f_FAST[iTurb].fy[i];
                velForceNodeData[iTurb][t].force_resid += (velForceNodeData[iTurb][t].force[i*3+2] - extinfw_i_f_FAST[iTurb].fz[i])*(velForceNodeData[iTurb][t].force[i*3+2] - extinfw_i_f_FAST[iTurb].fz[i]);
                velForceNodeData[iTurb][t].force[i*3+2] = extinfw_i_f_FAST[iTurb].fz[i];
            }
        }

        if(turbineData[iTurb].sType == EXTLOADS) {

            int nBlades = turbineData[iTurb].numBlades;
            int iRunTot = 0;
            for (int i=0; i < nBlades; i++) {
                int nPtsBlade = turbineData[iTurb].nBRfsiPtsBlade[i];
                for (int j=0; j < nPtsBlade; j++) {
                    for (int k=0; k < 6; k++) {
                        brFSIData[iTurb][t].bld_def[iRunTot*6+k] = extld_i_f_FAST[iTurb].bldDef[iRunTot*12+k];
                        brFSIData[iTurb][t].bld_vel[iRunTot*6+k] = extld_i_f_FAST[iTurb].bldDef[iRunTot*12+6+k];
                    }
                    iRunTot++;
                }
            }
            
            int nPtsTwr = turbineData[iTurb].nBRfsiPtsTwr;
            for (int i=0; i < nPtsTwr; i++) {
                for (int j = 0; j < 6; j++) {
                    brFSIData[iTurb][t].twr_def[i*6+j] = extld_i_f_FAST[iTurb].twrDef[i*12+j];
                    brFSIData[iTurb][t].twr_vel[i*6+j] = extld_i_f_FAST[iTurb].twrDef[i*12+6+j];
                }
            }

            //TODO: May be calculate the residual here as well
        }
        
    }
}

void fast::OpenFAST::readRestartFile(int iTurbLoc, int n_t_global) {
    
    int nvelpts = get_numVelPtsLoc(iTurbLoc);
    int nfpts = get_numForcePtsLoc(iTurbLoc);
    int iTurbGlob = turbineMapProcToGlob[iTurbLoc];
    
    hid_t restartFile = H5Fopen(("of_cpp_" + std::to_string(n_t_global) + "_T" + std::to_string(iTurbGlob) + ".chkp.h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    
    int nvelpts_file, nfpts_file ;
    
    {
        hid_t attr = H5Aopen(restartFile, "nvelpts", H5P_DEFAULT);
        herr_t ret = H5Aread(attr, H5T_NATIVE_INT, &nvelpts_file) ;
        H5Aclose(attr);
        
        attr = H5Aopen(restartFile, "nfpts", H5P_DEFAULT);
        ret = H5Aread(attr, H5T_NATIVE_INT, &nfpts_file) ;
        H5Aclose(attr);
    }
    
    if( (nvelpts != nvelpts_file) || (nfpts != nfpts_file))
        throw std::runtime_error("Number of velocity or force nodes from restart file does not match input.");
    
    if (nvelpts > 0) {  
        for (int j=0; j < 4; j++) {  // Loop over states - NM2, NM1, N, NP1

            std::string gName = "/data/" + std::to_string(j);
            hid_t group_id = H5Gopen2(restartFile, gName.c_str(), H5P_DEFAULT);

            hid_t dataSet = H5Dopen2(group_id, "x_vel", H5P_DEFAULT);
            herr_t status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].x_vel.data());
            status = H5Dclose(dataSet);
            
            dataSet = H5Dopen2(group_id, "xdot_vel", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].xdot_vel.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(group_id, "vel_vel", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].xdot_vel.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(group_id, "x_force", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].xdot_vel.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(group_id, "xdot_force", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].xdot_force.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(group_id, "vel_force", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].xdot_force.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(group_id, "force", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].force.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(group_id, "orient_force", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].orient_force.data());
            status = H5Dclose(dataSet);
        }
        
    }
    
    herr_t status = H5Fclose(restartFile);
    
}


void fast::OpenFAST::writeRestartFile(int iTurbLoc, int n_t_global) {
    
    /* // HDF5 stuff to write states to restart file or read back from it */
    
    int iTurbGlob = turbineMapProcToGlob[iTurbLoc];
    int nvelpts = get_numVelPtsLoc(iTurbLoc);
    int nfpts = get_numForcePtsLoc(iTurbLoc);
    hid_t restartFile = H5Fcreate(("of_cpp_" + std::to_string(n_t_global) + "_T" + std::to_string(iTurbGlob) + ".chkp.h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    {
        hsize_t dims[1];
        dims[0] = 1;
        hid_t dataSpace = H5Screate_simple(1, dims, NULL);
        hid_t attr = H5Acreate2(restartFile, "nvelpts", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
        herr_t status = H5Awrite(attr, H5T_NATIVE_INT, &nvelpts);
        status = H5Aclose(attr);
        status = H5Sclose(dataSpace);
        
        dataSpace = H5Screate_simple(1, dims, NULL);
        attr = H5Acreate2(restartFile, "nfpts", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
        status = H5Awrite(attr, H5T_NATIVE_INT, &nfpts);
        status = H5Aclose(attr);
        status = H5Sclose(dataSpace);
        
    }
    
    if (nvelpts > 0) {

        /* Create groups */
        hid_t group_id = H5Gcreate2(restartFile, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        herr_t status = H5Gclose(group_id);
        for (int j=0; j < 4; j++) { // Loop over states - NM2, NM1, N, NP1
            std::string gName = "/data/" + std::to_string(j);
            group_id = H5Gcreate2(restartFile, gName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            hsize_t dims[1];
            dims[0] = nvelpts*3;
            hid_t dataSpace = H5Screate_simple(1, dims, NULL);
            hid_t dataSet = H5Dcreate2(group_id, "x_vel", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
            herr_t status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].x_vel.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(group_id, "xdot_vel", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
            status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].xdot_vel.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(group_id, "vel_vel", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
            status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].vel_vel.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
            dims[0] = nfpts*3;
            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(group_id, "x_force", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
            status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].x_force.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(group_id, "xdot_force", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
            status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].xdot_force.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(group_id, "vel_force", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
            status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].vel_force.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(group_id, "force", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
            status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].force.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
            dims[0] = nfpts*9;
            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(group_id, "orient_force", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velForceNodeData[iTurbLoc][0].orient_force.data());
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
            
        }
        
    }
    
    herr_t status = H5Fclose(restartFile);
    
}


void fast::OpenFAST::loadSuperController(const fast::fastInputs & fi) {
    
    if(fi.scStatus) {
        
        scStatus = fi.scStatus;
        scLibFile = fi.scLibFile;
        
        // open the library
        scLibHandle = dlopen(scLibFile.c_str(), RTLD_LAZY);
        if (!scLibHandle) {
            std::cerr << "Cannot open library: " << dlerror() << '\n';
        }
        
        create_SuperController = (create_sc_t*) dlsym(scLibHandle, "create_sc");
        // reset errors
        dlerror();
        const char *dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'create_sc': " << dlsym_error << '\n';
            dlclose(scLibHandle);
        }
        
        destroy_SuperController = (destroy_sc_t*) dlsym(scLibHandle, "destroy_sc");
        // reset errors
        dlerror();
        const char *dlsym_error_us = dlerror();
        if (dlsym_error_us) {
            std::cerr << "Cannot load symbol 'destroy_sc': " << dlsym_error_us << '\n';
            dlclose(scLibHandle);
        }
        
        sc = create_SuperController() ;
        
        numScInputs = fi.numScInputs;
        numScOutputs = fi.numScOutputs;
        
        if ( (numScInputs > 0) && (numScOutputs > 0)) {
            scOutputsGlob.resize(nTurbinesGlob*numScOutputs) ;
            scInputsGlob.resize(nTurbinesGlob*numScInputs) ;
            for (int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
                for(int iInput=0; iInput < numScInputs; iInput++) {
                    scInputsGlob[iTurb*numScInputs + iInput] = 0.0 ; // Initialize to zero
                }
                for(int iOutput=0; iOutput < numScOutputs; iOutput++) {
                    scOutputsGlob[iTurb*numScOutputs + iOutput] = 0.0 ; // Initialize to zero
                }
            }
            
        } else {
            std::cerr <<  "Make sure numScInputs and numScOutputs are greater than zero" << std::endl;
        }
        
    } else {
        scStatus = false;
        numScInputs = 0;
        numScOutputs = 0;
    }
    
}

void fast::OpenFAST::fillScInputsGlob() {
    
    // Fills the global array containing inputs to the supercontroller from all turbines
    
    for(int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
        for(int iInput=0; iInput < numScInputs; iInput++) {
            scInputsGlob[iTurb*numScInputs + iInput] = 0.0; // Initialize to zero 
        }
    }
    
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        for(int iInput=0; iInput < numScInputs; iInput++) {
            scInputsGlob[turbineMapProcToGlob[iTurb]*numScInputs + iInput] = sc_i_f_FAST[iTurb].toSC[iInput] ;
        }
    }
    
    if (MPI_COMM_NULL != fastMPIComm) {
        MPI_Allreduce(MPI_IN_PLACE, scInputsGlob.data(), numScInputs*nTurbinesGlob, MPI_DOUBLE, MPI_SUM, fastMPIComm) ;
    }
    
}


void fast::OpenFAST::fillScOutputsLoc() {
    
    // Fills the local array containing outputs from the supercontroller to each turbine
    
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        for(int iOutput=0; iOutput < numScOutputs; iOutput++) {
            sc_o_t_FAST[iTurb].fromSC[iOutput] = scOutputsGlob[turbineMapProcToGlob[iTurb]*numScOutputs + iOutput] ;
        }
    }
    
}










