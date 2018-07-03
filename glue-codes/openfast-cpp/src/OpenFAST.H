#ifndef OpenFAST_h
#define OpenFAST_h
#include "FAST_Library.h"
#include "sys/stat.h"
#include <string>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include "dlfcn.h"
//TODO: The skip MPICXX is put in place primarily to get around errors in OpenFOAM. This will cause problems if the driver program uses C++ API for MPI.
#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "mpi.h"
#include "SC.h"


namespace fast {

//! An id to indicate the type of simulation for each turbine - Simple/Actuator with optional externally specified inflow or Blade-Resolved with externally specified loads
enum simType {
    EXTINFLOW = 0,
    EXTLOADS  = 1,
    simType_END
};    
    
//! A data structure to hold all turbine related information
struct turbineDataType {
    //!Integer id for every turbine
    int TurbID;
    //! The FAST Input file name. Typically a .fst file.
    std::string FASTInputFileName;
    //! The restart/checkpoint file name. 
    std::string FASTRestartFileName;
    //! The position of the base of the turbine in global coordinates
    std::vector<float> TurbineBasePos;
    //! The approximate position of the hub of the turbine in global coordinates
    std::vector<double> TurbineHubPos;
    //! Simulation type
    simType sType;
    //! Number of blades 
    int numBlades;
    //! Number of velocity nodes (AeroDyn) per blade
    int numVelPtsBlade;
    //! Number of velocity nodes (AeroDyn) on the tower
    int numVelPtsTwr;
    //! Total number of velocity nodes (AeroDyn)
    int numVelPts;
    //! Desired number of actuator points on each blade 
    int numForcePtsBlade;
    //! Desired number of actuator points on the tower 
    int numForcePtsTwr;
    //! Total number of actuator points
    int numForcePts;
    //! Inflow Type - 1 (InflowWind) or 2 (Externally specified)
    int inflowType;
    //! Drag coefficient of nacelle
    float nacelle_cd;
    //! Frontal area of the nacelle
    float nacelle_area;
    //! Air density around this turbine
    float air_density;
    //! Number of nodes at which the forces and deflections are computed for blade-resolved FSI on each blade
    std::vector<int> nBRfsiPtsBlade;
    //! Number of nodes at which the forces and deflections are computed for blade-resolved FSI on the tower
    int nBRfsiPtsTwr;
};


//! An id to indicate whether a particular actuator point is on the hub, node or tower    
enum ActuatorNodeType {
  HUB = 0,
  BLADE = 1,
  TOWER = 2,
  ActuatorNodeType_END
};

/** An id to indicate the start type of a simulation.
 * init - Start the simulation from scratch
 * trueRestart - Restart from a checkpoint file. Code expects checkpoint files for all parts of the simulation including the controller.
 * restartDriverInitFAST - Start all turbines from scratch and use the velocity data in 'velData.h5' file to run upto desired restart time, then continue the simulation like ''trueRestart'.
 */
enum simStartType {
  init = 0,
  trueRestart = 1,
  restartDriverInitFAST = 2,
  simStartType_END    
};

//! A data structure to hold all velocity and force node information
struct turbVelForceNodeDataType {
    //! Blade location at velocity nodes
    std::vector<double> x_vel;
    //! Blade velocity at velocity nodes
    std::vector<double> xdot_vel;
    //! Sampled velocity at velocity nodes
    std::vector<double> vel_vel;
    //! Blade location at force nodes
    std::vector<double> x_force;
    //! Blade velocity at force nodes
    std::vector<double> xdot_force;
    //! Blade orientation at force nodes
    std::vector<double> orient_force;
    //! Sampled velocity at force nodes
    std::vector<double> vel_force;
    //! Actuator force at force nodes
    std::vector<double> force;
    double x_vel_resid;
    double xdot_vel_resid;
    double vel_vel_resid;
    double x_force_resid;
    double xdot_force_resid;
    double orient_force_resid;
    double vel_force_resid;
    double force_resid;
};

//! A data structure to hold all loads and deflections information for blade-resolved FSI simulations
struct turbBRfsiDataType {
    //! Tower reference position
    std::vector<double> twr_ref_pos;
    //! Tower deflections 
    std::vector<double> twr_def;
    //! Tower velocity
    std::vector<double> twr_vel;
    //! Blade reference position
    std::vector<double> bld_ref_pos;
    //! Blade deflections
    std::vector<double> bld_def;
    //! Blade velocity
    std::vector<double> bld_vel;
    //! Hub reference position    
    std::vector<double> hub_ref_pos;
    //! Hub deflections 
    std::vector<double> hub_def;
    //! Hub velocity
    std::vector<double> hub_vel;
    //! Nacelle reference position    
    std::vector<double> nac_ref_pos;
    //! Nacelle deflections 
    std::vector<double> nac_def;
    //! Nacelle velocity
    std::vector<double> nac_vel;
    
    //! Tower loads
    std::vector<double> twr_ld;
    //! Blade loads
    std::vector<double> bld_ld;
    double twr_def_resid;
    double twr_vel_resid;
    double bld_def_resid;
    double bld_vel_resid;
    double twr_ld_resid;
    double bld_ld_resid;
};

//! An enum to keep track of information stored at different time steps
enum timeStep {
    nm2 = 0,
    nm1 = 1,
    n = 2,
    np1 = 3,
    timeStep_END
};

    
    
/**
 * A class to hold all input data for a simulation run through a OpenFAST C++ glue code
 */
class fastInputs {

 public:

  //! MPI Communicator
  MPI_Comm comm;
  //! Total number of turbines in the simulation
  int nTurbinesGlob;
  //! The simulation will not run if dryRun is set to true. However, the simulation will read the input files, allocate turbines to processors and prepare to run the individual turbine instances. This flag is useful to test the setup of the simulation before running it.
  bool dryRun;
  //! Enable debug outputs if set to true
  bool debug;
  //! Start time of the simulation
  double tStart;
  //! Start type of the simulation: 'init', 'trueRestart' or 'restartDriverInitFAST'.
  simStartType simStart;
  //!Restart files will be written every so many time stneps
  int nEveryCheckPoint;
  //! Max time of the simulation
  double tMax;
  //! Time step for FAST. All turbines should have the same time step.
  double dtFAST;
  //! Number of OpenFAST time steps per unit time step of the driver program
  int nSubsteps;

  //! Supercontroller status: True/False.
  bool scStatus;
  //! Name of the dynamic library containing the supercontroller implementation
  std::string scLibFile;
  //! Number of inputs and output to the supercontroller from/to each turbine
  int numScInputs, numScOutputs;

  //! Vector of turbine specific input data
  std::vector<turbineDataType>  globTurbineData;

  // Constructor 
  fastInputs() ;
 
  // Destructor
  ~fastInputs() {} ;

};


/**
 * A class to interface OpenFAST's fortran backend with a C++ driver program
 */    
class OpenFAST {

 private:

  //! MPI Communicator
  MPI_Comm mpiComm;
  //! The simulation will not run if dryRun is set to true. However, the simulation will read the input files, allocate turbines to processors and prepare to run the individual turbine instances. This flag is useful to test the setup of the simulation before running it.
  bool dryRun;        // If this is true, class will simply go through allocation and deallocation of turbine data
  //! Enable debug outputs if set to true
  bool debug;   // Write out extra information if this flags is turned on
  //! Number of turbines on this MPI rank
  int nTurbinesProc;
  //! Total number of turbines in the simulation
  int nTurbinesGlob;
  //! Start type of the simulation: 'init', 'trueRestart' or 'restartDriverInitFAST'.
  simStartType simStart;
  //! Is the time now zero: True/False
  bool timeZero;
  //! Time step for FAST. All turbines should have the same time step.
  double dtFAST;
  //! Number of OpenFAST time steps per unit time step of the driver program
  int nSubsteps_;  
  //! Is this the first pass through a time step
  bool firstPass_;  
  //! Max time of the simulation
  double tMax;
  //! Start time of the simulation
  double tStart;
  //! The current time step number 
  int nt_global;
  //! The starting time step number
  int ntStart;      
  //! Restart files will be written every so many time stneps
  int nEveryCheckPoint;

  std::vector<turbineDataType> globTurbineData; //All turbines
  std::vector<turbineDataType> turbineData;   // Only for turbines on the proc
  //! Velocity at force nodes - Store temporarily to interpolate to the velocity nodes
  std::vector<std::vector<std::vector<double> > > forceNodeVel; // Velocity at force nodes - Store temporarily to interpolate to the velocity nodes
  //! Position and velocity data at the velocity (aerodyn) nodes - (nTurbines, nTimesteps * nPoints * 6)
  std::vector<std::vector<double> > velNodeData; // Position and velocity data at the velocity (aerodyn) nodes - (nTurbines, nTimesteps * nPoints * 6)
  //! HDF-5 tag of file containing velocity (aerodyn) node data file
  hid_t velNodeDataFile; // HDF-5 tag of file containing velocity (aerodyn) node data file
  //! Array containing data at the velocity and force nodes    
  std::vector<std::vector<turbVelForceNodeDataType>> velForceNodeData;
  //! Array containing forces and deflections data for blade-resolved FSI simulations.
  std::vector<std::vector<turbBRfsiDataType>> brFSIData;

  //! Data structure to get forces and deflections from ExternalInflow module in OpenFAST
  std::vector<ExtInfw_InputType_t> extinfw_i_f_FAST; // Input from OpenFAST
  //! Data structure to send velocity information to ExternalInflow module in OpenFAST    
  std::vector<ExtInfw_OutputType_t> extinfw_o_t_FAST; // Output to OpenFAST

  //! Data structure to get deflections from ExternalLoads module in OpenFAST
  std::vector<ExtLdDX_InputType_t> extld_i_f_FAST; // Input from OpenFAST
  //! Data structure to send force information to ExternalLoads module in OpenFAST    
  std::vector<ExtLdDX_OutputType_t> extld_o_t_FAST; // Output to OpenFAST
    
  //! Data structure to get inputs to supercontroller from turbines in OpenFAST
  std::vector<SC_InputType_t> sc_i_f_FAST;
  //! Data structure to send outputs from supercontroller to turbines in OpenFAST    
  std::vector<SC_OutputType_t> sc_o_t_FAST;
    
  // Mapping of local turbine number to global turbine and processor number
  // Turbine Number is DIFFERENT from TurbID. Turbine Number simply runs from 0:n-1 locally and globally.
  //! Mapping global turbine number to processor number
  std::map<int, int> turbineMapGlobToProc;
  //! Mapping local to global turbine number
  std::map<int, int> turbineMapProcToGlob;
  //! Reverse Mapping global turbine number to local turbine number
  std::map<int, int> reverseTurbineMapProcToGlob;
  //! Set of processors containing atleast one turbine 
  std::set<int> turbineSetProcs;
  //! Same as the turbineSetProcs, but as an integer array
  std::vector<int> turbineProcs; 

  // Supercontroller stuff
  bool scStatus;
  std::string scLibFile;
  // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
  void *scLibHandle ; 
  typedef SuperController* create_sc_t(); 
  create_sc_t * create_SuperController;
  typedef void destroy_sc_t(SuperController *); 
  destroy_sc_t * destroy_SuperController;
  SuperController * sc;

  //! Number of inputs and output to the supercontroller from/to each turbine    
  int numScOutputs;  // # outputs from the supercontroller == # inputs to the controller == NumSC2Ctrl
  int numScInputs;   // # inputs to the supercontroller == # outputs from the controller == NumCtrl2SC
  std::vector<double> scOutputsGlob;  // # outputs from the supercontroller for all turbines
  std::vector<double> scInputsGlob;   // # inputs to the supercontroller for all turbines    

  // MPI related book keeping for all processors containing turbines
  //! Number of processors in a fastMPIGroup
  int fastMPIGroupSize;
  //! An MPI group created among all processors that simulate atleast one turbine
  MPI_Group fastMPIGroup;
  //! An MPI communicator for the MPI group created among all processors that simulate atleast one turbine
  MPI_Comm  fastMPIComm;
  //! MPI rank of processor on the fastMPIComm
  int fastMPIRank;

  //! Global MPI group
  MPI_Group worldMPIGroup;
  //! MPI rank of processor on global MPI Comm
  int worldMPIRank;

  //! Error status and Error message to communicate with OpenFAST
  int ErrStat;
  char ErrMsg[INTERFACE_STRING_LENGTH];  // make sure this is the same size as IntfStrLen in FAST_Library.f90
  static int AbortErrLev;

 public: 

  //! Constructor 
  OpenFAST() ;
  
  //! Destructor
  ~OpenFAST() {} ;

  //! Set inputs to OpenFAST through an object of the class fastInputs
  void setInputs(const fastInputs &);  

  //! Initialize the simulation - allocate memory for all data structures and initialize all turbines
  void init();
  //! Call FAST->solution0 for all turbines
  void solution0();
  //! Initialize velocity and force node data 
  void init_velForceNodeData();
  //! Set up before every OpenFAST time step
  void prework();
  //! Update states to next time step by calling FAST_AdvanceStates_T and CalcOutputs_And_SolveForInputs
  void update_states_driver_time_step();
  //! Copy the final predicted states from step t_global_next to actual states for that step
  void advance_to_next_driver_time_step();
  //! Set external inputs for OpenFAST modules by interpolating to substep
  void send_data_to_openfast(double ss_time);
  //! Set external inputs for OpenFAST modules at time step 't'
  void send_data_to_openfast(fast::timeStep t);
  //! Get ouput data from OpenFAST modules 
  void get_data_from_openfast(fast::timeStep t);
  //! Get reference positions of nodes from OpenFAST
  void get_refPositions(int iTurb);
  //! Extrapolate velocity and force node data to time step 'n+1' using data at 'n', 'n-1' and 'n-2'
  void predict_states();
  //! Advance all turbines by 1 OpenFAST timestep. Also write velocity data at the Aerodyn nodes for all turbines to 'velData.h5'
  void step();
  //! Step function to be used with sub-stepping fast between time steps of the driver program
  void step(double ss_time); 
  //! Advance all turbines by 1 OpenFAST timestep. Will not write any velocity data to file
  void stepNoWrite();
  //! Call FAST->end for all turbines
  void end();

  //! Read velocity and force node data at time step 'n', 'n-1' and 'n-2' to allow for a clean restart
  void readRestartFile(int iTurbLoc, int n_t_global);
  //! Write velocity and force node data at time step 'n', 'n-1' and 'n-2' to allow for a clean restart    
  void writeRestartFile(int iTurbLoc, int n_t_global);
    
  //! Compute the nacelle force
  void calc_nacelle_force(const float & u,
                          const float & v,
                          const float & w,
                          const float & cd,
                          const float & area,
                          const float & rho,
                          float & fx,
                          float & fy,
                          float & fz);
                        
  //! Open 'velData.h5' file to in preparation to write velocity data for all turbines
  hid_t openVelocityDataFile(bool createFile);
  //! Read velocity data at the Aerodyn nodes for all turbines from 'velData.h5'
  void readVelocityData(int nTimesteps);
  //! Write velocity data at the Aerodyn nodes for all turbines to 'velData.h5'  
  void writeVelocityData(hid_t h5file, int iTurb, int iTimestep, ExtInfw_InputType_t iData, ExtInfw_OutputType_t oData);
  //! Close 'velData.h5' file 
  herr_t closeVelocityDataFile(int nt_global, hid_t velDataFile);
  //! Backup the velocity data file
  void backupVelocityDataFile(int curTimeStep, hid_t & velDataFile);

  //! Allocate turbine number 'iTurbGlob' to the processor with global MPI rank 'procNo'
  void setTurbineProcNo(int iTurbGlob, int procNo) { turbineMapGlobToProc[iTurbGlob] = procNo; }
  //! Allocate all turbines to processors in a round-robin fashion
  void allocateTurbinesToProcsSimple();
  //! Get the approximate hub position for turbine number 'iTurbGlob'. This is the value specified in the input to OpenFAST. 
  void getApproxHubPos(std::vector<double> & currentCoords, int iTurbGlob);
  //! Get the exact hub position for turbine number 'iTurbGlob'. This is avaiable only after OpenFAST has been initialized for a given turbine.
  void getHubPos(std::vector<double> & currentCoords, int iTurbGlob, fast::timeStep t);
  //! Get a vector pointing downstream along the hub for turbine number 'iTurbGlob'
  void getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob, fast::timeStep t);

  //! Get the node type (HUB, BLADE, TOWER) of velocity node number 'iNode' for turbine number 'iTurbGlob'
  ActuatorNodeType getVelNodeType(int iTurbGlob, int iNode);
  //! Get the coordinates of velocity node number 'iNode' for turbine number 'iTurbGlob'
  void getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t);
  //! Set the velocity at velocity node 'iNode' for turbine number 'iTurbGlob'
  void setVelocity(std::vector<double> & velocity, int iNode, int iTurbGlob);
  //! Set the velocity at force node 'iNode' for turbine number 'iTurbGlob'
  void setVelocityForceNode(std::vector<double> & velocity, int iNode, int iTurbGlob);
  //! Map the velocity from the force nodes to the velocity nodes using linear interpolation along each blade and the tower.
  void interpolateVel_ForceToVelNodes();
  //! Get the node type (HUB, BLADE, TOWER) of force node number 'iNode' for turbine number 'iTurbGlob'
  ActuatorNodeType getForceNodeType(int iTurbGlob, int iNode);
  //! Get the coordinates of force node number 'iNode' for turbine number 'iTurbGlob'
  void getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t);
  //! Get the tensor orientation of force node number 'iNode' for turbine number 'iTurbGlob'
  void getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob, fast::timeStep t);
  //! Get the actuator force at force node 'iNode' for turbine number 'iTurbGlob'
  void getForce(std::vector<double> & force, int iNode, int iTurbGlob, fast::timeStep t);
  //! Get the chord at force node 'iNode' for turbine number 'iTurbGlob'
  double getChord(int iNode, int iTurbGlob);


  //! Get the blade reference positions array 'bldRefPos' of turbine number 'iTurbGlob' 
  void getBladeRefPositions(std::vector<double> & bldRefPos, int iTurbGlob);
  //! Get the blade deflections array 'bldDefl' of turbine number 'iTurbGlob' at time step 't'
    void getBladeDisplacements(std::vector<double> & bldDefl, std::vector<double> & bldVel, int iTurbGlob, fast::timeStep t);
  //! Get the tower reference positions array 'twrRefPos' of turbine number 'iTurbGlob' 
  void getTowerRefPositions(std::vector<double> & twrRefPos, int iTurbGlob);
  //! Get the tower deflections array 'twrDefl' of turbine number 'iTurbGlob' at time step 't'
  void getTowerDisplacements(std::vector<double> & twrDefl, std::vector<double> & twrVel, int iTurbGlob, fast::timeStep t);
  //! Get the hub reference position array 'hubRefPos' of turbine number 'iTurbGlob' 
  void getHubRefPosition(std::vector<double> & hubRefPos, int iTurbGlob);
  //! Get the hub deflections array 'hubDefl' of turbine number 'iTurbGlob' at time step 't'
  void getHubDisplacement(std::vector<double> & hubDefl, std::vector<double> & hubVel, int iTurbGlob, fast::timeStep t);
  //! Get the nacelle reference position array 'nacRefPos' of turbine number 'iTurbGlob' 
  void getNacelleRefPosition(std::vector<double> & nacRefPos, int iTurbGlob);
  //! Get the nacelle deflections array 'nacDefl' of turbine number 'iTurbGlob' at time step 't'
  void getNacelleDisplacement(std::vector<double> & nacDefl, std::vector<double> & nacVel, int iTurbGlob, fast::timeStep t);
   
  //! Set the blade forces array 'bldForce' for blade 'iBlade' of turbine number 'iTurbGlob' at time step 't'
  void setBladeForces(std::vector<double> & bldForce, int iTurbGlob, fast::timeStep t);
  //! Set the tower force array 'twrForce' of turbine number 'iTurbGlob' at time step 't'
  void setTowerForces(std::vector<double> & twrForce, int iTurbGlob, fast::timeStep t);
    

  //! Get all turbine parametric data
  void get_turbineParams(int iTurbGlob, turbineDataType & turbData);
  //! Get the starting time step of the simulation
  int get_ntStart() { return ntStart; }
  //! Return a boolean flag whether the simulation is dryRun
  bool isDryRun() { return dryRun; }
  //! Return a boolean flag whether the simulation is debug
  bool isDebug() { return debug; }
  //! Get an enum of type 'simStartType' indicating the start type of the simulation
  simStartType get_simStartType() { return simStart; }
  //! Is the simulation time zero right now?
  bool isTimeZero() { return timeZero; }
  //! Get the global MPI rank of the processor containing turbine number 'iTurbGlob'
  int get_procNo(int iTurbGlob) { return turbineMapGlobToProc[iTurbGlob] ; }
  //! Get the local turbine number of the turbine number 'iTurbGlob'
  int get_localTurbNo(int iTurbGlob) { return reverseTurbineMapProcToGlob[iTurbGlob]; }
  //! Get the total number of turbines in the simulation
  int get_nTurbinesGlob() { return nTurbinesGlob; } 

  //! Get the nacelle drag coefficient of turbine number 'iTurbGlob'
  float get_nacelleArea(int iTurbGlob) { return get_nacelleAreaLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the nacelle area of turbine number 'iTurbGlob'
  float get_nacelleCd(int iTurbGlob) { return get_nacelleCdLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the air density around turbine number 'iTurbGlob'
  float get_airDensity(int iTurbGlob) { return get_airDensityLoc(get_localTurbNo(iTurbGlob)); }
    
  //! Get the number of blades in turbine number 'iTurbGlob'
  int get_numBlades(int iTurbGlob) { return get_numBladesLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Aerodyn/velocity nodes on each blade in turbine number 'iTurbGlob'
  int get_numVelPtsBlade(int iTurbGlob) { return get_numVelPtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Aerodyn/velocity nodes on the tower in turbine number 'iTurbGlob'
  int get_numVelPtsTwr(int iTurbGlob) { return get_numVelPtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the total number of Aerodyn/velocity nodes in turbine number 'iTurbGlob'
  int get_numVelPts(int iTurbGlob) { return get_numVelPtsLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Actuator/force nodes on each blade in turbine number 'iTurbGlob'
  int get_numForcePtsBlade(int iTurbGlob) { return get_numForcePtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Actuator/force nodes on the tower in turbine number 'iTurbGlob'
  int get_numForcePtsTwr(int iTurbGlob) { return get_numForcePtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the total number of Actuator/force nodes in turbine number 'iTurbGlob'
  int get_numForcePts(int iTurbGlob) { return get_numForcePtsLoc(get_localTurbNo(iTurbGlob)); }

  //! Compute the torque and thrust for turbine number 'iTurbGlob'
  void computeTorqueThrust(int iTurGlob, std::vector<double> &  torque, std::vector<double> &  thrust);

  //! An example function to set velocities at the Aerodyn nodes using a power law wind profile using an exponent of 0.2 and a reference wind speed of 10 m/s at 90 meters.
  void setExpLawWindSpeed() ; 

  //! An example function to set a uniform X force at all blade nodes
  void setUniformXBladeForces();

 private:

  //! Check whether the error status is ok. If not quit gracefully by printing the error message
  void checkError(const int ErrStat, const char * ErrMsg);
  //! Check whether a file with name "name" exists
  inline bool checkFileExists(const std::string& name);

  //! Allocate memory for data structures for all turbines on this processor
  void allocateMemory();
  //! Allocate more memory for each turbine after intialization/restart    
  void allocateMemory2(int iTurbLoc);  

  //! Get the nacelle drag coefficient of local turbine number 'iTurbLoc'
  float get_nacelleCdLoc(int iTurbLoc) { return turbineData[iTurbLoc].nacelle_cd; }
  //! Get the nacelle area of local turbine number 'iTurbLoc'    
  float get_nacelleAreaLoc(int iTurbLoc) { return turbineData[iTurbLoc].nacelle_area; }
  //! Get the air density around local turbine number 'iTurbLoc'    
  float get_airDensityLoc(int iTurbLoc) { return turbineData[iTurbLoc].air_density; }
    
  //! Get the number of blades in local turbine number 'iTurbLoc'
  int get_numBladesLoc(int iTurbLoc) { return turbineData[iTurbLoc].numBlades; }
  //! Get the number of Aerodyn/velocity nodes on each blade in local turbine number 'iTurbLoc'
  int get_numVelPtsBladeLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPtsBlade; }
  //! Get the number of Aerodyn/velocity nodes on the tower in local turbine number 'iTurbLoc'
  int get_numVelPtsTwrLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPtsTwr; }
  //! Get the total number of Aerodyn/velocity nodes in local turbine number 'iTurbLoc'
  int get_numVelPtsLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPts; }
  //! Get the number of Actuator/force nodes on each blade in local turbine number 'iTurbLoc'
  int get_numForcePtsBladeLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePtsBlade; }
  //! Get the number of Actuator/force nodes on the tower in local turbine number 'iTurbLoc'
  int get_numForcePtsTwrLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePtsTwr; }
  //! Get the total number of Actuator/force nodes in local turbine number 'iTurbLoc'
  int get_numForcePtsLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePts; }

  void loadSuperController(const fastInputs & fi);
  void fillScInputsGlob() ;
  void fillScOutputsLoc() ;

  //! Apply the velocity data at the Aerodyn nodes in 'velData' to turbine number 'iTurb' at time step 'iPrestart' through the data structure 'cDriver_Output_to_FAST'
  void applyVelocityData(int iPrestart, int iTurb, ExtInfw_OutputType_t o_t_FAST, std::vector<double> & velData) ;

};

}

#endif