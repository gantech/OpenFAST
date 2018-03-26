#ifndef OpenFAST_h
#define OpenFAST_h
#include "FAST_Library.h"
#include "sys/stat.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include "dlfcn.h"
#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "mpi.h"
#include "SC.h"


namespace fast {

    struct turbineDataType {
        int TurbID;
        std::string FASTInputFileName;
        std::string FASTRestartFileName;
        std::vector<float> TurbineBasePos;
        std::vector<double> TurbineHubPos;
        int numBlades;
        int numVelPtsBlade;
        int numVelPtsTwr;
        int numVelPts;
        int numForcePtsBlade;
        int numForcePtsTwr;
        int numForcePts;
        int inflowType;
    };
    
enum ActuatorNodeType {
  HUB = 0,
  BLADE = 1,
  TOWER = 2,
  ActuatorNodeType_END
};

enum simStartType {
  init = 0,
  trueRestart = 1,
  restartDriverInitFAST = 2,
  simStartType_END    
};

struct turbVelForceNodeDataType {
    std::vector<double> x_vel;        // Blade location at velocity nodes
    std::vector<double> xdot_vel;     // Blade velocity at velocity nodes
    std::vector<double> vel_vel;      // Sampled velocity at velocity nodes
    std::vector<double> x_force;      // Blade location at force nodes
    std::vector<double> xdot_force;   // Blade velocity at force nodes
    std::vector<double> orient_force; // Blade orientation at force nodes
    std::vector<double> vel_force;    // Sampled velocity at force nodes
    std::vector<double> force;        // Actuator force at force nodes
    double x_vel_resid;
    double xdot_vel_resid;
    double vel_vel_resid;
    double x_force_resid;
    double xdot_force_resid;
    double orient_force_resid;
    double vel_force_resid;
    double force_resid;
};

enum timeStep {
    nm1 = 0,
    n = 1,
    np1 = 2,
    timeStep_END
};

    
    
class fastInputs {

 public:

  MPI_Comm comm;
  int nTurbinesGlob;  
  bool dryRun;
  bool debug; 
  double tStart;
  simStartType simStart;
  int nEveryCheckPoint;  
  double tMax;
  double dtFAST;
  int nSubsteps;

  bool scStatus;
  std::string scLibFile;
  int numScInputs, numScOutputs;

  std::vector<turbineDataType>  globTurbineData;

  // Constructor 
  fastInputs() ;
 
  // Destructor
  ~fastInputs() {} ;

};


class OpenFAST {

 private:

  MPI_Comm mpiComm;
  bool dryRun;        // If this is true, class will simply go through allocation and deallocation of turbine data
  bool debug;   // Write out extra information if this flags is turned on
  int nTurbinesProc;
  int nTurbinesGlob;
  simStartType simStart;
  bool timeZero;
  double dtFAST;
  int nSubsteps_;  // Number of OpenFAST time steps per unit time step of the driver program
  double tMax;
  double tStart;
  int nt_global;
  bool firstPass_;  //Is this the first pass through a time step 
  int ntStart;      // The time step to start the FAST simulation
  int nEveryCheckPoint;    // Check point files will be written every 'nEveryCheckPoint' time steps

  std::vector<turbineDataType> globTurbineData; //All turbines
  std::vector<turbineDataType> turbineData;   // Only for turbines on the proc
  std::vector<std::vector<std::vector<double> > > forceNodeVel; // Velocity at force nodes - Store temporarily to interpolate to the velocity nodes
  std::vector<std::vector<double> > velNodeData; // Position and velocity data at the velocity (aerodyn) nodes - (nTurbines, nTimesteps * nPoints * 6)
  hid_t velNodeDataFile; // HDF-5 tag of file containing velocity (aerodyn) node data file
  std::vector<std::vector<turbVelForceNodeDataType>> velForceNodeData; // Array containing data at the velocity and force nodes

  std::vector<OpFM_InputType_t> i_f_FAST; // Input from OpenFAST
  std::vector<OpFM_OutputType_t> o_t_FAST; // Output to OpenFAST

  std::vector<SC_InputType_t> sc_i_f_FAST;
  std::vector<SC_OutputType_t> sc_o_t_FAST;

  // Mapping of local turbine number to global turbine and processor number
  // Turbine Number is DIFFERENT from TurbID. Turbine Number simply runs from 0:n-1 locally and globally.
  std::map<int, int> turbineMapGlobToProc; // Mapping global turbine number to processor number
  std::map<int, int> turbineMapProcToGlob; // Mapping local to global turbine number
  std::map<int, int> reverseTurbineMapProcToGlob; // Reverse Mapping global turbine number to local turbine number
  std::set<int> turbineSetProcs; // Set of processors containing atleast one turbine 
  std::vector<int> turbineProcs; // Same as the turbineSetProcs, but as an integer array

  //Supercontroller stuff
  bool scStatus;
  std::string scLibFile;
  // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
  void *scLibHandle ; 
  typedef SuperController* create_sc_t(); 
  create_sc_t * create_SuperController;
  typedef void destroy_sc_t(SuperController *); 
  destroy_sc_t * destroy_SuperController;
  SuperController * sc;

  int numScOutputs;  // # outputs from the supercontroller == # inputs to the controller == NumSC2Ctrl
  int numScInputs;   // # inputs to the supercontroller == # outputs from the controller == NumCtrl2SC
  std::vector<double> scOutputsGlob;  // # outputs from the supercontroller for all turbines
  std::vector<double> scInputsGlob;   // # inputs to the supercontroller for all turbines    

  // MPI related book keeping for all processors containing turbines
  int fastMPIGroupSize;
  MPI_Group fastMPIGroup;
  MPI_Comm  fastMPIComm;
  int fastMPIRank;

  MPI_Group worldMPIGroup;
  int worldMPIRank;

  int ErrStat;
  char ErrMsg[INTERFACE_STRING_LENGTH];  // make sure this is the same size as IntfStrLen in FAST_Library.f90

 public: 

  // Constructor 
  OpenFAST() ;
  
  // Destructor
  ~OpenFAST() {} ;

  void setInputs(const fastInputs &);  

  void init();
  void solution0();
  void prework();
  void update_states_driver_time_step();
  void advance_to_next_driver_time_step();
  void send_data_to_openfast(double ss_time);
  void send_data_to_openfast(fast::timeStep t);
  void get_data_from_openfast(fast::timeStep t);
  void predict_states();
  void step();
  void step(double ss_time); // Step function to be used with sub-stepping fast between time steps of the driver program
  void stepNoWrite();
  void end();

  void readRestartFile(int iTurbLoc, int n_t_global);
  void writeRestartFile(int iTurbLoc, int n_t_global);
    
  hid_t openVelocityDataFile(bool createFile);
  void readVelocityData(int nTimesteps);
  void writeVelocityData(hid_t h5file, int iTurb, int iTimestep, OpFM_InputType_t iData, OpFM_OutputType_t oData);
  herr_t closeVelocityDataFile(int nt_global, hid_t velDataFile);

  void setTurbineProcNo(int iTurbGlob, int procNo) { turbineMapGlobToProc[iTurbGlob] = procNo; }
  void allocateTurbinesToProcsSimple();
  void getApproxHubPos(std::vector<double> & currentCoords, int iTurbGlob);
  void getHubPos(std::vector<double> & currentCoords, int iTurbGlob, fast::timeStep t);
  void getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob, fast::timeStep t);

  ActuatorNodeType getVelNodeType(int iTurbGlob, int iNode);
  void getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t);
  void setVelocity(std::vector<double> & velocity, int iNode, int iTurbGlob);
  void setVelocityForceNode(std::vector<double> & velocity, int iNode, int iTurbGlob);
  void interpolateVel_ForceToVelNodes();
  ActuatorNodeType getForceNodeType(int iTurbGlob, int iNode);
  void getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t);
  void getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob, fast::timeStep t);
  void getForce(std::vector<double> & force, int iNode, int iTurbGlob, fast::timeStep t);
  double getChord(int iNode, int iTurbGlob);

  int get_ntStart() { return ntStart; }
  bool isDryRun() { return dryRun; }
  bool isDebug() { return debug; }
  simStartType get_simStartType() { return simStart; }
  bool isTimeZero() { return timeZero; }
  int get_procNo(int iTurbGlob) { return turbineMapGlobToProc[iTurbGlob] ; } // Get processor number of a turbine with global id 'iTurbGlob'
  int get_localTurbNo(int iTurbGlob) { return reverseTurbineMapProcToGlob[iTurbGlob]; }
  int get_nTurbinesGlob() { return nTurbinesGlob; } 

  int get_numBlades(int iTurbGlob) { return get_numBladesLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPtsBlade(int iTurbGlob) { return get_numVelPtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPtsTwr(int iTurbGlob) { return get_numVelPtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPts(int iTurbGlob) { return get_numVelPtsLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePtsBlade(int iTurbGlob) { return get_numForcePtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePtsTwr(int iTurbGlob) { return get_numForcePtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePts(int iTurbGlob) { return get_numForcePtsLoc(get_localTurbNo(iTurbGlob)); }

  void computeTorqueThrust(int iTurGlob, std::vector<double> &  torque, std::vector<double> &  thrust);

 private:

  void checkError(const int ErrStat, const char * ErrMsg);
  inline bool checkFileExists(const std::string& name);

  void allocateMemory();
  void allocateMemory2(int iTurbLoc);  //Allocate more memory for each turbine after intialization/restart
  
  int get_numBladesLoc(int iTurbLoc) { return turbineData[iTurbLoc].numBlades; }
  int get_numVelPtsBladeLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPtsBlade; }
  int get_numVelPtsTwrLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPtsTwr; }
  int get_numVelPtsLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPts; }
  int get_numForcePtsBladeLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePtsBlade; }
  int get_numForcePtsTwrLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePtsTwr; }
  int get_numForcePtsLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePts; }

  void loadSuperController(const fastInputs & fi);
  void fillScInputsGlob() ;
  void fillScOutputsLoc() ;

  void setExpLawWindSpeed(int iTurbLoc) ; // An example to set velocities at the Aerodyn nodes
  void applyVelocityData(int iPrestart, int iTurb, OpFM_OutputType_t o_t_FAST, std::vector<double> & velData) ;

};

}

#endif
