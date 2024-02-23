# ORCh
ORCh (Optimised and Reduced Chemistry) is a fully automated method to reduce detailed chemical schemes.

ORCh is a preprocessing tool designed to automatically generate and optimise reduced chemistry for specific CFD conditions, including fuel spray, heat losses, etc.

ORCh combines stochastic methods with graph-analysis together with genetic algorithms. Chemical kinetics are probed in laminar flames and turbulent micro-mixing (pairwise interaction or euclidean minimum spanning tree) canonical problems.

ORCh is coupled with HYPE (HYbrid stochastic sectional method for solving Population balance Equations) and machine learning to deal with carbon nanoparticles formation (nanotube and soot) and non-inertial solid particles (crystallisation).

The development of ORCh is reported in two archival papers (make sure to cite those papers when using this program):

- N. Jaouen, L. Vervisch, P. Domingo, G. Ribert (2017) Automatic reduction and optimisation of chemistry for turbulent combustion modeling: Impact of the canonical problem, Combust. Flame, 175: 60-79.
- N. Jaouen, L. Vervisch, P. Domingo (2017) Auto-thermal reforming (ATR) of natural gas: An automated derivation of optimised reduced chemical schemes, Proc. Combust. Inst., 36(3): 3321-3330.



ORCh has been applied to various reactive flow problems:
- K. Bioche, G. Ribert, L. Vervisch (2019) Simulating upstream flame propagation in a narrow channel after wall preheating: Flame analysis and chemistry reduction strategy, Combust. Flame. 200: 219-231.
- C. Locci, L. Vervisch, B. Farcy, P. Domingo, N. Perret (2018) Selective Non-Catalytic Reduction (SNCR) of nitrogen oxide emissions: A perspective from numerical modeling, Flow Turbulence and Combust. 100(2): 301-340.
- K. Bioche, L. Vervisch, G. Ribert (2018) Premixed flame-wall interaction in a narrow channel: Impact of wall thermal conductivity and heat losses, J. Fluid Mech. 856: 5-35
- A. Bouaniche, N. Jaouen, P. Domingo, L. Vervisch (2019) Vitiated high Karlovitz n-decane/air turbulent flames: Scaling laws and micro-mixing modeling analysis, Flow Turbulence and Combust. 102(1): 235–252.
- K. Wan, L. Vervisch, C. Jianga, P. Domingo, Z. Gao, J. Xia, Z. Wang (2020) Development of reduced and optimized reaction mechanism for potassium emissions during biomass combustion based on genetic algorithms, Energy 211: 118565.

ORCh can include the training of ANN-CNN to reduce chemistry and pre-integrate stiff chemical systems or to solve for soot:

- K. Wan, C. Barnaud, L. Vervisch, P. Domingo (2020) Chemistry reduction using machine learning trained from non-premixed micro-mixing modeling: Application to DNS of a syngas turbulent oxy-flame with side-wall effects, Combust. Flame 220: 119-129.
- A. Seltz, P. Domingo, L. Vervisch (2021) Solving the population balance equation for non-inertial particles dynamics using PDF and neural networks: Application to a sooting flame, Phys. Fluids. 33, 013311.
- H.-T. Nguyen, C. Barnaud, P. Domingo, P.-D. Nguyen, L. Vervisch (2023) Large-Eddy Simulation of flameless combustion with neural-network driven chemistry, Application Energy Combust. Sci. 14:100126.

Most recent and on-going applications:
- Furnaces with Urea DeNOx
- High-pressure partial oxydation of methane (hydrogen production)
- Combustion in narrow channels with wall heat losses (micro- and meso-scale combustion, fire propagation between battery cells)
- Heavy hydrocarbon liquid fuel (kerosene surrogate or pyrolyse)
- Chlorofluorocarbon chemistry
- Soot precursors
- Coal combustion
- Combustion of blast furnace gases



# Pre installation for Cantera 2.4 and ORCh
ORCh uses the 2.4 version of Cantera https://cantera.org/blog/cantera-240-released.html, it needs:
- C++ Boost librairies (Version 1.68, the newer don't work with Cantera) to run properly. 
- SCONS compiler

- SCONS compiler
  - https://sourceforge.net/projects/scons/files/scons/3.0.3/scons-3.0.3.tar.gz/download
  - cd scons-3.0.3
  - python setup.py install
  - (this scons version requires python >= 3.5 and <= 3.10)

- BOOST LIBRARIES
  - https://www.boost.org/users/history/version_1_68_0.html
  - cd boost_1_68_0 
  - sh bootstrap.sh
  - ./b2  install 
  - ... patience ... the boost installation may take some time.

- Eigen
  - visit https://eigen.tuxfamily.org/index.php?title=Main_Page
  - Download the 3.2.10 tar
  - untar it

- Sundials version 2.6.2
  - visit https://computing.llnl.gov/projects/sundials/sundials-software
  - download the 2.6.2 archive
  - create empty directory : mkdir SUNDIALS
  - move downloaded sundials directory to SUNDIALS : mv sundials-2.6.2 path/to/SUNDIALS/
  - cd SUNDIALS
  - mkdir builddir
  - mkdir instdir
  - cd builddir
  - cmake -DCMAKE_INSTALL_PREFIX=/home/example/path/to/SUNDIALS/instdir -DEXAMPLE_INSTALL_PATH=/home/example/path/to/SUNDIALS/instdir/examples ../sundials-2.6.2
  - make install

- Opencv 3.4.12
  - visit https://opencv.org/releases/
  - download the 3.4.12 release
  - create empty directory : mkdir OPENCV
  - move downloaded opencv directory to OPENCV : mv opencv-3-4-12 path/to/OPENCV/
  - mkdir builddir
  - cd builddir
  - cmake -DCMAKE_INSTALL_PREFIX='/home/path/to/OPENCV/' ../opencv-3-4-12
  - make (you can use make -j8 to speed up the compilation)
  - make install


- Tensorflow C api 1.15
  - download 1.15 version cpu only : https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-linux-x86_64-1.15.0.tar.gz
  - untar it


# MPI and hdf5 librairies will also be needed
if it is not already installed on your computer
- visit https://www.open-mpi.org/software/ompi/v1.8/
  - download openmpi 1.8.8
  - cd openmpi-1.8.8
  - ./configure --prefix=$HOME/opt/openmpi (or any other location you want it installed)
  - make all
  - make install
- visit https://www.hdfgroup.org/packages/hdf5-1107-source/
  - download hdf5 1.10.7
  - cd hdf5-1.10.7
  - ./configure --prefix=$HOME/opt/hdf5 (or any other location you want it installed)
  - make all
  - make install



# Installation
The following steps must be completed before you can use ORCh:
- In the directory workdir/ORCh_IFPEN-ORCh
  - fill in the path to the directories of the newly installed librairies in set_orch_env.sh
  - bash set_orch_env.sh
  - source set_orch_env.sh
- In the directory workdir/ORCh_IFPEN-ORCh/Cantera
  - build with SCONS : "scons build"
  - If the build fail, for any reason, use "scons clean" before doing again the previous step, after coping with problems indicated in the output.
  - Then : "scons install"



If you downloaded openmpi and hdf5, you need to add these lines in the set_orch_env.sh file, putting the right paths (then you need to re do a "source set_orch_env.sh") :

- <code>export MPI_HOME="/home/username/opt/openmpi"</code>
- <code>export MPI_INC="$MPI_HOME/include"</code>
- <code>export MPI_LIB="$MPI_HOME/lib"</code>
- <code>export MPI_BIN="$MPI_HOME/bin"</code>
- <code>export PATH="$MPI_BIN:$PATH"</code>
- <code>export LIBRARY_PATH="$MPI_LIB:$LIBRARY_PATH"</code>
- <code>export LD_LIBRARY_PATH="$MPI_LIB:$LD_LIBRARY_PATH"</code>

- <code>export HDF5_HOME="/home/username/opt/hdf5"</code>
- <code>export HDF5_INC="$HDF5_HOME/include"</code>
- <code>export HDF5_LIB="$HDF5_HOME/lib"</code>
- <code>export LIBRARY_PATH="$HDF5_LIB:$LIBRARY_PATH"</code>
- <code>export LD_LIBRARY_PATH="$HDF5_LIB:$LD_LIBRARY_PATH"</code>

In the directory workdir/ORCh_IFPEN-ORCh/ORCh:
  - make clean
  - make

If any library problem occurs during make:
- does no find library ifcore -lifcore
   - if libifcore is not on your computer, you can install intel fortran compiler "Stand-alone version" from https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.2hxff2
   - bash l_fortran-compiler_p_2024.0.1.31.sh
   - add to set_orch_env.sh and re do "source set_orch_env.sh" in workdir/ORCh_IFPEN-ORCh
      -  <code>export INTEL_LIB="/home/username/intel/oneapi/compiler/2024.0/lib"</code>
      -  <code>export LIBRARY_PATH="${LIBRARY_PATH}:${INTEL_LIB}"</code>
      -  <code>export LD_LIBRARY_PATH="${INTEL_LIB}:${LD_LIBRARY_PATH}"</code>
- undefined reference to functions in libuser.a
   - go to workdir/ORCh_IFPEN-ORCh/Cantera/user
   - make clean
   - make
- c++ errors
   - try to add the compilation option -std=c++11 in all Makefile files from ORCh_IFPEN-ORCh/ORCh/cantera, drgep, lumping, main, exe, optimisation and tools


# Running ORCh
All the information necessary to run the subsequent steps of the reduction process are to be entered in the file: "conditions.cpp"

<code> new_mixing = 'False' or 'True' </code> 
- if 'False': the stochastic mixing will follow the same pattern during trajectories (i.e. the same time history of random numbers are used along the composition space trajectories);
- if 'True': the time history of the stochastic mixing differs from the previous run;
This parameter must be 'True' for the first simulation of the composition space trajectories, then one may decide to perform the reduction and optimisation with the same, or a new, time history of random numbers controlling the stochastic mixing. So far, ORCh produced good results using 'False', but for the first run where it must be set to 'True'.

Chemical scheme files:
- The reference detailed scheme is entered here:
 <code> mech_ref = "mechanisms/gri12.xml"; </code> 

- Intermediate scheme from which the reduction should be pursued (after completing one of the steps), it must be updated at every stage:
 <code>  mech = "./outputs/mechanisms/drgepSpecies17.xml"; </code> 

- If first run, use the same as the reference detailed mechanism:
 <code>    mech = "mechanisms/gri12.xml";</code> 
- Always provide a cantera id for the chemical scheme (you can choose any, but it must be there and agrees between files, i.e. do not change it during the whole optimisation process):
 <code>  mech_desc = "gri12";</code> 

For instance, you may find it in the "gri12.xml" file
- <code> phase dim="3" id="gri12" </code>

Targeted species: The species that you wish to reproduce with accuracy in the reduced scheme should be given in 'conditions.cpp' following the format:
- <code> listTargets.push_back("H2");</code>
- <code> listTargets.push_back("O2");</code>
- <code> listTargets.push_back("H2O");</code>
- <code> listTargets.push_back("CH4");</code>
- <code> listTargets.push_back("CO");</code>
- <code> listTargets.push_back("CO2"); </code>

These target species will have the highest DRGEP coefficient (=1) and the hight weight during optimisation
(note that for the DRGEP steps the species into reactant are automatically preserved, i.e. their DRGEP coefficient is also set to 1 so they can't be suppressed).

Inlet of composition space trajectories:
- You can prescribe as many inlets as you wish. For every inlet, you need to provide all these parameters: Temperature, Pressure, Mass flow rate (note that the number of particles per inlet is 1000 times the mass flow inlet in kg/s), and then directly the composition, either in mole fractions (first set, see below), or mass fractions if you prefer those, the code will understand anyway. 

- The last inlet must be composed of burnt gases, in which you need to give in addition: the mixing time, the time step and the number of iterations. 

- Of course, it will then be needed to check that, for these parameters, the solution is well captured and that you have reached the chemical equilibrium condition.

 <code> 
  //-- Temp, Pressure, number of particles, "moles fractions" or "mass fractions" leave empty the other one --//
  listInlets.push_back(new MultipleInlet(320, 1E+05, NbParticles1,
                        "O2:0.15, N2:0.52, H2O:0.0029, CH4:0.33", ""));
   //Inlet 2 //
   listInlets.push_back(new MultipleInlet(1350, 1E+05, NbParticles2,
                        "O2:0.12, N2:0.73, H2O:0.15, CH4:0.0003", ""));
   //-- The very last inlet must be burnt gases --//
   //-- Same as other inlet + mixing time tau_t, delta_t, number of iterations --//
   //BurnedGases
   double tau_t = 2e-03;
   double delta_t = 5e-05;
   bool BurnedGases = true;
   listInlets.push_back(new Characteristics_MultipleInlet(2290, 1E+05, NbParticlesBurnedGases,
                        "N2:0.682, H2O:0.236, CO2:0.0539, CO:0.00823, H2:0.00647, O2:0.00516, OH:0.00431", "",
                        tau_t, delta_t, 400, BurnedGases)); </code>

After any modification in 'conditions.cpp', the code must be re-compile, therefore this must be done between every step and sub-steps of the method
- <code> make clean </code>
-<code> make </code>

'''DRGEP''' is the first stage to run: it includes the computation of the trajectories with detailed chemistry and the DRGEP analysis to reduce the number of species and reactions:
- '''First sub-step''' - select in 'conditions.cpp': 
 <code> step = "DRGEP_Species"; </code>
- This first step provides a series of chemical schemes whose number of species ranges between the initial number of species of the mechanism and 10 species (this limit can be changed in main/mainDRGEPSpecies.cpp).
- Here schemes will be built with a total number of species ranging between 20 and 15. Notice that the control of the threshold, T<sub>S</sub>, of Jaouen et al. (see references above) is automatically computed from  'nbSpeciesToKeep' and 'nbSpeciesToKeep>' in 'ORCh/main.cpp'.

Run
-  <code> ./main </code>

Once done, you need to go into the 'outputs' directory, there you can examine the trajectories (Premixed or Stochastic) of the species you previously set in conditions.cpp:
-  <code> speciesToPlot.push_back("H2O") </code>

To find the positions of species in the files, you need to go into the reference and reduced mechanisms, which are stored in '../mechanisms/' and './outputs/mechanisms'

Select the best mechanism you like and put his name in 'conditions.cpp':
- <code>  mech = "./outputs/mechanisms/yourbestmechname.xml" </code>
- 'yourbestmechname.xml' is found in './outputs/mechanisms' and it looks like './outputs/mechanisms/drgepSpecies17.xml', here 17 means that this is a mechanism with 17 species.

Then, it is time for further reducing the number of reactions. Make sure that in 'conditions.cpp':
 - <code> new_mixing = 'False' </code> 
'''Second sub-step''' - select in 'conditions.cpp': 
 - <code> step = "DRGEP_Reactions"; </code>

This second step provides a series of chemical schemes whose number of reactions ranges between the initial reaction number and 10 reactions (limit that can be changed in main/mainDRGEPReactions.cpp).

Once done, you need to go into the 'outputs' directory, there you can examine the trajectories:
- To find the positions of species in the files, you need to go into the reference and reduced mechanisms, which are stored in '../mechanisms/' and './outputs/mechanisms'
- Select the best mechanism you like and put his name in 'conditions.cpp':
 <code>  mech = "./outputs/mechanisms/yourbestmechname.xml"; </code>
- 'yourbestmechname.xml' is found in './outputs/mechanisms' and it looks like './outputs/mechanisms/drgepReactions36.xml', here 36 means that this is a mechanism with 36 reactions.
- '''Automated quasi-steady state''' is the second stage to run: it includes the determination of species in quasi-steady state and the building of the needed relations:

'''First sub-step''' - select in 'conditions.cpp' 
-  <code> step = "computeQSSCriteria"; </code>

The code will rerun trajectories from the specified mechanism and compute the QSS criteria, to provide in an online output an evaluation of the QSS criterion for every species, i.e. the ratio between the integral over the trajectory of the net variation rate of a species and the integral of the maximum between the production rate and the consumption rate of that given species.

In addition, the output provides the list of species which are related in a non-linear manner, to decide which species can be put in QSS.
'''Second sub-step''' - select in 'conditions.cpp' 
 - <code> step = "getQSSfile"; </code>

The user choose a set of QSS scenarios he wants to try, this information should be put in 'conditions.cpp' (information which can be left there or not for all runs, but which is mandatory for this step):
- <code> string array1[1] {"CH3O"};
   vector<string> vec(array1, array1 + sizeof(array1) /sizeof(array1[0]));
   string array2[3] {"CH3O", "HCO", "HO2"};  
  vector<string> vec2(array2, array2 + sizeof(array2) /sizeof(array2[0])); 
  listQSSscenarios.push_back(new QSSscenario(vec2));</code>

 For every QSS scenarios you wish to examine, the code will do two things for you:
- Generate the analytical relations between the species and write them ready for their use in a CFD software in the file 'analytic_schemes/RefQSSAnalysis0/mech_QSS.h', for the conditions you entered.
- Rerun the trajectories for these conditions, so that you can have a look and check which QSS hypotheses are best. The results for visualisation are in 'analytic_schemes/RefQSSAnalysis0', 'analytic_schemes/RefQSSAnalysis1', etc.

If you need the QSS relations in fortran, after running as above (step = "getQSSfile"), you can run with the option:
- <code> step = "getQSSfileFORTRAN"; </code>

Notice that the trajectories will not be computed and that only the first QSS scenario is applied.
- '''Optimisation of the rates''' is the last stage to run: 
'''Preliminary step''' - select in 'conditions.cpp' 
 <code> step = "Optimisation"; </code>

The optimisation will proceed with the QSS conditions with the first QSS scenario given in 'conditions.cpp'. If one scenario exists, a QSS-step (above) must have been completed before. If none, the optimisation proceed for the full scheme without QSS.

The corresponding mechanism mut be given in 'conditions.cpp' (if you have run the previous run with various scenario, make sure to put here the mechanism corresponding to the one you prefer)
-  <code>  mech = "./outputs/mechanisms/yourbestmechname.xml"; </code>
- 'yourbestmechname.xml' is found in './outputs/mechanisms' and it looks like './outputs/mechanisms/drgepReactions36.xml', here 36 means that this is a mechanism with 36 reactions.

The parameters of the Genetic Algorithm must be given:
-  <code>  int PopSize = 100;</code>
- <code>  int MaxAllowableGenerations = 150;</code>
- <code>  int NbElitism = 1;</code>
-  <code> double CrossoverRate = 0.75;</code>
- <code>  double MutationRate = 0.02;</code>
- <code>  double AllowedVariation_A = 0.03;</code>
- <code>  double AllowedVariation_b = 0.03;</code>
- <code>  double AllowedVariation_E = 0.03; </code>
- 'PopSize' is the number of chromosomes >= 100 is advised
- 'MaxAllowableGenerations' number of generation to be examined >= 150 is advised
- 'NbElitism' number of solution copied from previous generation to secure elitism, 1 is best choice
- 'CrossoverRate' % of crossover, 0.75 is best choice
- 'MutationRate' % of mutation, 0.02 is best choice
- 'AllowedVariation_A' % of variation for pre-exponential constant 
- 'AllowedVariation_b' % of variation for temperature coefficient
- 'AllowedVariation_E' % of variation for activation energy

This last step may be run in parallel, the number of processors used should be chosen so that PopSize/n_proc is an integer:
- <code> mpirun -np n_proc main </code>
- During the run, temporary directories are created in 'tool/analytic_schemes/Ref#nbproc' to store all information concerning one chromosome (i.e. a trajectory obtained for one set of chemical parameters), with possibility to check corresponding trajectories.
- In 'analytic_schemes/PLOTS' .eps files are automatically generated to visualise these trajectories, for each chromosomes. The evolution of the best scheme fitness at each generation is also computed in the "fitnessEvolution.eps" file.
- These .eps file read 'GEN#nb_POP0_FIT#value.eps', where #nb is the  generation number and #value is the fitness value. Looking at the latest generation gives an idea of the quality of the parameters.
- At each iteration (generation) in 'analytic_schemes/Ref', the best chemical scheme is in 'chemistry_for_restart.xml'.
- In case of QSS the relations are in: 'best_mech_QSS.h'

The procedure to run the GA is as follow:
- Run for a given number of generation
- Check 'GEN#nb_POP0_FIT#value.eps' to evaluate the quality of the result. The evolution of the best fitness also indicates if the optimisation reached an optimum.
- If more generation are needed (i.e. continue the optimisation):
--  Move the 'analytic_schemes/Ref/chemistry_for_restart.xml' into 'outputs/mechanism' and modify the name of 'mech=' in 'conditions.cpp'
-- Decrease the 'Allowed variation' to narrow the research parameters area
-- Relaunch

This can be repeated up to a solution that suits you.
Then, you may rerun the trajectory only, with the optimised mechanism, in order to plot all the variables, in 'conditions.cpp':
- '''Computation of the trajectories''' is an independant stage: 
 <code> step = "ComputeTrajectories"; </code>
and provides the computation of the mechanism defined in 'mech =' in 'conditions.cpp' in the current configuration. This is not a part of the reduction process.
in 'outputs/Premixed' or 'outputs/Stochastic' the name of the file is either Premixed_.dat or Trajectory_i,  i : number of the inlet for the stochastic configuration.

# Common Errors

If you get on your screen this kind of error :
 - <code> Cantera Error </code>
or 
-  <code>...terminate called after throwing an instance of 'Cantera::CanteraError'  what():  </code>
-  <code>***********************************************************************</code>
- <code>CanteraError thrown by ct2ctml_string:</code>
-  <code>Error converting input file "./" to CTML.</code>
-  <code>Python command was: 'python'</code>
-  <code>The exit code was: 1</code>

Usually, it is due to the input_ini file : Cantera isn't able to find and read your scheme.
To make sure:
- Check every path in condition.cpp and in your Makefile
- Check and spell every name of folder/files you entered in the program (condition.cpp).
- Don't let  "mech_ref" in input_ini empty, write "mech_ref = None" instead
- Check the path of  "trajectory_ref =", if you don't want to use it, leave it empty.
