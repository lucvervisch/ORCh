
#include "optimisation.h"

#include <Cantera.h>
#include <IdealGasMix.h>    // defines class IdealGasMix
#include <equilibrium.h>    // chemical equilibrium
#include <transport.h>      // transport properties
#include <zerodim.h>
#include <user.h>

#include "../read_write/read.h"


#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <ctime>



Chromosome::Chromosome(string bits, double fitness) //Constructeur
   :m_bits(bits), m_fitness(fitness)
{}

Chromosome::~Chromosome() //Destructeur
{}




//---Optim---

Optim::Optim() //Constructeur
{}




int comparison (const void *x, const void *y)
{
   double xx = *(double*)x;
   double yy = *(double*)y;
   if (xx < yy) return -1;
   if (xx > yy) return 1;
   return 0;
}








void Optim::Optimise(string mech_ref, string mech, string mech_desc, vector<bool> QSS_Species, OptimScenario* listOptimScenarios, int nbInlets_nbFlames, vector<string> listTargets, string configuration, vector<string> trajectory_ref) 
{

   
         


   //Treat parallel stuff
   int rank, nb_processus;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nb_processus);

   int PopSize =  listOptimScenarios->m_PopSize;
   int MaxAllowableGenerations = listOptimScenarios->m_MaxAllowableGenerations;

   double CrossoverRate = listOptimScenarios->m_CrossoverRate;
   double MutationRate = listOptimScenarios->m_MutationRate;

   int Ifi = (PopSize/nb_processus)*rank;
   int Ila = (PopSize/nb_processus)*rank + (PopSize/nb_processus);


   string Dim;
   if (configuration == "MultipleInlet")
      Dim = "0D";
      
   if (configuration == "PremixedFlames")
      Dim = "1D";





   vector<Reaction_ORCh*> listReactions;
   vector<Species_ORCh*> listSpecies;

   Read *r = new Read();
   r->Read_reactions(mech, listReactions);
   r->Read_species(mech, listSpecies);

   vector<vector<double> > min_max_A (listReactions.size(), vector<double> (2));
   vector<vector<double> > min_max_b (listReactions.size(), vector<double> (2));
   vector<vector<double> > min_max_Ea_R (listReactions.size(), vector<double> (2));

   Min_Max(listOptimScenarios ,listReactions, min_max_A, min_max_b, min_max_Ea_R);

   IdealGasMix *mixture  = new IdealGasMix(mech,mech_desc);

   int nreac = mixture->nReactions();
   int nsp = mixture->nSpecies();

   vector<int> NBits_A (listReactions.size());
   vector<int> NBits_b (listReactions.size());
   vector<int> NBits_Ea_R (listReactions.size());

   vector<double> A (listReactions.size());
   vector<double> b (listReactions.size());
   vector<double> E (listReactions.size());



   for (unsigned int j=0; j<listReactions.size(); j++)
   {
         if (log10(listReactions[j]->m_A) < 1e-10 || (min_max_A[j][0] == min_max_A[j][1]))
            NBits_A[j] = 0;
         else
            NBits_A[j] = int(log((min_max_A[j][1]-min_max_A[j][0])/accuracy_A)/log(2)) +1;
         
         if (min_max_b[j][0] == min_max_b[j][1])
            NBits_b[j] = 0;
         else
            NBits_b[j] = int(log((min_max_b[j][1]-min_max_b[j][0])/accuracy_b)/log(2)) +1;

         if (min_max_Ea_R[j][0] == min_max_Ea_R[j][1])
            NBits_Ea_R[j] = 0;
         else
            NBits_Ea_R[j] = int(log((min_max_Ea_R[j][1]-min_max_Ea_R[j][0])/accuracy_Ea_R)/log(2)) +1;
   }


   //seed the random number generator
   srand((int)time(NULL));

   vector<Chromosome*> Population;
   for (int k=0; k<PopSize; k++)
   {
      Population.push_back(new Chromosome("", 0.0));
      Population[k]->m_bits = "a";
   }



   int NbTotBits = 0;
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      NbTotBits += NBits_A[j] + NBits_b[j] + NBits_Ea_R[j];
   }

   if (rank == 0)
   {
      for (int k=0; k<PopSize; k++)
      {
         string bits_random = GetRandomBits(NbTotBits);
         Population[k]->m_bits = bits_random;
      }
   }

   int count = 0;
   vector<float> mean_fitness;
   vector<float> best_fitness;
   float mean(0);

   

   //we will set this flag if a solution has been found
   bool bFound = false;
   bool conv = false;

   while(!bFound)
   {

      //Send the chromosomes computed by the root processus 0 to the other processus
      char msg[NbTotBits];

      for (int i=0; i<PopSize; i++)
      {
         if (rank == 0)
         {
            string message = Population[i]->m_bits;
            strcpy(msg, message.c_str());
         }
         MPI_Bcast(msg, NbTotBits, MPI_CHAR, 0, MPI_COMM_WORLD);
         Population[i]->m_bits = string(msg).substr(0, NbTotBits);
      }

      vector<bool> Species_to_add (nsp, true);
      int nbQSS = 0;
      for (int k=0; k<nsp; k++)
      {
         if (QSS_Species[k])
         {
            Species_to_add[k] = false;
            nbQSS += 1;
         }
      }

     for (int i=Ifi; i<Ila; i++)
     {
         Translate_Bits(A, b, E, NBits_A, NBits_b, NBits_Ea_R, Population[i]->m_bits,
                        min_max_A, min_max_b, min_max_Ea_R);



         stringstream s_rank;
         s_rank << rank;

         string path = "./analytic_schemes/Ref";


         path.append(s_rank.str());



         CreateAnalyticDirectory(configuration, path, nbInlets_nbFlames, "Optimisation", trajectory_ref);

         string chem_file = path;
         chem_file.append("/mech_QSS.h");

         Write *w = new Write();
         Write_QSS *w_qss = new Write_QSS();

         w_qss->Write_QSS_file(Dim, mech, chem_file.c_str(), QSS_Species, true, A, b, E);

         string scheme_file = path;
         scheme_file.append("/scheme.xml");
     
     
         // Commented by Huu-Tri Nguyen - 2020.02.17
       //HT  w->Write_xml_for_Analytic_Applications(mech, mech_desc, scheme_file.c_str(), Species_to_add);
         
         // If QSS, write QSS scheme (without reaction). If Optim, write OPTIM scheme
         // Added by Huu-Tri Nguyen - Written by Kaidi Wan - 2020.02.17
         if (nbQSS > 0)
         {
             w->Write_xml_for_Analytic_Applications(mech, mech_desc, scheme_file.c_str(), Species_to_add);
         }
         else
         {
             vector<bool> Reactions_to_add (nreac, true);
             for (int k=0; k<nreac; k++)
             Reactions_to_add[k] = true;
             w->Write_xml_file(mech, mech_desc, scheme_file.c_str(),
                               Species_to_add, Reactions_to_add, true, A, b, E, vector<Species_ORCh*> (), vector<Reaction_ORCh*> ());
         }
	// End add	
	
         string execute = "cd ";
         execute.append(path).append("; sh launch_cantera.sh");

         cout << rank << "      path " << execute << endl;
         //getchar();

         int ret = system(execute.c_str());
         if (ret == -1) cout << "ERROR with launch calculation";
         
           
         double fitness(0);
         if (configuration == "MultipleInlet")
           fit_function_0D(mech_ref, mech, fitness, nbInlets_nbFlames-1, listTargets, trajectory_ref,  nsp-nbQSS, "Optim", rank);

 
         if (configuration == "PremixedFlames")
           fit_function_1D(mech_ref, mech, fitness, nbInlets_nbFlames, listTargets, trajectory_ref, nsp-nbQSS, "Optim", rank);

         Population[i]->m_fitness = fitness;

       
       
         stringstream POPULATION;
         stringstream GENERATION;
         stringstream FITNESS;


//Plot only the best trajectory at every generation
         if (i == 0)
         {
            POPULATION << i;
            GENERATION << count;
            FITNESS << fitness;
            for (int inlet=0; inlet < nbInlets_nbFlames; inlet++)
            { 
               stringstream INLET;
               INLET << inlet;

               string name = "GEN";
               name.append(GENERATION.str()).append("_POP").append(POPULATION.str()).append("_FIT").append(FITNESS.str());

               if (configuration == "MultipleInlet")
                  name.append("_INLET").append(INLET.str());
               else  if (configuration == "PremixedFlames")              
                  name.append("_Flame").append(INLET.str());
 
               name.append(".eps");
         
    
               string mv = "mv ";
               mv.append(path).append("/QSS");
               stringstream s_nsp;
               s_nsp << nsp-nbQSS;
               mv.append(s_nsp.str());
               if (configuration == "MultipleInlet")         
                  mv.append("_Inlet").append(INLET.str());
               else if (configuration == "PremixedFlames")   
                  mv.append("_Flame").append(INLET.str());

               mv.append(".eps ");

               mv.append(path).append("/").append(name);
               ret = system(mv.c_str());
               if (ret == -1) cout << "ERROR system";

               mv = "mv "; 
               mv.append(path).append("/").append(name).append(" analytic_schemes/PLOTS/");
               ret = system(mv.c_str());
               if (ret == -1) cout << "ERROR system";
            }
         }
      }


      for (unsigned int i=0; i<PopSize; i++)
      {
         double send = Population[i]->m_fitness;
         double receive = 0.0;
         MPI_Reduce(&send, &receive, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         if (rank == 0)
         {
            Population[i]->m_fitness = receive;
         }
      }



      if (rank == 0)
      {
         cout << "-----------------------------------------------------------------------------------------------------" << endl;
         cout << "         New generation " << count << endl;
         cout << "-----------------------------------------------------------------------------------------------------" << endl;

         ofstream log("Optimisation.log", ios::app);
         log << "Generation n°" << count << endl;

         for (int i=0; i<PopSize; i++)
         {

            log <<  " Population[" << i << "]->m_fitness " << Population[i]->m_fitness << endl;

            cout << " Population[" << i << "]->m_fitness " << Population[i]->m_fitness << endl;
            mean +=  Population[i]->m_fitness;
	    
         }
         cout << endl << endl;
         log << endl;

         log.close();
         mean /= PopSize;
            

         mean_fitness.push_back(fabs(mean));  
         best_fitness.push_back(fabs(Population[0]->m_fitness));

         //define some temporary storage for the new population we are about to create
         vector<Chromosome*> Temporary;
         vector<Chromosome*> Sort;

         for (int k=0; k<PopSize; k++)
         {
            Temporary.push_back(new Chromosome("", 0.0));
            Sort.push_back(new Chromosome("", 0.0));
         }

         Rank(Sort, Population, PopSize);

         double Total_fitness_linear_ranking = 0.0;

         Linear_Ranking(Sort, PopSize);

         for (int i=0; i<PopSize; i++)
            Total_fitness_linear_ranking += Sort[i]->m_fitness;

         int cPop = 0;
         Elitism(cPop, Temporary, Sort, PopSize);
         cPop = NB_ELITISM;


         //loop until we have created PopSize new chromosomes
         while (cPop < PopSize)
         {
            // we are going to create the new population by grabbing members of the old population
            // two at a time via roulette wheel selection.
            string offspring1 = Roulette(int(Total_fitness_linear_ranking), Sort, PopSize);
            string offspring2 = Roulette(int(Total_fitness_linear_ranking), Sort, PopSize);

            //add crossover dependent on the crossover rate
            Crossover(offspring1, offspring2, CrossoverRate);



            //now mutate dependent on the mutation rate
            Mutate(offspring1, MutationRate);
            Mutate(offspring2, MutationRate);

           //add these offspring to the new population. (assigning zero as their
           //fitness scores)
           Temporary[cPop]->m_bits = offspring1;
           Temporary[cPop]->m_fitness = 0.0;
           cPop += 1;

           if (cPop < PopSize)
           {
              Temporary[cPop]->m_bits = offspring2;
              Temporary[cPop]->m_fitness = 0.0;
              cPop += 1;
           }
           
         }//end loop



         //copy temp population into main population array
         for (int i=0; i<PopSize; i++)
         {
            Population[i] = Temporary[i];
         }

         //Store the best solution
         
         //-----------Create new Ref ... directory----------//
 
         string createFolder = "mkdir ";
         createFolder.append("analytic_schemes/Ref/");
 
         int ret = system(createFolder.c_str());
         if (ret == -1) cout << "ERROR system";
 
             




         Translate_Bits(A, b, E, NBits_A, NBits_b, NBits_Ea_R, Population[0]->m_bits,
                        min_max_A, min_max_b, min_max_Ea_R);



         vector<bool> Species_to_add (nsp, true);
         for (int k=0; k<nsp; k++)
         {
            if (QSS_Species[k])
            {
               Species_to_add[k] = false;
               nbQSS += 1;
            }
         }

         vector<bool> Reactions_to_add (nreac, true);
         Write *w = new Write();
         Write_QSS *w_qss = new Write_QSS();

         w_qss->Write_QSS_file(Dim, mech, "analytic_schemes/Ref/best_mech_QSS.h", QSS_Species, true, A, b, E);
         
         w->Write_xml_for_Analytic_Applications(mech, mech_desc, "analytic_schemes/Ref/best_scheme.xml",
                                                   Species_to_add);

         for (int k=0; k<nsp; k++)
            Species_to_add[k] = true;
         for (int k=0; k<nreac; k++)
            Reactions_to_add[k] = true;

       

	// Original write xml file - Only write the last generation - Commented by Huu-Tri NGUYEN - 17/02/2020
	//HT  w->Write_xml_file(mech, mech_desc, "analytic_schemes/Ref/chemistry_for_restart.xml",
        //HT                             Species_to_add, Reactions_to_add, true, A, b, E, vector<Species_ORCh*> (), vector<Reaction_ORCh*> ());

	// Write xml file for each GEN - Add 17/02/2020
	stringstream str_scheme;
	stringstream str_count;  //Which Generation
	str_count << count;
	str_scheme << "analytic_schemes/Ref/GEN_" << str_count.str() << "_chemistry_for_restart.xml";
	w->Write_xml_file(mech, mech_desc, str_scheme.str(),
                                     Species_to_add, Reactions_to_add, true, A, b, E, vector<Species_ORCh*> (), vector<Reaction_ORCh*> ());
	// End Add

//write a gnuplot script for the mean and min fitness

// mean and best fitness stored

         string data;
         data = "analytic_schemes/PLOTS/fitness.dat";

         ofstream data_fitness(data.c_str());
         for (int i=1;i<=count;i++)
            data_fitness << i  << "   " << mean_fitness[i] << "   "  << best_fitness[i]  << endl;
         data_fitness.close();


//gnuplot fitness
 
          string path = "analytic_schemes/PLOTS/";    
 
          plotFitness(path);


          string graph = "cd ";
          graph.append(path).append("; ").append("gnuplot make.gnu");
          int retg = system(graph.c_str());
          if (retg == -1) cout << "ERROR system";


//Arrêt si pendant 20 générations la fitness ne diminue pas
/*      if (count >20)
      {
         if (best_fitness[count] == best_fitness[count-20])
         {  
            cout << "Optimisation is converged." << endl;
            conv = true;
         }
      }

*/
      } //end if rank == 0



      ++count;



//Ou arrêt quand toutes les générations ont été calculées
      if ( count > MaxAllowableGenerations)
         bFound = true;



 }
















}


//----------Translate_Bits----------

void Optim::Translate_Bits(vector<double>& A, vector<double>& b, vector<double>& E, 
   vector<int> NBits_A, vector<int> NBits_b, vector<int> NBits_Ea_R, string Bits,
   vector<vector<double> > min_max_A, 
   vector<vector<double> > min_max_b,
   vector<vector<double> > min_max_Ea_R) 
{

  int cut = 0;
  double Variable_dec;

  for (unsigned int j=0; j<NBits_A.size(); j++)
  {

        //----------A----------
        if (NBits_A[j] != 0)
        {
           Variable_dec = BinToDec(Bits.substr(cut, NBits_A[j]));
   
           A[j] = pow(10, (min_max_A[j][0] + 
              (Variable_dec/(pow(2.0,NBits_A[j])-1))*(min_max_A[j][1]-min_max_A[j][0])));
   
           cut += NBits_A[j];
        }
        else
        {
           A[j] = pow(10, (min_max_A[j][1]+min_max_A[j][0])/2);
        }
   
        //----------b----------
        if (NBits_b[j] != 0)
        {
           Variable_dec = BinToDec(Bits.substr(cut, NBits_b[j]));
   
           b[j] = min_max_b[j][0] +
              (Variable_dec/(pow(2.0,NBits_b[j])-1))*(min_max_b[j][1]-min_max_b[j][0]);
   
           cut += NBits_b[j];
        }
        else
        {
           b[j] = (min_max_b[j][1]+min_max_b[j][0])/2;
        }

   
        //----------Ea_R----------
        if (NBits_Ea_R[j] != 0)
        {
           Variable_dec = BinToDec(Bits.substr(cut, NBits_Ea_R[j]));
   
           E[j] = (min_max_Ea_R[j][0] +
              (Variable_dec/(pow(2.0,NBits_Ea_R[j])-1))*(min_max_Ea_R[j][1]-min_max_Ea_R[j][0]))*8.314;
   
           cut += NBits_Ea_R[j];
        }
        else
        {
           E[j] = (min_max_Ea_R[j][1]+min_max_Ea_R[j][0])*8.314/2;
        }
  }


}



//----------Min_Max----------

void Optim::Min_Max(OptimScenario* listOptimScenarios, vector<Reaction_ORCh*>& listReactions, 
   vector<vector<double> >& min_max_A, 
   vector<vector<double> >& min_max_b, 
   vector<vector<double> >& min_max_Ea_R) 
{

   double AllowedVariation_A = listOptimScenarios->m_AllowedVariation_A;
   double AllowedVariation_b = listOptimScenarios->m_AllowedVariation_b;
   double AllowedVariation_E = listOptimScenarios->m_AllowedVariation_E;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
//   if (j==212||j==135||j==189||j==7||j==14||j==13||j==208||j==142||j==214||j==210||j==211||j==5||j==199||j==6||j==209||j==140||j==237||j==240||j==192||j==221||j==15||j==8||j==141||j==21){
 //   if (j==2){

         min_max_A[j][0] = (1-AllowedVariation_A)*log10(listReactions[j]->m_A);
         min_max_A[j][1] = (1+AllowedVariation_A)*log10(listReactions[j]->m_A);
                  
         
         if (log10(listReactions[j]->m_A) < 0.0)
         {
            min_max_A[j][0] = log10(listReactions[j]->m_A);
            min_max_A[j][1] = log10(listReactions[j]->m_A);
         }

         if (listReactions[j]->m_E > 0)
         {
            min_max_Ea_R[j][0] = (1-AllowedVariation_E)*listReactions[j]->m_E/8.314;
            min_max_Ea_R[j][1] = (1+AllowedVariation_E)*listReactions[j]->m_E/8.314;
         }
         else
         {
            min_max_Ea_R[j][0] = (1+AllowedVariation_E)*listReactions[j]->m_E/8.314;
            min_max_Ea_R[j][1] = (1-AllowedVariation_E)*listReactions[j]->m_E/8.314;
         }

         if (listReactions[j]->m_E == 0.0)
         {
            min_max_Ea_R[j][0] = 0.0;
            min_max_Ea_R[j][1] = 0.0;
         }

         if (listReactions[j]->m_b > 0)
         {
            min_max_b[j][0] = (1-AllowedVariation_b)*listReactions[j]->m_b;
            min_max_b[j][1] = (1+AllowedVariation_b)*listReactions[j]->m_b;
         }
         else
         {
            min_max_b[j][0] = (1+AllowedVariation_b)*listReactions[j]->m_b;
            min_max_b[j][1] = (1-AllowedVariation_b)*listReactions[j]->m_b;
         }

         if (listReactions[j]->m_b == 0.0)
         {
            min_max_b[j][0] = 0.0;
            min_max_b[j][1] = 0.0;
         }

   }


// 
/*else {

 min_max_A[j][0] = log10(listReactions[j]->m_A);
 min_max_A[j][1] = log10(listReactions[j]->m_A);
 min_max_b[j][0] = listReactions[j]->m_b;
 min_max_b[j][1] = listReactions[j]->m_b;
 min_max_Ea_R[j][0] = listReactions[j]->m_E/8.314;
 min_max_Ea_R[j][1] = listReactions[j]->m_E/8.314;

     }
   }*/
//



}


//----------Rank----------

void Optim::Rank(vector<Chromosome*>& sort, vector<Chromosome*> Population, int PopSize) 
{
   double sort_ref[PopSize][2];

   for (int i=0; i<PopSize; i++)
   {
      sort_ref[i][0] = Population[i]->m_fitness;
      sort_ref[i][1] = i;
   }

   qsort (sort_ref, sizeof(sort_ref)/sizeof(sort_ref[0]), sizeof(sort_ref[0]), comparison);

   int Reference;
   for (int j=0; j<PopSize; j++)
   {
      Reference = int(sort_ref[j][1]);
      sort[j]->m_bits = Population[Reference]->m_bits;
      sort[j]->m_fitness = Population[Reference]->m_fitness;
   }

}


//----------Linear_Ranking----------

void Optim::Linear_Ranking(vector<Chromosome*>& Sort, int PopSize) 
{
   double SP = 2.0;
   for (int i=0; i<PopSize; i++)
   {
      Sort[i]->m_fitness =  2.0-SP+2.0*(SP-1.0)*(i)/(PopSize-1.0);
   }
}


//----------Roulette----------

string Optim::Roulette(int TotalFitness, vector<Chromosome*> Population, int PopSize) 
{
    //generate a random number between 0 & total fitness count
    double Rand = Get_random();
    float Slice = (float)(Rand * TotalFitness);

    //go through the chromosones adding up the fitness so far
    float FitnessSoFar = 0.0f;

    for (int i=0; i<PopSize; i++)
    {
        FitnessSoFar += Population[i]->m_fitness;

        //if the fitness so far > random number return the chromo at this point
        if (FitnessSoFar >= Slice)

            return Population[i]->m_bits;
    }

    return "";
}


//----------Get_random----------

double Optim::Get_random () 
{
    int random = rand() % 10000;
    double random_number = double(random)/10000.0;
    return random_number;
}


//----------GetRandomBits----------

string  Optim::GetRandomBits(int length) 
{
    string bits;

    for (int i=0; i<length; i++)
    {
    double Rand = Get_random();
        if (Rand > 0.5f)
            bits += "1";
        else
            bits += "0";
    }

    return bits;
}


//----------BinToDec----------

double Optim::BinToDec(string bits) 
{
    double val = 0;
    double value_to_add = 1;

    for (int i = bits.length(); i > 0; i--)
    {
        if (bits.at(i-1) == '1')
            val += value_to_add;
        value_to_add *= 2;
    }

    return val;
}


//----------DecToBin----------

string Optim::DecToBin(int number) 
{
    string bits;
    //cout<<"Le nombre proposé :: "<<number<<endl;
    bool loop = true;

    while (loop)
    {
       if (number%2 == 1)
          bits = "1" + bits;
       else
          bits = "0" + bits ;

       number = number/2;
       //cout<<"Le nombre après division entière par 2 :: "<<number<<endl;

       if (number == 1)
       {
          loop = false;
          bits = "1" + bits;
       }
       else if (number == 0)
       {
          loop = false;
          bits = "0" + bits;
       }
    }
    //cout<<"Le retour binaire obtenu  :: "<<bits<<endl;
    return bits;
}


//----------Crossover----------

void Optim::Crossover(string &offspring1, string &offspring2, double CrossoverRate) 
{

  double Rand1 = Get_random();

  if (Rand1 < CrossoverRate)
  {
     //Choix aléatoire du nombre de points de coupure
     double random_Crossover = Get_random();
     int NB_CROSSOVER = floor(random_Crossover*10) + 1;

     string t1;
     string t2;

     double Rand [NB_CROSSOVER];
     for (int k=0; k<NB_CROSSOVER; k++)
     {
        Rand[k] = Get_random();
     }

     qsort (Rand, sizeof(Rand)/sizeof(Rand[0]), sizeof(Rand[0]), comparison);

     vector<int> crossover (NB_CROSSOVER);
     for (int k=0; k<NB_CROSSOVER; k++)
     {
        crossover[k] = (int) (Rand[k]*offspring1.size());
     }


     t1 = offspring1.substr(0, crossover[0]);
     t2 = offspring2.substr(0, crossover[0]);

     for (int k=1; k<NB_CROSSOVER; k++)
     {
        if (k%2 == 0 /*pair*/)
        {
           t1 += offspring1.substr(crossover[k-1], crossover[k]-crossover[k-1]);
           t2 += offspring2.substr(crossover[k-1], crossover[k]-crossover[k-1]);
        }

        if (k%2 == 1 /*impair*/)
        {
           t1 += offspring2.substr(crossover[k-1], crossover[k]-crossover[k-1]);
           t2 += offspring1.substr(crossover[k-1], crossover[k]-crossover[k-1]);
        }
     }


     if (NB_CROSSOVER%2 == 0 /*pair*/)
     {
        t1 += offspring1.substr(crossover[NB_CROSSOVER-1], offspring1.size()-crossover[NB_CROSSOVER-1]);
        t2 += offspring2.substr(crossover[NB_CROSSOVER-1], offspring1.size()-crossover[NB_CROSSOVER-1]);
     }

     if (NB_CROSSOVER%2 == 1 /*impair*/)
     {
        t1 += offspring2.substr(crossover[NB_CROSSOVER-1], offspring1.size()-crossover[NB_CROSSOVER-1]);
        t2 += offspring1.substr(crossover[NB_CROSSOVER-1], offspring1.size()-crossover[NB_CROSSOVER-1]);
     }

     offspring1 = t1; offspring2 = t2;
  }


}


//----------Mutate----------

void Optim::Mutate(string &bits, double MutationRate) 
{
    for (unsigned int i=0; i<bits.length(); i++)
    {
        double Rand = Get_random();
        if (Rand < MutationRate)
        {
            if (bits.at(i) == '1')
                bits.at(i) = '0';

            else
                bits.at(i) = '1';
        }
    }

    return;
}


//----------Elitism----------

void Optim::Elitism(int cPop, vector<Chromosome*>& Temporary, vector<Chromosome*> Sort, int PopSize) 
{
             for (int j=0; j<NB_ELITISM; j++)
             {
                Temporary[cPop++] = Sort[PopSize-(1+j)];
             }
}


//----------replace----------

bool Optim::replace(std::string& str, const std::string& from, const std::string& to) 
{
   size_t start_pos = str.find(from);
   if (start_pos == std::string::npos)
      return false;
   str.replace(start_pos, from.length(), to);
   return true;
}


Optim::~Optim() //Destructeur
{}



