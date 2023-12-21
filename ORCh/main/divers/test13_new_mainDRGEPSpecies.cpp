#include "mainDRGEPSpecies.h"

//----------------------------------------
//   <mainDRGEPSpecies> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis by
//            1- Computing reference trajectories
//            2- Computing the species DRGEP coefficients (R_AD)
//            3- 
//   IN: 
//         configuration: Either "Premixed" or "MultipleInlet" case
//   OUT:
//
//----------------------------------------
void drgepSpecies(int debuglevel, vector<string> speciesToPlot, vector<string> listTargets,  string configuration, string initial_mech, string mech_desc, vector<MultipleInlet*> listInlets, vector<PremixedFlames*> listFlames, vector<bool> Targets, bool new_mixing, bool plot_T, bool plot_U, vector<string> trajectory_ref,string mech_ref)
{

   cout << endl;
                        cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl << endl;
                        cout << "------------------------------------------------ STARTING DRGEP SPECIES STEP -----------------------------------------------------" << endl << endl;
                        cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl << endl;



   print_to_screen_with_newline("MECHANISM:------------------------------------------------------------------------------------------------------------------------", 0, debuglevel);
   print_to_screen("            Reading initial mechanism \"", 0, debuglevel);
   print_to_screen(initial_mech, 0, debuglevel);
   print_to_screen("\" with description \"", 0, debuglevel);
   print_to_screen(mech_desc, 0, debuglevel);
   print_to_screen("\" ----------> ", 0, debuglevel);
 


   IdealGasMix *mixture  = new IdealGasMix(initial_mech,mech_desc);
   print_to_screen_with_newline("OK", 0, debuglevel);

   int nsp_init = mixture->nSpecies();
   int nreac_init = mixture->nReactions();
   int nbFlame = listFlames.size()-1;   

   stringstream s_nbFlame;
   s_nbFlame << nbFlame;

   print_to_screen("               Number of species: ", 0, debuglevel);
   print_to_screen_with_newline(nsp_init, 0, debuglevel);
   print_to_screen("               Number of reactions: ", 0, debuglevel);
   print_to_screen_with_newline(nreac_init, 0, debuglevel);
   print_to_screen_with_newline("", 0, debuglevel);

   vector<Species*> listSpecies_init;
   vector<Reaction*> listReactions_init;

   Read *r = new Read();
   r->Read_species(initial_mech, listSpecies_init);
   r->Read_reactions(initial_mech, listReactions_init);

   int nbIterations = 0;
   int nbInlets = 0;

   if (listInlets.size() > 1)
   {
      nbInlets = listInlets.size()-1; //minus 1 because last inlet is always for burned gases
      nbIterations =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets])->m_NbIterations;
   }

   vector<vector<vector<double> > > R_AD_Trajectories (nbInlets, vector<vector<double> > (nsp_init, vector<double>(nsp_init, 0.0)));
   vector<vector<double> > R_AD_Premixed (nsp_init, vector<double>(nsp_init, 0.0));
   vector<vector<double> > max_j_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<double> > max_jf_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<double> > max_jr_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<vector<double> > > Production_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<vector<vector<double> > > Consumption_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<vector<double> > T_Trajectories_ref (nbInlets, vector<double> (nbIterations, 0.0));
   vector<vector<vector<double> > > Ym_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<double> time_store (nbIterations, 0.0);
   vector<bool> SpeciesIntoReactants (nsp_init, false);
   vector<vector<double> > QSS_Criteria (nsp_init, vector<double> (2, 0.0));

   double sort_R_AD [nsp_init][2]; //sort_R_AD[k][0] -> interaction coefficient //sort_R_AD[k][1] -> species number
   int e = 1; //epsilon
   double T_final_ref(0);

   ofstream log("DRGEP_Species.log");

   if (configuration == "MultipleInlet")
   {
      for (int k=0; k<nsp_init; k++)
      {
         sort_R_AD[k][0] = 0.0;
         sort_R_AD[k][1] = 0;
      }

      string outputName = "./outputs/Stochastic/Ref_DRGEP_Species";
      stringstream s_nbSpeciesInit;
      s_nbSpeciesInit << nsp_init;
      outputName.append(s_nbSpeciesInit.str()).append("_");

      computeMultipleInlet *c = new computeMultipleInlet();
      c->getMultipleInlet(initial_mech, 
                          mech_desc, 
                          listInlets, 
                          Targets, 
                          new_mixing, 
                          "DRGEP_Species", 
                          R_AD_Trajectories, 
                          max_j_on_Target, 
                          Ym_Trajectories_ref, 
                          Production_Trajectories_ref, 
                          Consumption_Trajectories_ref, 
                          T_Trajectories_ref, 
                          time_store, 
                          SpeciesIntoReactants);

      OutputDeterministicTrajectories (listInlets.size(), nbIterations, listSpecies_init, outputName, Ym_Trajectories_ref, T_Trajectories_ref, time_store);

         T_final_ref = T_Trajectories_ref[0][nbIterations-1];

         double r1(0);
         double r2(0);
         double list_real [nsp_init][2]; //list_real[k][0] -> interaction coefficient //list_real[k][1] -> species number
         double list_real0 [nsp_init][2]; //list_real0[k][0] -> interaction coefficient //list_real0[k][1] -> species number

         for (int k=0; k<nsp_init; k++)
         {
            list_real0[k][0] = 0;
            list_real0[k][1] = 0;
         } 


         for (int ka=0; ka<nsp_init; ka++)
         {
   double interm_sort_R_AD [nsp_init][2]; //interm_sort_R_AD[k][0] -> interaction coefficient //interm_sort_R_AD[k][1] -> species number
      for (int k=0; k<nsp_init; k++)
      {
         interm_sort_R_AD[k][0] = 0.0;
         interm_sort_R_AD[k][1] = k;
      }
            if (Targets[ka])
            {
               for (int kb=0; kb<nsp_init; kb++)
               {
                  for (unsigned int n=0; n<nbInlets; n++)
                  {
                     if (interm_sort_R_AD[kb][0] < R_AD_Trajectories[n][ka][kb])
                        interm_sort_R_AD[kb][0] = R_AD_Trajectories[n][ka][kb];

                  }
               }
               //If the species is within the targets or within the initial reactants its R_AD coefficient is set to one to ensure it is kept while reducing the mechanism
               int s = nsp_init; 
               for (int k=0; k<nsp_init; k++)
               {
                  if (SpeciesIntoReactants[k] && !Targets[k])
                  {    
                     interm_sort_R_AD[k][0] = 1;
                     s--;
                  }
                  if (Targets[k] && !SpeciesIntoReactants[k])
                  {    
                     interm_sort_R_AD[k][0] = 1;
                     s--;
                  }
                  if (Targets[k] && SpeciesIntoReactants[k])
                  {    
                     interm_sort_R_AD[k][0] = 1;
                     s--;
                  }
               }
               //Sort the interm_R_AD vector
               qsort (interm_sort_R_AD, sizeof(interm_sort_R_AD)/sizeof(interm_sort_R_AD[0]), sizeof(interm_sort_R_AD[0]), compare_numbers);


               //Keep in memory the original coefficients
               for (int k=0; k<nsp_init; k++)
               {
                  list_real[k][0] = interm_sort_R_AD[k][0];
                  list_real[k][1] = interm_sort_R_AD[k][1];
               } 



               //Ranking of the species depending on their position for each target
               for (int kb=0; kb<nsp_init; kb++)
               {
                  double a = s + kb;
                  double b = s;
                  if (interm_sort_R_AD[kb][0]!=1)
                     interm_sort_R_AD[kb][0] = abs(a/b -1);
               }
   for (int k=0; k<nsp_init; k++)
   {
      cout << interm_sort_R_AD[k][0] << "  " << listSpecies_init[interm_sort_R_AD[k][1]]->m_Name << endl;
      log << interm_sort_R_AD[k][0] << "  " << listSpecies_init[interm_sort_R_AD[k][1]]->m_Name << endl;
   }

//               for (int kb=0; kb<nsp_init; kb++)
//               {
//                  if (sort_R_AD[kb][1] <= interm_sort_R_AD[kb][1])
//                  {
//                     sort_R_AD[kb][1] = interm_sort_R_AD[kb][1];
//                     sort_R_AD[kb][0] = interm_sort_R_AD[kb][0];
//                  }
//               }



//               //marche bien pour 1 seule target
//               for (int kb=nsp_init-1; kb>=0; kb--)
//               {
//                  if (interm_sort_R_AD[kb][0] == 1)
//                  {
//                  sort_R_AD[kb][1] = interm_sort_R_AD[kb][1];
//                  sort_R_AD[kb][0] = interm_sort_R_AD[kb][0];
//                  }
//                  else if (interm_sort_R_AD[kb][0] != 1)
//                  {
//                     if (abs(sort_R_AD[kb][1] - interm_sort_R_AD[kb][1]) < e)
//                     {
//                        sort_R_AD[kb][1] = interm_sort_R_AD[kb][1];
//                        sort_R_AD[kb][0] = interm_sort_R_AD[kb][0];
//                     }     
//                     else if (abs(sort_R_AD[kb][1] - interm_sort_R_AD[kb][1]) >= e)
//                     {    
//                        for (int kr=0; kr<nsp_init; kr++)
//                        {
//                           if (abs(interm_sort_R_AD[kb][1] -  list_real[kr][1]) < e)
//                           {
//                              r1 = list_real[kr][0];
//                           }   
//                           if (abs(list_real0[kr][1] - sort_R_AD[kb][1]) < e)
//                           {
//                              r2 = list_real0[kr][0];
//                           }   
//                        }
//                        if (r2 <= r1)
//                        {
//                           sort_R_AD[kb][1] = interm_sort_R_AD[kb][1];
//                           sort_R_AD[kb][0] = interm_sort_R_AD[kb][0];
//                           list_real0[kb][0] = list_real[kb][0];
//                           list_real0[kb][1] = list_real[kb][1];
//                        }
//                     }
//                  }
//               }
               //test 2 seule target
               for (int kb=nsp_init-1; kb>=0; kb--)
               {
log  << "kb : " << kb << endl;
                     if (abs(sort_R_AD[kb][1] - interm_sort_R_AD[kb][1]) < e)
                     {
log  << "egalité"  << endl;
                        sort_R_AD[kb][1] = interm_sort_R_AD[kb][1];
                        sort_R_AD[kb][0] = interm_sort_R_AD[kb][0];
                     }
                     else
                     {
                        for (int k=0; k<kb; k++)
                        {
                           if (sort_R_AD[kb][1] == sort_R_AD[k][1])
                           {
log  << "sort deja pris"  << endl;
                              sort_R_AD[kb][1] = interm_sort_R_AD[kb][1];
                              sort_R_AD[kb][0] = interm_sort_R_AD[kb][0];
                           }
                           else if (sort_R_AD[k][1] != interm_sort_R_AD[kb][1] && sort_R_AD[kb][1] != sort_R_AD[k][1])
                           {
log  << "nouveaux et diff"  << endl;
                              for (int kr=0; kr<nsp_init; kr++)
                              {
                                 if (abs(interm_sort_R_AD[kb][1] -  list_real[kr][1]) < e)
                                 {
                                    r1 = list_real[kr][0];
                                 }
                                 if (abs(list_real0[kr][1] - sort_R_AD[kb][1]) < e)
                                 {
                                    r2 = list_real0[kr][0];
                                 }
                              }
                              if (r2 <= r1)
                              {
                                 sort_R_AD[kb][1] = interm_sort_R_AD[kb][1];
                                 sort_R_AD[kb][0] = interm_sort_R_AD[kb][0];
                              }

                           } 
if (sort_R_AD[k][1] == interm_sort_R_AD[kb][1])
log << "interm pris" << endl;
                        }
                     }
               }



               for (int kb=0; kb<nsp_init; kb++)
               {
                  list_real0[kb][0] = sort_R_AD[kb][0];
                  list_real0[kb][1] = sort_R_AD[kb][1];

               }



            }
         }



   } //end if (configuration == "MultipleInlet")

   if (configuration == "PremixedFlames")
   {
      for (int k=0; k<nsp_init; k++)
      {
         sort_R_AD[k][0] = 0.0;
         sort_R_AD[k][1] = k;
      }

      string outputName = "./outputs/Premixed/Ref_DRGEP_Species";
      stringstream s_nbSpeciesInit;
      s_nbSpeciesInit << nsp_init;
      outputName.append(s_nbSpeciesInit.str());

      computePremixedFlames(initial_mech, 
                            outputName, 
                            mech_desc, 
                            listFlames, 
                            Targets, 
                            "DRGEP_Species", 
                            SpeciesIntoReactants, 
                            R_AD_Premixed, 
                            max_j_on_Target, 
                            max_jf_on_Target, 
                            max_jr_on_Target, 
                            QSS_Criteria, 
                            true);


      for (int ka=0; ka<nsp_init; ka++) //BE CAREFUL, FOR NOW DRGEP ANALYSIS ON PREMIXED FLAMES IS CODED FOR ONLY ONE FLAME, JUST NEED TO LOOP ON THE NFLAMES
      {
         if (Targets[ka])
         {
            for (int kb=0; kb<nsp_init; kb++)
            {
               if (sort_R_AD[kb][0] < R_AD_Premixed[ka][kb])
                  sort_R_AD[kb][0] = R_AD_Premixed[ka][kb];
            }
         }
      }
   } //end if (configuration == "PremixedFlames")

  



  
//   qsort (sort_R_AD, sizeof(sort_R_AD)/sizeof(sort_R_AD[0]), sizeof(sort_R_AD[0]), compare_numbers);

   log << "mech = " << initial_mech << endl;
   log << "mech_ref = " << mech_ref << endl;
   log << "trajectory_ref = " << trajectory_ref[0] << endl;
   log << "DRGEP coefficients ---------->" << endl; 

   for (int k=0; k<nsp_init; k++)
   {
      cout << sort_R_AD[k][0] << "  " << listSpecies_init[sort_R_AD[k][1]]->m_Name << endl;
      log << sort_R_AD[k][0] << "  " << listSpecies_init[sort_R_AD[k][1]]->m_Name << endl;
   }



   cout << endl;

   for (int nbSpeciesToKeep=nsp_init-1; nbSpeciesToKeep>9; nbSpeciesToKeep--)
   {
      vector<bool> Species_to_add (nsp_init, false);
      vector<bool> Reactions_to_add (nreac_init, true);


      string outputSchemeName = "./outputs/mechanisms/drgepSpecies";
      stringstream s_nbSpeciesToKeep;
      s_nbSpeciesToKeep << nbSpeciesToKeep;
      outputSchemeName.append(s_nbSpeciesToKeep.str()).append(".xml");

      Write *w = new Write();
      w->Write_xml_file(initial_mech, mech_desc, outputSchemeName, Species_to_add, Reactions_to_add, false, 
                        vector<double> (nreac_init, 0.0), vector<double> (nreac_init, 0.0), vector<double> (nreac_init, 0.0), vector<Species*> (), vector<Reaction*> ());
   }

   ofstream NbSpecies_fitness("./outputs/NbSpecies_fitness.dat");

double T_final(0);
double Te = 50; //50K limit between the final ref temperature and the final reduced temperature
double fitness_final = 0.0;


  
   for (int nbSpeciesToKeep=nsp_init-1; nbSpeciesToKeep>9; nbSpeciesToKeep--)
   {
   if (fitness_final < Te)
   {
      string outputSchemeName = "./outputs/mechanisms/drgepSpecies";
      stringstream s_nbSpeciesToKeep;
      s_nbSpeciesToKeep << nbSpeciesToKeep;
      outputSchemeName.append(s_nbSpeciesToKeep.str()).append(".xml");

      vector<Species*> listSpecies_drgep;
      Read *r = new Read();
      r->Read_species(outputSchemeName, listSpecies_drgep);
      int nsp_drgep = listSpecies_drgep.size();


      if (configuration == "MultipleInlet")
      {


         string outputName = "./outputs/Stochastic/Reduced_DRGEP_Species";
         outputName.append(s_nbSpeciesToKeep.str()).append("_");
 
         vector<vector<vector<double> > > Ym_Trajectories (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_drgep, 0.0)));
         vector<vector<double> > T_Trajectories (nbInlets, vector<double> (nbIterations, 0.0));

         
         computeMultipleInlet *c = new computeMultipleInlet();
         c->getMultipleInlet(outputSchemeName, 
                             mech_desc, 
                             listInlets, 
                             Targets, 
                             false, 
                             "ComputeTrajectories", 
                             R_AD_Trajectories, 
                             max_j_on_Target, 
                             Ym_Trajectories, 
                             Production_Trajectories_ref, 
                             Consumption_Trajectories_ref, 
                             T_Trajectories, 
                             time_store, 
                             SpeciesIntoReactants);


         T_final = T_Trajectories[0][nbIterations-1];
         fitness_final = abs(T_final - T_final_ref);

         OutputDeterministicTrajectories (listInlets.size(), nbIterations, listSpecies_drgep, outputName, Ym_Trajectories, T_Trajectories, time_store);

         double fitness = 0.0;
         double fitness_equilibrium = 0.0;
         for (int kinit=0; kinit<nsp_init; kinit++)
         {
            if (Targets[kinit])
            {
               for (int k=0; k<nsp_drgep; k++)
               {
                  if (listSpecies_drgep[k]->m_Name == listSpecies_init[kinit]->m_Name)
                  {
                     for (unsigned int n=0; n<nbInlets; n++)
                     {
                        for (int i=0; i<nbIterations; i++)
                        {
                           if (Ym_Trajectories_ref[n][i][kinit] > 1e-08)
                              fitness -= fabs(Ym_Trajectories[n][i][k]-Ym_Trajectories_ref[n][i][kinit])/Ym_Trajectories_ref[n][i][kinit];
                        }
                        fitness_equilibrium -= fabs(Ym_Trajectories[n][nbIterations-1][k]-Ym_Trajectories_ref[n][nbIterations-1][kinit])/Ym_Trajectories_ref[n][nbIterations-1][kinit];
                     }
                  }
               }
            }
         }

         NbSpecies_fitness << nbSpeciesToKeep << "species, fitness :" << fitness << " and fitness at equilibrium " << fitness_equilibrium << endl;



      
      } //end if configuration == "MultipleInlet"

      if (configuration == "PremixedFlames")
      {

         string outputName = "./outputs/Premixed/Reduced_DRGEP_Species";
         outputName.append(s_nbSpeciesToKeep.str()).append("_").append(s_nbFlame.str());

         computePremixedFlames(outputSchemeName, 
                               outputName, 
                               mech_desc, 
                               listFlames, 
                               Targets, 
                               "ComputeTrajectories", 
                               SpeciesIntoReactants, 
                               R_AD_Premixed, 
                               max_j_on_Target, 
                               max_jf_on_Target, 
                               max_jr_on_Target, 
                               QSS_Criteria, 
                               false);
      }
      
  



        Script_gnuplot("DRGEP_Species", initial_mech, speciesToPlot, configuration, mech_desc, nbSpeciesToKeep, plot_U, plot_T, trajectory_ref, mech_ref, nbInlets, listTargets ) ;




      string path="outputs/";
      if (configuration == "MultipleInlet")
         path.append("Stochastic/");
      if (configuration == "PremixedFlames")
         path.append("Premixed/");
      string graph = "cd ";
      graph.append(path).append("; ").append("gnuplot makeGnu.gnu");
      int ret = system(graph.c_str());
      if (ret == -1) cout << "ERROR system";
	

   } //end if
   else 
   {
   cout << "DRGEP Species step terminated. For next step, use the scheme \"drgepSpecies" << nbSpeciesToKeep+1  <<"\" located in output/mechanisms." << endl; 
   break;
   }
 
  }// end for nbTokeep
  NbSpecies_fitness.close();





   log.close();




}





