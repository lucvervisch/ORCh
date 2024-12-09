#include <fstream>
#include <sstream>

#include "math.h"
#include "write_QSS.h"

//---Write_QSS---

Write_QSS::Write_QSS() //Constructeur
{}


void Write_QSS::Write_QSS_file(string Dimension, string mech, string write_mech, vector<bool> &QSS_Species, bool optim, vector<double> A_Arrh, vector<double> b_Arrh, vector<double> E_Arrh) const
{





   vector<Species_ORCh*> listSpecies;
   vector<Reaction_ORCh*> listReactions;

   Read *r = new Read();
   r->Read_species(mech, listSpecies);
   r->Read_reactions(mech, listReactions);

   int nsp = listSpecies.size();
   int nreac = listReactions.size();


   vector<double> Keep_fwd_reaction (nreac, true);
   vector<double> Keep_rev_reaction (nreac, true);

   if (optim)
   {
      for (unsigned int j=0; j<listReactions.size(); j++)
      {
         listReactions[j]->m_A = A_Arrh[j];
         listReactions[j]->m_b = b_Arrh[j];
         listReactions[j]->m_E = E_Arrh[j];
     }
   }
















   vector<int> Index_QSS_Species_number;
   int count_nQSS = 0;
   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k] == false)
      {
         Index_QSS_Species_number.push_back(count_nQSS);
         count_nQSS += 1;
      }
      else
         Index_QSS_Species_number.push_back(-1);
   }


   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k])
         listSpecies[k]->m_QSS = true;
   }

   // Count the number of QSS species
   // Added by Huu-Tri Nguyen - Written by Kaidi Wan - 2020.02.17
   int nQSS_Species = 0;
   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k])
      {
         nQSS_Species += 1;
      }
   }
   // End add

   cout << "Analytic scheme 1/2 : " <<  write_mech << endl;

   ofstream write(write_mech.c_str(), ofstream::out);
   write.precision(7);


   write << endl << endl;
   write << "void Methane(vector_fp& get_ropnet_4steps){ }" << endl;
   write << "int Kerosene(int x){ return x+1; }" << endl;

   write << endl << endl;
   write << "extern IdealGasMix* mixture;" << endl;
//   write << "extern string mech;" << endl;


   if (Dimension == "1D")
   {
      write << endl << endl;
      write << "extern StFlow* flow;" << endl;

      write << endl << endl << endl << endl;

   // By pass QSS step if there is no QSS specie
   // Added by Huu-Tri Nguyen - Written by Kaidi Wan - 2020.02.17
      if (nQSS_Species > 0)
      {
         write << "int Mechanism = 1;" << endl;	//Analytical calculation oneD (Cantera)
      }
      else
      {
         write << "int Mechanism = 2;" << endl;	// There is no 2D case >> By pass analytical calculation
      }
   // End add

      write << endl << endl << endl;
      write << "extern int optimisation_count;" << endl;
      write << endl << endl << endl;


      write << "void Reduced(doublereal *w)" << endl;
      write << "{" << endl;
      write << endl << endl;
      write << "   int nsp = flow->phase().nSpecies();" << endl;
      write << "   doublereal T = flow->phase().temperature();" << endl;
      write << endl;
      write << "   vector_fp C;" << endl;
      write << "   C.resize(nsp);" << endl;
      write << "   flow->phase().getActivityConcentrations(&C[0]);" << endl;
   }

   if (Dimension == "0D")
   {
      write << endl << endl << endl << endl;

   // By pass QSS step if there is no QSS specie
   // Added by Huu-Tri Nguyen - Written by Kaidi Wan - 2020.02.17
      if (nQSS_Species > 0)
      {
         write << "int Mechanism = 0;" << endl;	//Analytical calculation zeroD (Cantera)
      }
      else
      {
         write << "int Mechanism = 2;" << endl;	// There is no 2D case >> By pass analytical calculation
      }
   // End add

      write << endl << endl << endl;
      write << "//----------Function for 1D----------" << endl;
      write << "void Reduced(doublereal *w)" << endl;
      write << "{" << endl;
      write << "}" << endl;
      write << endl << endl << endl;
      write << "      //----------Function for 0D----------" << endl;
      write << "void Reduced(doublereal *w, doublereal T, doublereal *mass, doublereal density, const doublereal *mw)" << endl;
      write << "{" << endl;
      write << "   int nsp = mixture->nSpecies();" << endl;
      write << "   vector<Species_ORCh*> listSpecies;;" << endl;
      write << "   " << endl;
//      write << "  Read *r = new Read(); " << endl;
//      write << "  r->Read_species(mech, listSpecies); " << endl;
      write << "   " << endl;
      write << "   " << endl;
      write << "   vector<double> Y (nsp, 0.0);" << endl;
      write << "   vector<double> C (nsp, 0.0);" << endl;
      write << "" << endl;
      write << "   double total_mass = 0.0;" << endl;
      write << "   for (int k=0; k<nsp; k++)" << endl;
      write << "   {" << endl;
      write << "      total_mass += mass[k];" << endl;
      write << "   }" << endl;
      write << endl << endl;
      write << "   for (int k=0; k<nsp; k++)" << endl;
      write << "   {" << endl;
      write << "      C[k] = density*mass[k]/(mw[k]*total_mass);" << endl;
      write << "   }" << endl;
      write << "   " << endl;
      write << "   " << endl;
      write << "   for (int k=0; k<nsp; k++)" << endl;
      write << "   {" << endl;
      write << "      Y[k] = mass[k]/total_mass;" << endl;
      write << "   }" << endl;




   }
   

   write << endl;
   write << endl;

   int nSimple = 0;
   int nFalloff = 0;
   int nThreeBody = 0;
   int nTroe = 0;
   int nLindemann = 0;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <Simple *> (listReactions[j]))
         nSimple += 1;

      if (dynamic_cast <ThreeBody *> (listReactions[j]))
         nThreeBody += 1;

      if (dynamic_cast <FalloffR *> (listReactions[j]))
         nFalloff += 1;

      if (dynamic_cast <Troe *> (listReactions[j]))
         nTroe += 1;

      if (dynamic_cast <Lindemann *> (listReactions[j]))
         nLindemann += 1;
   }


   //write << "   optimisation_count += 1;" << endl;
   write << endl << endl;
   write << "   double CTB [" << nThreeBody << "];" << endl; //Falloff are seen as particular ThreeBody reactions

   write << endl << endl;
   write << "//Required inputs are : Temperature, Transported species concentration, ..." << endl;
   write << "//Output is : Transported species source terms " << endl;
   write << endl << endl;

   write << "//parameters" << endl;
   write << "   double RU = 8.31447215e3;" << endl;
   write << "   double PATM = 1.01325e5;" << endl;
   write << "   double SMALL = 1.0e-50;" << endl;
   write << endl << endl;

   write << "//variables " << endl;
   write << "   double Qf ["<< listReactions.size() << "];" << endl;
   write << "   double Qf_lowP ["<< nFalloff << "];" << endl;
   write << "   double Qr ["<< listReactions.size() << "];" << endl;

   write << endl;
   write << "   double SMH [" << nsp << "]; //Computed from NASA coefficients" << endl;
   write << "   double EG [" << nsp << "]; //Species equilibrium" << endl;

   write << endl;
   write << "   double EQ ["<< listReactions.size() << "]; //Reactions equilibrium" << endl;
   write << "   double ROP ["<< listReactions.size() << "]; //Reactions rate of progress" << endl;

   write << "   double WDOT [" << nsp << "]; //Transported species source terms" << endl;

   double R = 8.31447215;

   write << endl << endl;
   write << "   double lnT = log(T);" << endl;
   write << "   double TI = 1/T;" << endl;
   write << endl;


   write << endl << endl;

   write.precision(7);

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
    write << "   Qf[" << j << "] = ";
      //If non-zero activation energy
      if (listReactions[j]->m_E != 0.0)
      {
         write << "exp( " << log(listReactions[j]->m_A);
         if (listReactions[j]->m_E > 0.0)
            write << " -" << listReactions[j]->m_E*4.184/R << "*TI";
         else
            write << " +" << fabs(listReactions[j]->m_E*4.184/R) << "*TI";

         if (listReactions[j]->m_b != 0.0)
         {
            if (listReactions[j]->m_b > 0.0)
               write << " +" << listReactions[j]->m_b << "*lnT";
            else
               write << " -" << fabs(listReactions[j]->m_b) << "*lnT";
         }
         write << ");" << endl;
      }
      else
      {
         if (listReactions[j]->m_b != 0.0)
         {
            write << "exp( " <<  log(listReactions[j]->m_A);

            if (listReactions[j]->m_b > 0.0)
               write << " +" << listReactions[j]->m_b << "*lnT);" << endl;
            else
               write << " -" << fabs(listReactions[j]->m_b) << "*lnT);" << endl;
         }
         else
            write  << listReactions[j]->m_A << ";" << endl;
      }

   }

   write << endl << endl;


   int count_falloff = 0;
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <FalloffR *> (listReactions[j]))
      {
         write << "   Qf_lowP[" << count_falloff << "] = ";

         //If non-zero activation energy
         if (dynamic_cast <FalloffR *> (listReactions[j])->m_E_low != 0.0)
         {
            write << "exp(" << log(dynamic_cast <FalloffR *> (listReactions[j])->m_A_low);
            if (dynamic_cast <FalloffR *> (listReactions[j])->m_E_low > 0.0)
               write << " -" << dynamic_cast <FalloffR *> (listReactions[j])->m_E_low*4.184/R << "*TI";
            else
               write << " +" << fabs(dynamic_cast <FalloffR *> (listReactions[j])->m_E_low*4.184/R) << "*TI";

            if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low != 0.0)
            {
               if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low > 0.0)
                  write << " +" << dynamic_cast <FalloffR *> (listReactions[j])->m_b_low << "*lnT";
               else
                  write << " -" << fabs(dynamic_cast <FalloffR *> (listReactions[j])->m_b_low) << "*lnT";
            }
            write << ");" << endl;
         }
         else
         {
            if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low != 0.0)
            {
               write << "exp(" << log(dynamic_cast <FalloffR *> (listReactions[j])->m_A_low);

               if (dynamic_cast <FalloffR *> (listReactions[j])->m_b_low > 0.0)
                  write << " +" << dynamic_cast <FalloffR *> (listReactions[j])->m_b_low << "*lnT);" << endl;
               else
                  write << " -" << fabs(dynamic_cast <FalloffR *> (listReactions[j])->m_b_low) << "*lnT);" << endl;
            }
            else
               write << dynamic_cast <FalloffR *> (listReactions[j])->m_A_low << ";" << endl;
         }

         count_falloff += 1; 
      }

   }

   //Estimation of the equilibrium through the Gibbs free energy calculation
   write << endl;
   write << "   //Thermal data " << endl << endl;
   write << "   double Tn1, Tn2, Tn3, Tn4, Tn5; " << endl;
   write << "   Tn1 = lnT - 1.0; " << endl;
   write << "   Tn2 = T; " << endl;
   write << "   Tn3 = Tn2 * T; " << endl;
   write << "   Tn4 = Tn3 * T; " << endl;
   write << "   Tn5 = Tn4 * T; " << endl;

   write << endl << endl;

   //Check for the different Temperatures Tmax in the low T region
   vector<double> T;
   bool Check_new_T = true;

   T.push_back(listSpecies[0]->m_NASA_lowT_max);
   int nVarious_T = 1;

   for (int k=0; k<nsp; k++)
   {
      Check_new_T = true;
      for (int n=0; n<nVarious_T; n++)
      {
         if (listSpecies[k]->m_NASA_lowT_max == T[n])
         {
            Check_new_T = false;
         }
      }

      if (Check_new_T == true)
      {
         T.push_back(listSpecies[k]->m_NASA_lowT_max);
         nVarious_T += 1;
      }
   }

   write << endl << endl;
   write << "//Low temperature Gibbs free energy" << endl;
   //Low temperature Gibbs values
   for (int n=0; n<nVarious_T; n++)
   {
      write << "   if (T < " << T[n] << ") " << endl;
      write << "   { " << endl;

      for (int k=0; k<nsp; k++)
      {
         if (listSpecies[k]->m_NASA_lowT_max == T[n])
         {

            write << "      SMH[" << k << "] = " << print_sign(listSpecies[k]->m_NASACoeffs_lowT[0]) << "*Tn1 "
                                              << print_sign(0.5*listSpecies[k]->m_NASACoeffs_lowT[1]) << "*Tn2 " << endl
                                              << "                " << print_sign(listSpecies[k]->m_NASACoeffs_lowT[2]/6) << "*Tn3 "
                                              << print_sign(listSpecies[k]->m_NASACoeffs_lowT[3]/12) << "*Tn4 " << endl
                                              << "                " << print_sign(0.05*listSpecies[k]->m_NASACoeffs_lowT[4]) << "*Tn5 "
                                              << print_sign(-listSpecies[k]->m_NASACoeffs_lowT[5]) << "*TI " << endl
                                              << "                " << print_sign(listSpecies[k]->m_NASACoeffs_lowT[6])
                                              << ";"
                                              << endl;
         }
      }
      write << "   } " << endl;
   }

   write << endl << endl;
   write << "//High temperature Gibbs free energy" << endl;
   //High temperature Gibbs values
   for (int n=0; n<nVarious_T; n++)
   {
      write << "   if (T > " << T[n] << ") " << endl;
      write << "   { " << endl;

      for (int k=0; k<nsp; k++)
      {
         if (listSpecies[k]->m_NASA_lowT_max == T[n])
         {

            write << "      SMH[" << k << "] = " << print_sign(listSpecies[k]->m_NASACoeffs_highT[0]) << "*Tn1 "
                                              << print_sign(0.5*listSpecies[k]->m_NASACoeffs_highT[1]) << "*Tn2 " << endl
                                              << "                " << print_sign(listSpecies[k]->m_NASACoeffs_highT[2]/6) << "*Tn3 "
                                              << print_sign(listSpecies[k]->m_NASACoeffs_highT[3]/12) << "*Tn4 " << endl
                                              << "                " << print_sign(0.05*listSpecies[k]->m_NASACoeffs_highT[4]) << "*Tn5 "
                                              << print_sign(-listSpecies[k]->m_NASACoeffs_highT[5]) << "*TI " << endl
                                              << "                " << print_sign(listSpecies[k]->m_NASACoeffs_highT[6])
                                              << ";"
                                              << endl;
         }
      }
      write << "   } " << endl;
   }

   write << endl << endl;

   write << "   for (int k=0; k<" << nsp << "; k++) " << endl;
   write << "   { " << endl;
   write << "      EG[k] = exp(SMH[k]); " << endl;
   write << "   } " << endl;

   write << endl;
   write << "   double PFAC = PATM / (RU*T); " << endl;
   write << "   double PFAC2 = PFAC*PFAC; " << endl;
   write << "   double PFAC3 = PFAC2*PFAC; " << endl;

   bool Add_multiplication;


   //Equilibrium on the SimpleReactions
   write << endl << endl;
   write << "//Equilibrium constants " << endl;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      //If the reaction is not reversible, it is not necessary to estimate the equilibrium
      if (listReactions[j]->m_reversible && Keep_rev_reaction[j])
      {
         write << "   EQ[" << j <<"] = ";

         Add_multiplication = false;
         for (unsigned int nj=0; nj<listReactions[j]->m_ProductSpecies.size(); nj++)
         {
            if (Get_species_number(listReactions[j]->m_ProductSpecies[nj], listSpecies) != -1)
            {
               if (Add_multiplication)
               {
                  write << "*";
               }
               write << "EG[" << Get_species_number(listReactions[j]->m_ProductSpecies[nj], listSpecies) << "]";

               if (listReactions[j]->m_ProductStoichCoeffs[nj] == 2)
               {
                   write << "*EG[" << Get_species_number(listReactions[j]->m_ProductSpecies[nj], listSpecies) << "]";
               }

               Add_multiplication = true;
            }
         }
         for (unsigned int nj=0; nj<listReactions[j]->m_ReactantSpecies.size(); nj++)
         {
            if (Get_species_number(listReactions[j]->m_ReactantSpecies[nj], listSpecies) != -1)
            {
               write << "/EG[" << Get_species_number(listReactions[j]->m_ReactantSpecies[nj], listSpecies) << "]";

               if (listReactions[j]->m_ReactantStoichCoeffs[nj] == 2)
               {
                   write << "/EG[" << Get_species_number(listReactions[j]->m_ReactantSpecies[nj], listSpecies) << "]";
               }
            }
         }

         double Delta_n = 0.0;
         for (unsigned int k=0; k<listReactions[j]->m_ReactantSpecies.size(); k++)
         {
            Delta_n -= listReactions[j]->m_ReactantStoichCoeffs[k];
         }
         for (unsigned int k=0; k<listReactions[j]->m_ProductSpecies.size(); k++)
         {
            Delta_n += listReactions[j]->m_ProductStoichCoeffs[k];
         }

         if (Delta_n == 1)
         {
            write << "*PFAC";
         }
         else if (Delta_n == -1)
         {
            write << "/PFAC";
         }
         else if (Delta_n != 0)
         {
            cout << "ERROR - Delta n is different from 1 or -1" << "  value: " << Delta_n << endl;
         }

         write << ";" << endl;
      }

   }


   write << "   for (int j=0; j<" << listReactions.size() << "; j++)" << endl;
   write << "   {" << endl;
   write << "      Qr[j] = 0.0;" << endl;
   write << "   }" << endl;



   write << endl << endl << endl;
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (listReactions[j]->m_reversible && Keep_rev_reaction[j])
      {
         write << "   Qr[" << j << "] = Qf[" << j << "] / max(EQ[" << j << "], SMALL);" << endl;
      }
      else
         write << "   //Reaction " << j << " is not reversible" << endl;
   }


   write.precision(5);

write << endl << endl << "/*" << endl;
write << "//------------------------------------------------------------" << endl;
   write << "   cout << \"Before multiplying by TB, Falloff and conc \" << endl;" << endl;

   write << "   for (int j=0; j<" << listReactions.size() << "; j++)" << endl;
   write << "   {" << endl;
   write << "      cout << j << \"  \" << Qf[j] << endl; " << endl;
   write << "      cout << j << \"  \" << Qr[j] << endl; " << endl;
   write << "      cout << endl; " << endl;
   write << "   }" << endl;

   write << endl << "   // getchar(); " << endl;
write << "//------------------------------------------------------------" << endl;
write << "*/" << endl << endl;



   write << endl << endl << endl;
   write << "//Three Body concentration " << endl << endl;

   write << "   double Ctot = 0.0;" << endl;
   write << "   for (int k=0; k<nsp; k++)" << endl;
   write << "   { " << endl;
   write << "      Ctot += C[k]; " << endl;
   write << "   } " << endl << endl;

   write << endl << endl;
   write << "//Three body efficiencies" << endl;


   int count_threebody = 0;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <ThreeBody *> (listReactions[j]))
      {
         write << "   CTB[" << count_threebody << "] = Ctot";
         for (unsigned int k=0; k<(dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc.size(); k++)
         {

            if (Index_QSS_Species_number
                  [Get_species_number((dynamic_cast <ThreeBody *> (listReactions[j]))->
                                                  m_NameTBconc[k], listSpecies)] != -1)
            {
               if ((dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc[k]-1 > 0.0)
                  write << " +" << (dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc[k]-1;
               else
                  write << " -" << fabs((dynamic_cast <ThreeBody *> (listReactions[j]))->m_TBconc[k]-1);
               write << " *C[";
               write << Index_QSS_Species_number[Get_species_number((dynamic_cast <ThreeBody *> (listReactions[j]))->m_NameTBconc[k], listSpecies)];
               write << "]";
            }

         }
         write << ";" << endl;
         count_threebody += 1;
      }
   }

   write << endl << endl;


   write << endl << endl << endl;
   write << "   double Pr;" << endl;
   write << "   double Pcor;" << endl;
   write << "   double PrLog;" << endl;
   write << "   double Fcent;" << endl;
   write << "   double FcLog;" << endl;
   write << "   double Xn;" << endl;
   write << "   double CprLog;" << endl;
   write << "   double FLog;" << endl;
   write << "   double Fc;" << endl;


   count_falloff = 0;
   count_threebody = 0;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (dynamic_cast <ThreeBody *> (listReactions[j]))
      {
         if (dynamic_cast <Lindemann *> (listReactions[j]))
         {
            write << endl;
            write << "//Lindemann type" << endl;
            write << "   Pr = Qf_lowP[" << count_falloff << "] * CTB[" << count_threebody <<"]" 
                              << " / Qf[" << j << "];" << endl;
            write << "   Pcor = Pr / (1.0 + Pr);" << endl;
            if (Keep_fwd_reaction[j])
               write << "   Qf[" << j << "] = Qf[" << j << "] * Pcor;" << endl;
            if (listReactions[j]->m_reversible && Keep_rev_reaction[j])
               write << "   Qr[" << j << "] = Qr[" << j << "] * Pcor;" << endl;

            count_falloff += 1;
         }

         else if (dynamic_cast <Troe *> (listReactions[j]))
         {
            write << endl;
            write << "//Troe type" << endl;
            write << "   Pr = Qf_lowP[" << count_falloff << "] * CTB[" << count_threebody << "]" 
                              << " / Qf[" << j << "];" << endl;
            write << "   Pcor = Pr / (1.0 + Pr);" << endl;
            write << "   PrLog = log10(max(Pr, SMALL));" << endl;
            write << "   Fcent = " << (1-(dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[0]) 
                         << "*exp(-T/" << (dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[1] << ")+" 
                         << (dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[0];
            write << "*exp(-T/" << (dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[2] << ")" << endl;
            write << "           +exp(" << -(dynamic_cast <Troe *> (listReactions[j]))->m_TroeCoeffs[3] << "/T);" << endl;
            write << "   FcLog = log10(max(Fcent, SMALL));" << endl;
            write << "   Xn = 0.75 - 1.27*FcLog;" << endl;
            write << "   CprLog = PrLog - (0.4 + 0.67*FcLog);" << endl;
            write << "   FLog = FcLog/(1.0 + pow(CprLog/(Xn-0.14*CprLog),2));" << endl;
            write << "   Fc = pow(10, FLog);" << endl;
            write << "   Pcor *= Fc;" << endl;
            if (Keep_fwd_reaction[j])
               write << "   Qf[" << j << "] = Qf[" << j << "] * Pcor;" << endl;
            if (listReactions[j]->m_reversible && Keep_rev_reaction[j])
               write << "   Qr[" << j << "] = Qr[" << j << "] * Pcor;" << endl;

            count_falloff += 1;
         }
         count_threebody += 1;
      }
   }

   write << endl << endl << endl;


//Multiply by the concentration of the known species (not the one of the QSS species)

   count_threebody = 0;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (Keep_fwd_reaction[j])
      {
         write << "   Qf[" << j << "] = Qf[" << j << "]";

         for (unsigned int k=0; k<listReactions[j]->m_ReactantSpecies.size(); k++)
         {
            if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)] != true)
            {
                 write << "*C[" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)] << "]";

                 if (listReactions[j]->m_ReactantStoichCoeffs[k] != 1.0)
                 {
                   if (listReactions[j]->m_ReactantStoichCoeffs[k] == 2.0)
                   {
                      write << "*C[" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)] << "]";
                   }
                   else
                   {
                      cout << "Value of the StoichCoeff " << listReactions[j]->m_ReactantStoichCoeffs[k] << endl;
                      cout << "Error" << endl;
                   }
                 }
             } //end if QSS species
         }

         if ((dynamic_cast <ThreeBody *> (listReactions[j])))
         {
            if (!(dynamic_cast <FalloffR *> (listReactions[j])))
            {
               write << " *CTB[" << count_threebody << "]";
            }
         }


         write << ";" << endl;
      }
      if ((dynamic_cast <ThreeBody *> (listReactions[j])))
      {
         count_threebody += 1;
      }
   }

   write << endl << endl;

   count_threebody = 0;

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (listReactions[j]->m_reversible && Keep_rev_reaction[j])
      {
         write << "   Qr[" << j << "] = Qr[" << j << "]";
   
         for (unsigned int k=0; k<listReactions[j]->m_ProductSpecies.size(); k++)
         {
            if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)] != true)
            {
                 write << "*C[" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)] << "]";
   
                 if (listReactions[j]->m_ProductStoichCoeffs[k] != 1.0)
                 {
                   if (listReactions[j]->m_ProductStoichCoeffs[k] == 2.0)
                   {
                      write << "*C[" << Index_QSS_Species_number[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)] << "]";
                   }
                   else
                   {
                      cout << "Value of the StoichCoeff " << listReactions[j]->m_ProductStoichCoeffs[k] << endl;
                      cout << "Error" << endl;
                   }
                 }
             } //end if QSS species
         }

         if ((dynamic_cast <ThreeBody *> (listReactions[j])))
         {
            if (!(dynamic_cast <FalloffR *> (listReactions[j])))
            {
               write << " *CTB[" << count_threebody << "]";
            }
         }
         write << ";" << endl;
      }
      if ((dynamic_cast <ThreeBody *> (listReactions[j])))
      {
         count_threebody += 1;
      }
   }

write << endl << endl << "/*" << endl;
write << "//------------------------------------------------------------" << endl;
   write << "   cout << \"Just before creating the matrix \" << endl;" << endl;

   write << "   for (int j=0; j<" << listReactions.size() << "; j++)" << endl;
   write << "   {" << endl;
   write << "      cout << j << \"  \" << Qf[j] << endl; " << endl;
   write << "      cout << j << \"  \" << Qr[j] << endl; " << endl;
   write << "   }" << endl;
   write << "   cout << endl; " << endl;

   write << endl << "   // getchar(); " << endl;

write << "//------------------------------------------------------------" << endl;
write << "*/" << endl << endl;





   write << endl << endl << endl;

   vector<bool> Reaction_fwd_to_zero (listReactions.size(), false);
   vector<bool> Reaction_rev_to_zero (listReactions.size(), false);

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (Keep_fwd_reaction[j])
      {
         for (unsigned int k=0; k<listReactions[j]->m_ReactantSpecies.size(); k++)
         {
            if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[k], listSpecies)])
            {
               if (listReactions[j]->m_ReactantStoichCoeffs[k] == 2)
               {
                  Reaction_fwd_to_zero[j] = true;
               }
   
               for (unsigned int ksec=0; ksec<listReactions[j]->m_ReactantSpecies.size(); ksec++)
               {
                  if (ksec != k)
                  {
                     if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[ksec], listSpecies)])
                     {
                         Reaction_fwd_to_zero[j] = true;
                     }
                  }
               } //end for ksecond
            }
         }
      }

      if (Keep_rev_reaction[j])
      {
         for (unsigned int k=0; k<listReactions[j]->m_ProductSpecies.size(); k++)
         {
            if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[k], listSpecies)])
            {
               if (listReactions[j]->m_ProductStoichCoeffs[k] == 2)
               {
                  Reaction_rev_to_zero[j] = true;
               }

               for (unsigned int ksec=0; ksec<listReactions[j]->m_ProductSpecies.size(); ksec++)
               {
                  if (ksec != k)
                  {
                     if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[ksec], listSpecies)])
                     {
                         Reaction_rev_to_zero[j] = true;
                     }
                  }
               } //end for ksecond
            }
         }
      }
   }




   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      if (Reaction_fwd_to_zero[j])
         write << "   Qf[" << j << "] = 0.0;" << endl;
      if (Reaction_rev_to_zero[j])
         write << "   Qr[" << j << "] = 0.0;" << endl;
   }

   write << endl << endl << endl;


   //Solve for the QSS species
   write << endl << endl << endl;
   write << "//Solve for the QSS species " << endl;
 
   // Commented by Huu-Tri Nguyen - 2020.02.17
//HT   int nQSS_Species = 0;
//HT   for (int k=0; k<nsp; k++)
//HT   {
//HT      if (QSS_Species[k])
//HT      {
//HT         nQSS_Species += 1;
//HT      }
//HT   }
   // End comment

if (nQSS_Species > 0)
{
   vector<vector<bool> > QSS_matrix (nQSS_Species, vector<bool> (nQSS_Species, false));

   for (int k=0; k<nQSS_Species; k++)
   {
         for (int kb=0; kb<nQSS_Species; kb++)
         {
            if (k == kb)
               QSS_matrix[k][kb] = true;
            else
               QSS_matrix[k][kb] = false;
         }
   }


    for (int i=0; i<nQSS_Species; i++)
    {
       write << "   double x" << i << " = 0.0;" << endl;
    }
    write << endl;

    for (int i=0; i<nQSS_Species; i++)
    {
       write << "   double RHS" << i << " = 0.0;" << endl;
    }
    write << endl;

    for (int i=0; i<nQSS_Species; i++)
    {
       for (int j=0; j<nQSS_Species; j++)
       {
          write << "   double A" << i << "_" << j << " = 0.0;" << endl;
       }
    }

   write << "   double Den;" << endl;


   vector <int> QSSindex (nsp, -1);

   int i = 0;
   for (int k=0; k<nsp; k++)
   {
      if (QSS_Species[k])
      {
         QSSindex[k] = i;
         i += 1;
      }
   }

   string Den;
   int count_den;
   int count_RHS;
   bool look_on_the_other_side;
   bool write_into_RHS;

   vector<string> QSS_SpeciesCoeffs (nQSS_Species+1, "");

   for (int k=0; k<nsp; k++)
   {
      for (int ka=0; ka<nQSS_Species+1; ka++)
         QSS_SpeciesCoeffs[ka] = "";

      Den = "";

      if (QSS_Species[k])
      {
         count_den = 1;
         count_RHS = 1;

         for (unsigned int j=0; j<listReactions.size(); j++)
         {
            look_on_the_other_side = false;

            if (Keep_fwd_reaction[j])
            {
               for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
               {
                  if (k == Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies))
                  {
                     stringstream sstm;
                     sstm << " +Qf[" << j << "]";
                     if (count_den == 5)
                     {
                        sstm << endl << "       ";
                        count_den = 0;
                     }
                     Den.append(sstm.str());
                     count_den += 1;
                     look_on_the_other_side = true;
                  }
               }
            }


            if (look_on_the_other_side)
            {
               write_into_RHS = true;
               //Look for the other QSS species involved (which can only be present on the other side of the reaction)

               if (Keep_rev_reaction[j])
               {
                  for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
                  {
                     if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies)])
                     {
                        stringstream sstm;
                        sstm << " +Qr[" << j << "]";

                        QSS_SpeciesCoeffs[QSSindex[Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies)]]
                                         .append(sstm.str());

                        QSS_matrix[QSSindex[k]]
                                  [QSSindex[Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies)]] = true;
                        write_into_RHS = false;

                     }
                  }

                  if (write_into_RHS == true)
                  {
                         stringstream sstm;
                         sstm << " +Qr[" << j << "]";
                         if (count_RHS == 5)
                         {
                            sstm << endl << "       ";
                            count_RHS = 0;
                         }
                         QSS_SpeciesCoeffs[nQSS_Species].append(sstm.str());
                         count_RHS += 1;
                  }
               }
            } //end if (look_on_the_other_side)


      //---Repeat this for a QSS species into the products

            look_on_the_other_side = false;
            if (Keep_rev_reaction[j])
            {
               for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
               {
                  if (k == Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies))
                  {
                     stringstream sstm;
                     sstm << " +Qr[" << j << "]";
                     if (count_den == 5)
                     {
                        sstm << endl << "       ";
                        count_den = 0;
                     }
                     Den.append(sstm.str());
                     count_den += 1;
                     look_on_the_other_side = true;
                  }
               }
            }


            if (look_on_the_other_side)
            {
               write_into_RHS = true;
               //Look for the other QSS species involved (which can only be present on the other side of the reaction)

               if (Keep_fwd_reaction[j])
               {
                  for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
                  {
                     if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies)])
                     {
                        stringstream sstm;
                        sstm << " +Qf[" << j << "]";

                        QSS_SpeciesCoeffs[QSSindex[Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies)]]
                                         .append(sstm.str());

                        QSS_matrix[QSSindex[k]]
                                  [QSSindex[Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies)]] = true;
                        write_into_RHS = false;
                     }
                  }
                  if (write_into_RHS == true)
                  {
                         stringstream sstm;
                         sstm << " +Qf[" << j << "]";
                         if (count_RHS == 5)
                         {
                            sstm << endl << "       ";
                            count_RHS = 0;
                         }
                         QSS_SpeciesCoeffs[nQSS_Species].append(sstm.str());
                         count_RHS += 1;
                  }
               }
            } //end if (look_on_the_other_side)
         } //end for (unsigned int j=0; j<listReactions.size(); j++)


         size_t remove_first_plus = Den.find("+");
         Den.erase(0, remove_first_plus+1);

         write << endl << endl << "//QSS species " << listSpecies[k]->m_Name << endl;

         if (Den != "")
         {
            write << "   Den = " << Den << ";" << endl;

write << "//------------------------------------------------------------" << endl;
    write << "  // cout << \"Den \" << Den << endl;" << endl; 
write << "//------------------------------------------------------------" << endl << endl;

         }
         else
            //should be removed to something else
            write << "   Den = 0.0;" << endl;
           

         for (int ka=0; ka<nQSS_Species; ka++)
         {
            if (QSS_SpeciesCoeffs[ka] != "")
            {
               remove_first_plus = QSS_SpeciesCoeffs[ka].find("+");
               QSS_SpeciesCoeffs[ka].erase(0, remove_first_plus+1);
               write << "   A" << QSSindex[k] << "_" << ka << " = ";
               write << "(" << QSS_SpeciesCoeffs[ka] << ") /max(Den, SMALL);" << endl;

write << "//------------------------------------------------------------" << endl;
    write << "  // cout << \"A" << QSSindex[k] << "_" << ka << " \" << A" << QSSindex[k] << "_" << ka << " << endl;" << endl; 
write << "//------------------------------------------------------------" << endl << endl;
            }

         }

         remove_first_plus = QSS_SpeciesCoeffs[nQSS_Species].find("+");
         QSS_SpeciesCoeffs[nQSS_Species].erase(0, remove_first_plus+1);

         if (QSS_SpeciesCoeffs[nQSS_Species] != "")
         {
            write << "   RHS" << QSSindex[k] << " = ";
            write << "-(" << QSS_SpeciesCoeffs[nQSS_Species] << ") /max(Den, SMALL);" << endl;

write << "//------------------------------------------------------------" << endl;
    write << "  // cout << \"RHS" << QSSindex[k] << " \" << RHS" << QSSindex[k] << " << endl;" << endl; 
write << "//------------------------------------------------------------" << endl << endl;
         }
         else
         {
            write << "   RHS" << QSSindex[k] << " = 0.0;";
         }

      } //end if QSS_Species[k]
   } //end for k=0; k<nsp; k++

   write << endl << endl;

   for (int k=0; k<nQSS_Species; k++)
   {
      write << "   A" << k << "_" << k << " = -1;" << endl;
   }

//---------------------------------------------------------------------
//   "Solve the matrix"
//---------------------------------------------------------------------


   write << endl << endl << endl;
   write << "//Solve the matrix" << endl;

   write << endl << endl;
   write << "   double coefficient;" << endl;

   vector<double> b (nQSS_Species);
   vector<vector<double> > A (nQSS_Species,vector<double> (nQSS_Species));

   for (int i=0; i<nQSS_Species; i++)
   {
      b[i] = 2.0*(i+1); //Give ordinary values to the b vector
      for (int j=0; j<nQSS_Species; j++)
      {
         A[i][j] = QSS_matrix[i][j];
      }
   }

   vector<double> x (nQSS_Species);

   //Elimination
   write << endl << "//-----Elimination-----" << endl;

    for(int i=0; i<nQSS_Species-1; i++)
    {
        for(int j=i+1; j<nQSS_Species; j++)
        {
            double Coeff = -(A[j][i]/A[i][i]);
            bool storeCoeff = true;

            for(int k = i; k < nQSS_Species; k++)
            {
                A[j][k] = A[j][k] + Coeff*A[i][k];

                if (Coeff != 0.0)
                {
                   if (A[i][k] != 0.0)
                   {
                      if (storeCoeff)
                      {
                         write << "   coefficient = -(A" << j << "_" << i << "/A" << i << "_" << i << ");" << endl;
                         storeCoeff = false;
                      }
                      write << "   A" << j << "_" << k << " += coefficient*A" << i << "_" << k << ";" << endl;
                   }
                }
            }

            b[j] = b[j] + Coeff*b[i];

            if (storeCoeff == false)
               write << "   RHS" << j << " += coefficient*RHS" << i << ";" << endl;
        }
    }




    //Back-substitution
    write << endl << "//-----Back substitution-----" << endl;

    x[nQSS_Species-1] = b[nQSS_Species-1]/A[nQSS_Species-1][nQSS_Species-1];
    write << "   x" << nQSS_Species-1 << " = " 
          << "RHS" << nQSS_Species-1 << " /A" << nQSS_Species-1 << "_" << nQSS_Species-1 << ";" << endl;

write << "//------------------------------------------------------------" << endl;
    write << "   // cout << \"x" << nQSS_Species-1 << " \" << " << "x" << nQSS_Species-1 << " << endl;" << endl; 
write << "//------------------------------------------------------------" << endl << endl;


    write << endl << "   double Sum;" << endl;
    for(int i = nQSS_Species-2; i>=0; i--)
    {
        stringstream sstm;
        sstm << "   Sum = ";
        bool Coeff = false;

        double sum = 0;
        for(int j = i; j < nQSS_Species; j++)
        {
            sum += x[j]*A[i][j];

           if ((x[j]*A[i][j]) != 0.0)
           {
              sstm << "+ x" << j << "*A" << i << "_" << j << "";
              Coeff = true;
           }

        }
        sstm << ";" << endl;

        if (Coeff)
        {
          write << sstm.str();
          write << "   x" << i << " = (RHS" << i << "-Sum)/A" << i << "_" << i << ";" << endl << endl << endl;
        }
        else
        {
          write << "   x" << i << " = RHS" << i << "/A" << i << "_" << i << ";" << endl << endl << endl;
        }

write << "//------------------------------------------------------------" << endl;
    write << "  // cout << \"x" << i << " \" << " << "x" << i << " << endl;" << endl; 
write << "//------------------------------------------------------------" << endl << endl;

        x[i] = (b[i]-sum)/A[i][i];
    }

   write << " // getchar(); " << endl;



//   End of "Solve the matrix"
//---------------------------------------------------------------------


   write << endl << endl << endl << endl;
   write << "   //Multiply the forward rates by the concentration of the QSS species " << endl;
   write << endl << endl;
   for (unsigned int j=0; j<listReactions.size(); j++)
   {
       for (int k=0; k<nsp; k++)
       {
          if (QSS_Species[k])
          {
            if (Keep_fwd_reaction[j])
            {
               for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
               {
                  if (k == Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies))
                  {
                      //Within the forward part of this Simple reaction there is a QSS species
                      write << "   Qf[" << j << "] *= x" << QSSindex[k] << ";" << endl;
                  }
               }
            }

            if (Keep_rev_reaction[j])
            {
               for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
               {
                  if (k == Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies))
                  {
                      //Within the forward part of this Simple reaction there is a QSS species
                      write << "   Qr[" << j << "] *= x" << QSSindex[k] << ";" << endl;
                  }
               }
            }
          }
       }
   }


} //End if (nQSS_Species > 0) 
else
{
   write << endl << "//---0 QSS expressions---//" << endl;
}


  write << endl << endl << endl;
 // write << "   for (int j=0; j<" << listReactions.size() << "; j++)" << endl;

  for (unsigned int j=0; j<listReactions.size(); j++)
  {
     if (Keep_fwd_reaction[j] && Keep_rev_reaction[j])
        write << "      ROP[" << j << "] = Qf[" << j << "] - Qr[" << j << "];" << endl;
     if (Keep_fwd_reaction[j] && Keep_rev_reaction[j] == false)
        write << "      ROP[" << j << "] = Qf[" << j << "];" << endl;
     if (Keep_fwd_reaction[j] == false && Keep_rev_reaction[j])
        write << "      ROP[" << j << "] = -Qr[" << j << "];" << endl;
  }



  write << endl << endl << endl;
  for (int k=0; k<nsp; k++)
  {
     if (QSS_Species[k] == false)
     {
        write << "//Species " << Index_QSS_Species_number[k] << endl;
        write << "   WDOT[" << Index_QSS_Species_number[k] << "] =" << endl;

        bool term = false;

        for (unsigned int j=0; j<listReactions.size(); j++)
        {
           for (unsigned int ka=0; ka<listReactions[j]->m_ReactantSpecies.size(); ka++)
           {
              if (k == Get_species_number(listReactions[j]->m_ReactantSpecies[ka], listSpecies))
              {
//                   write << "-" << listReactions[j]->m_ReactantStoichCoeffs[ka]  << "*ROP[" << j << "] ";
//                   term = true;
                     
                 if (listReactions[j]->m_ReactantStoichCoeffs[ka] == 2)
                 {
                    write << "-2*ROP[" << j << "] ";
                    term = true;
                 }
                 else if (listReactions[j]->m_ReactantStoichCoeffs[ka] == 1)
                 {
                    write << "-ROP[" << j << "] ";
                    term = true;
                 }
                 else if (listReactions[j]->m_ReactantStoichCoeffs[ka] == 3)
                 {
                    write << "-3*ROP[" << j << "] ";
                    term = true;
                 }
                 else
                 {
                    cout << "Warning stoichiometric value different from 1 or 2 or 3 for reaction " << j << endl;
                   write << "-"  << "0*ROP[" << j << "] ";
                   //write << "-" << listReactions[j]->m_ReactantStoichCoeffs[ka]  << "*ROP[" << j << "] ";
                   term = true;
                 }

              }
           }

           for (unsigned int ka=0; ka<listReactions[j]->m_ProductSpecies.size(); ka++)
           {
              if (k == Get_species_number(listReactions[j]->m_ProductSpecies[ka], listSpecies))
              {
//                   write << "+" << listReactions[j]->m_ReactantStoichCoeffs[ka]  << "*ROP[" << j << "] ";
  //                 term = true;

                 if (listReactions[j]->m_ProductStoichCoeffs[ka] == 2)
                 {
                    write << "+2*ROP[" << j << "] ";
                    term = true;
                 }
                 else if (listReactions[j]->m_ProductStoichCoeffs[ka] == 1)
                 {
                    write << "+ROP[" << j << "] ";
                    term = true;
                 }
                 else if (listReactions[j]->m_ProductStoichCoeffs[ka] == 3)
                 {
                    write << "+3*ROP[" << j << "] ";
                    term = true;
                 }
                 else
                 {
                    cout << "Warning stoichiometric value different from 1 or 2 or 3 for reaction " << j << endl;
                   write << "+" << listReactions[j]->m_ReactantStoichCoeffs[ka]  << "*ROP[" << j << "] ";
                   term = true;
                 }
              }
           }
        }

        if (term == false)
           write << "0.0" << endl;


        write << ";" << endl << endl;
     }
  }


   write << endl << endl;
   write << "   for (int k=0; k<nsp; k++)" << endl;
   write << "   {" << endl;
   write << "       w[k] = WDOT[k];" << endl;
   write << "   }" << endl;
   write << endl;

   write << "}" << endl;

   write.close();



}



void Write_QSS::Check_Non_Linearity(string mech, vector<vector<double> > QSS_Criteria, int rank) const
{

   ofstream log("computeQSSCriteria.log", ios::app);

   vector<Species_ORCh*> listSpecies;
   vector<Reaction_ORCh*> listReactions;

   Read *r = new Read();
   r->Read_species(mech, listSpecies);
   r->Read_reactions(mech, listReactions);

   int nsp = listSpecies.size();

   vector<vector<int> > interactions (nsp, vector<int> (nsp, 0));

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
      for (unsigned int kb=0; kb<listReactions[j]->m_ReactantSpecies.size(); kb++)
      {
         if (listReactions[j]->m_ReactantStoichCoeffs[kb] == 2)
         {
            interactions[Get_species_number(listReactions[j]->m_ReactantSpecies[kb], listSpecies)][Get_species_number(listReactions[j]->m_ReactantSpecies[kb], listSpecies)] += 1;
         }

         for (unsigned int ks=0; ks<listReactions[j]->m_ReactantSpecies.size(); ks++)
         {
            if (ks != kb)
               interactions[Get_species_number(listReactions[j]->m_ReactantSpecies[kb], listSpecies)][Get_species_number(listReactions[j]->m_ReactantSpecies[ks], listSpecies)] += 1;
         }

      }

      if (listReactions[j]->m_reversible)
      {
         for (unsigned int kb=0; kb<listReactions[j]->m_ProductSpecies.size(); kb++)
         {
            if (listReactions[j]->m_ProductStoichCoeffs[kb] == 2)
            {
               interactions[Get_species_number(listReactions[j]->m_ProductSpecies[kb], listSpecies)][Get_species_number(listReactions[j]->m_ProductSpecies[kb], listSpecies)] += 1;
            }

            for (unsigned int ks=0; ks<listReactions[j]->m_ProductSpecies.size(); ks++)
            {
               if (ks != kb)
                  interactions[Get_species_number(listReactions[j]->m_ProductSpecies[kb], listSpecies)][Get_species_number(listReactions[j]->m_ProductSpecies[ks], listSpecies)] += 1;
            }
         }
      }
   } 

   if (rank==0)
   {
      log << endl << endl;
      for (int k=0; k<nsp; k++)
      {
         if (QSS_Criteria[k][0] < 0.2)
         {
            cout << endl;
            cout << "Interactions with species " << listSpecies[int(QSS_Criteria[k][1])]->m_Name << "  with QSS Criteria " << QSS_Criteria[k][0] << endl; 
            log <<  "Interactions with species " << listSpecies[int(QSS_Criteria[k][1])]->m_Name << "  with QSS Criteria " << QSS_Criteria[k][0] << endl; 
            for (int ks=0; ks<nsp; ks++)
            {
               if ((interactions[QSS_Criteria[k][1]][ks] > 0) && (QSS_Criteria[ks][0] < 0.2))
               {
                  cout << listSpecies[ks]->m_Name << ":" << interactions[QSS_Criteria[k][1]][ks] << "  ";
                  log << listSpecies[ks]->m_Name << ":" << interactions[QSS_Criteria[k][1]][ks] << "  ";
               }
            }
            cout << endl;
            log << endl  <<  endl;
         }   
      }
      
      log.close();
   }


}



/*
void Write_QSS::Check_Non_Linearity(string mech, vector<bool> &QSS_Species) const
{

   vector<Species*> listSpecies;
   vector<Reaction_ORCh*> listReactions;

   Read *r = new Read();
   r->Read_species(mech, listSpecies);
   r->Read_reactions(mech, listReactions);

   int nsp = listSpecies.size();

   int total_counter;
   bool loop = true;


   while(loop)
   {
      total_counter = 0;

      vector<vector<int> > Number_of_QSS_interactions (nsp, vector<int> (nsp, 0));
      vector<int> Number_of_QSS_stoichiometric_coefficient_problem (nsp, 0);

   for (unsigned int j=0; j<listReactions.size(); j++)
   {
         for (unsigned int kbis=0; kbis<listReactions[j]->m_ReactantSpecies.size(); kbis++)
         {
            if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[kbis], listSpecies)])
            {
               if (listReactions[j]->m_ReactantStoichCoeffs[kbis] == 2)
               {
                  Number_of_QSS_stoichiometric_coefficient_problem[
                                Get_species_number(listReactions[j]->m_ReactantSpecies[kbis], listSpecies)] += 1;
               }

               for (unsigned int ksec=0; ksec<listReactions[j]->m_ReactantSpecies.size(); ksec++)
               {
                  if (ksec != kbis)
                  {
                     if (QSS_Species[Get_species_number(listReactions[j]->m_ReactantSpecies[ksec], listSpecies)])
                     {
                         Number_of_QSS_interactions[Get_species_number(listReactions[j]->m_ReactantSpecies[kbis], listSpecies)]
                                                   [Get_species_number(listReactions[j]->m_ReactantSpecies[ksec], listSpecies)] += 1;
                     }
                  }
               }
            } //if it is one QSS species
         } //loop on the reactants


         for (unsigned int kbis=0; kbis<listReactions[j]->m_ProductSpecies.size(); kbis++)
         {
            if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[kbis], listSpecies)])
            {
               if (listReactions[j]->m_ProductStoichCoeffs[kbis] == 2)
               {
                  Number_of_QSS_stoichiometric_coefficient_problem[
                                Get_species_number(listReactions[j]->m_ProductSpecies[kbis], listSpecies)] += 1;
               }

               for (unsigned int ksec=0; ksec<listReactions[j]->m_ProductSpecies.size(); ksec++)
               {
                  if (ksec != kbis)
                  {
                     if (QSS_Species[Get_species_number(listReactions[j]->m_ProductSpecies[ksec], listSpecies)])
                     {
                         Number_of_QSS_interactions[Get_species_number(listReactions[j]->m_ProductSpecies[kbis], listSpecies)]
                                                   [Get_species_number(listReactions[j]->m_ProductSpecies[ksec], listSpecies)] += 1;
                     }
                  }
               }
            } //if it is one QSS species
         } //loop on the reactants
      } //loop over all the reactions

      int species_with_max_interaction = 0;
      int max_interaction = 0;









      for (unsigned int k=0; k<listSpecies.size(); k++)
      {
         int counter = 0;

         if (Number_of_QSS_stoichiometric_coefficient_problem[k] != 0)
         {
            cout << listSpecies[k]->m_Name;
            cout << "Number of non-unity stoichiometric coefficients "
                 << Number_of_QSS_stoichiometric_coefficient_problem[k] << endl;
         }

         for (unsigned int kbis=0; kbis<listSpecies.size(); kbis++)
         {
            if (Number_of_QSS_interactions[k][kbis] != 0)
            {
               cout << listSpecies[k]->m_Name << " number of interactions with "
                    << listSpecies[kbis]->m_Name << "  " << Number_of_QSS_interactions[k][kbis] << endl;
               counter += Number_of_QSS_interactions[k][kbis];
            }
         }

         counter += Number_of_QSS_stoichiometric_coefficient_problem[k];

         if (counter != 0)
         {
            cout << "Total number of interactions " << counter << endl;
            total_counter += counter;
            cout << endl;

            if (counter > max_interaction)
            {
                max_interaction = counter;
                species_with_max_interaction = k;
            }
         }
      }

      if (total_counter <= 2)
         loop = false;
      else
      {
         QSS_Species[species_with_max_interaction] = false;
         cout << "-------------------------------------------------------------" << endl << endl;
      }

      getchar();

   } //end while(loop)


}
*/


string Write_QSS::print_sign(double input) const
{
   ostringstream output;

   output.precision(18);

   if (input > 0.0)
      output << "+" << fabs(input);
   else
      output << "-" << fabs(input);

   string string_output = output.str();

   return string_output;
}


int Write_QSS::Get_species_number(string species, vector<Species_ORCh*> listSpecies) const
{
   int species_number = -1;

   for (unsigned int k=0; k<listSpecies.size(); k++)
   {
      if (listSpecies[k]->m_Name == species)
         species_number = k;
   }

   return species_number;
}


Write_QSS::~Write_QSS() //Destructeur
{}


