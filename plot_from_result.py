import numpy as np

import os, shutil
import warnings
import sys, getopt
import matplotlib.pyplot as plt


Benchmark_Fold={"UBQ":"protg", "DEP":"myb_DNA-binding", "igvar-h":"Desulfoferrodox",
                "Rib_hydrolayse":"TIR", "Lum_binding":"EFTU_C", "DnaJ":"ATP-synt_DE_N",
                "Cohesin":"Fimbrial", "PAC":"GAF", "Agglutinin":"intb", "hemery":"SRP54", "LRR":"Recep_L_domain",
                "SSB":"TIMP", "TBCA": "BAG", "RRF":"RRF"}
Benchmark_SF={"His_biosynth":"OMPdecase", "svmp":"Astacin", "Lipoprotein_4":"mofe", "FAD-oxidase_NC": "MurB_C",
              "ETF_alpha":"Arginosuc_synth", "histone": "Arch_histone", "GA":"B", "PCNA":"DNA_PPF"} 

#Plot
def Affichage_Accuracy(List_Benchmark_Fold,  List_Benchmark_SF, Results, threshold, name):
    
    #stocke les accuracy
    acc=[0]*threshold
    
    
    for query in Results.keys():
        r=Results.get(query)
        rank_query=sorted(r, key=r.__getitem__) #Range les resultats du plus petit au plus grand
        rank_query.reverse() #list type
        if threshold>len(rank_query):
            threshold=len(rank_query)
            print("Pas assez de resultats possibles pour afficher le threshold demande dans le cas de "+query)
        accuracy_inter = 0  
        for i in range(threshold):     
            f=rank_query[i]
            #Fold (accorder un poids)
            if (f==List_Benchmark_Fold.get(query)) and (f!= None) and (List_Benchmark_Fold.get(query)!= None):
                accuracy_inter+=1
            #SF (accorder un poids)
            if (f==List_Benchmark_SF.get(query)) and (f!= None) and (List_Benchmark_SF.get(query)!= None):
                accuracy_inter+=1  
            acc[i]=acc[i]+accuracy_inter

    
    nb_queries=len(Results.keys())        
    for k in range(len(acc)):
        acc[k]=acc[k]/nb_queries


    #Area under curve
    Aire=np.trapz(acc,range(threshold), dx=1)
    Aire=round(Aire/threshold,2)
    #plot    
    fig=plt.figure()    
    #plt.plot(range(threshold), acc, 'ro')
    lw = 2
    plt.plot(range(threshold), acc, color='darkorange',
         lw=lw, label='Semi-ROC curve area '+str(Aire),linewidth=3)
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='--')
    plt.ylim((-0.01,1.01))
    plt.xlim((-1,405))
    plt.xlabel("Threshold")
    plt.ylabel("Number of positives")
    plt.title(name)
    plt.legend(loc="lower right")
    fig.savefig(name+".png")
    plt.show()
    
    return acc

def Accuracy_according_to_Fold_and_SF(List_Benchmark_Fold,  List_Benchmark_SF, Results, threshold):
        #stocke les accuracy
    accSF=[0]*threshold
    accF=[0]*threshold
    
    
    for query in Results.keys():
        r=Results.get(query)
        rank_query=sorted(r, key=r.__getitem__) #Range les resultats du plus petit au plus grand
        rank_query.reverse() #list type
        if threshold>len(rank_query):
            threshold=len(rank_query)
            print("Pas assez de resultats possibles pour afficher le threshold demande dans le cas de "+query)
        accuracy_SF = 0  
        accuracy_Fold = 0
        for i in range(threshold):     
            f=rank_query[i]
            #Fold (accorder un poids)
            if (f==List_Benchmark_Fold.get(query)) and (f!= None) and (List_Benchmark_Fold.get(query)!= None):
                accuracy_Fold+=1
            #SF (accorder un poids)
            if (f==List_Benchmark_SF.get(query)) and (f!= None) and (List_Benchmark_SF.get(query)!= None):
                accuracy_SF+=1  
            accSF[i]=accSF[i]+accuracy_SF
            accF[i]=accF[i]+accuracy_Fold

    
    nb_queries_Fold=len(List_Benchmark_Fold.keys()) 
    nb_queries_SF=len(List_Benchmark_SF.keys())   
    for k in range(len(accSF)):
        accSF[k]=accSF[k]/nb_queries_SF
    for k in range(len(accF)):       
        accF[k]=accF[k]/nb_queries_Fold
        

    #plot     
    fig=plt.figure()  
    plt.subplot(1, 2, 1)
   
    #Area under curve
    AireSF=np.trapz(accSF,range(threshold), dx=1)
    AireSF=round(AireSF/threshold,2)
         
    lw = 2
    plt.plot(range(threshold), accSF, color='darkorange',
         lw=lw, label='Semi-ROC curve area '+str(AireSF),linewidth=3)
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='--')
    plt.ylim((-0.01,1.01))
    plt.xlim((-1,405))
    plt.xlabel("Threshold")
    plt.ylabel("Number of positives")
    plt.title("Super-Family", fontsize=15)
    plt.legend(loc="lower right")

    plt.subplot(1, 2, 2)
    #Area under curve
    AireF=np.trapz(accF,range(threshold), dx=1)
    AireF=round(AireF/threshold,2)
         
    lw = 2
    plt.plot(range(threshold), accF, color='darkorange',
         lw=lw, label='Semi-ROC curve area '+str(AireF),linewidth=3)
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='--')
    plt.ylim((-0.01,1.01))
    plt.xlim((-1,405))
    plt.xlabel("Threshold")
    plt.ylabel("Number of positives")
    plt.title("Fold", fontsize=15)
    plt.legend(loc="lower right")
    
    
    plt.tight_layout()
    plt.show()   
    fig.savefig('Fold_and_SuperFamily.jpg')
    
    return

import seaborn as sns
def plot_score(List_Benchmark_Fold,  List_Benchmark_SF, Results):
    list_query=Results.keys()
    nombre_abscisse=len(list_query)

  
    x=range(nombre_abscisse)
    Top_one=[]
    for query in list_query:
        results_for_the_query=Results.get(query)
        for res in results_for_the_query.keys():
            if res==List_Benchmark_Fold.get(query) or List_Benchmark_SF.get(query)==res: 
                Top_one.append(results_for_the_query[res])
        
    fig=plt.figure()
    plt.plot(x, Top_one, '.', color='r', marker='s')
    
    compteur=0

    for query in list_query:
        results_for_the_query=Results.get(query)
        not_matched=[]
        for res in results_for_the_query.keys():
            not_matched.append(results_for_the_query[res])
        c = [0.0, float(np.std(not_matched))/float(80), float(80-np.std(not_matched))/float(80)] #R,G,B 
        print(np.std(not_matched))
        plt.scatter([compteur]*len(not_matched), not_matched, alpha=0.6, marker='.', color=c)
        compteur=compteur+1
    plt.xticks(x, list_query, rotation='vertical')
    plt.show()
    fig.savefig("details_on_each_query")
    
    
    
List_Full={"UBQ":["gln-synt_NC", "MAP1_LC3", 
                  "protg", "RA", "Staphylokinase", 
                  "UBX", "Band_41_M", "CIDE-N", 
                  "DUF170", "fer2"],
            "DEP":["hom", "HSF_DNA-bind", "HTH_AraC_2", "HTH_AraC", 
                   "Integrase_Zn", "IRF", "linker_histone", 
                   "Methyltransf_1", "Methyltransf_2", "myb_DNA-binding", 
                   "PAX", "recombinase", "Ribosomal_L11", "RNA_pol_N", "sigma54", 
                   "trans_reg_C", "Arg_repressor_N", "ARID", "ets", 
                   "Fe_dep_repress", "Fork_head", "GerE"],
            "igvar-h":["igC1", "igC2", "igcon", 
                       "igI", "igvar-l", "igV", "Invasin", 
                       "isoamylase_NC", "isoamylase_N", "malt_amylase_NC", 
                       "MHC_II_alpha_NC", "MHC_II_beta_NC", "MHC_II_C", 
                       "MSP_domain", "mycin Fold", "pili_assembly", "PKD", "RHD", 
                       "Rho_GDI", "sodcu", "TIG", "Transglutamin_C", "Transglutamin_NC", 
                       "Transglutamin_N", "Alpha_adaptinC12", "Alpha_adaptinC2", 
                       "alpha-amylase_N", "arrestin_C", "arrestin_NC", "arrestin", 
                       "cdh", "CopC", 
                       "Cyclodex_gly_tran", "Desulfoferrodox", "Filamin", "fn3", 
                       "Glyco_hydro_2"],          
    "His_biosynth": ["igps", "IMPDH", "MAAL_C", "MAAL", "Malate_synthase", 
                     "mle", "MM_CoA_mutase_N", "OMPdecase", "Orn_Arg_deC_N", 
                     "oxidored_FMN", "PEP-utilizers_C", "PEP-utilizers_NC", 
                     "phosphotriesterase", "piplc", "QRPTase", "RuBisCO_large_NC", 
                     "RuBisCO_large", "TGT", "tim", "Transaldolase", "trp_syntA", 
                     "URO-D", "xia", "ALAD", "Ala_racemase_N", "aldose", "aldosered", 
                     "alpha-amylase", "AP_endonuc_2", "bac_luciferase", "CBS", "DAHP_synth_1", 
                     "Dehydratase_LU", "DeoC", "DHDPS", "DHOdehase", "DHPS", "enolase", 
                     "F_bP_aldolase", "flavbb", "ghf10", "ghf14", "ghf17", "ghf18", "ghf1",
                     "ghf2", "ghf5", "Glu_syn_central", 
                    "Glyco_hydro_18_D1", "Glyco_hydro_18", 
                    "Glyco_hydro_20", "Glyco_hydro_67", "ICL"],
                     
      "RRF":["tRNA-synt_1d_N"],
      "Lum_binding":["EFG_IV", "EFTU_2", "FAD_binding", "reductases"],
      "Rib_hydrolayse":[ "AlaDh_PNT_D1", "DHquinase_II", "flav", 
                        "GATase", "gdh", "ligase-CoA_NC", "ligase-CoA", 
                        "MM_CoA_mutase", "NADHdh_2",
                        "ndk", "PAF-AH", "response_reg", "TIR"],
      "Fe_hyd_lg_C": ["CYCc", "E2_C", "EF1BD", "EFG_C", "EFG_IV", "fer4", 
                      "FTR", "GHMP_kinases_C", "HMA", "HPPK", 
                      "MCR_alpha_N", "MCR_beta_N", "MCR_gamma", 
                      "NTP_transf_2_2", "P-II", "PMSR", "Propep_M14", 
                      "PyrI_NC", "PyrI", "rrm", "RuBisCO_large_N", 
                      "Thr_dehydrat_C", "UreE", "ACT", "acyo"],
     "DnaJ": ["sodfe"],
     "Cohesin":["Fimbrial", "T-box", "A2M_A", "Apocytochrome_F", "CBD_2", "CBD_3"],
     "ETF_alpha":["tRNA-synt_1b", "tRNA-synt_1c", "tRNA-synt_1d_M", "Usp", 
                  "Arginosuc_synth", "Asn_synthase_NC", "Asn_synthase",
                  "Cytidylyltransf", "ETF_beta"],
                  
     "PAC":["GAF", "PAS", "profilin"],
     "Agglutinin": ["intb", "Ricin_B_lectin", "sti"],
      "hemery": ["Lipoprotein_6", "mcp2", "SRP54", 
                 "TMV_coat", "Yers_vir_YopE", "cytprime"],
      "LRR": ["Recep_L_domain", "Furin-like", "internalin"],
      "svmp": ["mmp", "Peptidase_M27", "tln", 
               "Astacin", "Glyco_hydro_20b"],
      "Lipoprotein_4": ["Peripla_BP_2", "Ferrochelatase", "mofe"],
      "SSB":["Pertussis_S2S3_NC", "Pertussis_S2S3", "ppase", 
             "Ribosomal_L2_N", "RuvA", "S1", "SLT_beta", "Stap_Strp_toxin", 
             "TIMP", "tRNA_bind", "asprs", "CcmE", "csp", "DNA_ligase_C", 
             "eIF-1a", "eIF-5a", "ltb", "ModE_C"],
     "TBCA": ["SPEC", "BAG", "FAD_binding_2"],
     "FAD-oxidase_NC": ["FAD-oxidase_C"],  
     "GA":["B"],
     "histone": ["Arch_histone"],
     "PCNA":["DNA_PPF"]
    }
    
    
def full_plot(Results, List_Full, name="Enrichment with more than one possible positive per query", threshold=405):
    #stocke les accuracy
    acc=[0]*threshold
   
    
    Number_of_positiv=0
    for query in Results.keys():
        r=Results.get(query)
        rank_query=sorted(r, key=r.__getitem__) #Range les resultats du plus petit au plus grand
        rank_query.reverse() #list type
        if threshold>len(rank_query):
            threshold=len(rank_query)
            print("Pas assez de resultats possibles pour afficher le threshold demande dans le cas de "+query)
        accuracy_inter = 0  
        for i in range(threshold):     
            f=rank_query[i]
            list_to_check=List_Full[query]
            if f in list_to_check:
                accuracy_inter=accuracy_inter+1
                Number_of_positiv+=1
            acc[i]=acc[i]+accuracy_inter
     
    for k in range(len(acc)):
        acc[k]=acc[k]/Number_of_positiv
        
    #Area under curve
    Aire=np.trapz(acc,range(threshold), dx=1)
    Aire=round(Aire/threshold,2)
    #plot    
    fig=plt.figure()    
    #plt.plot(range(threshold), acc, 'ro')
    lw = 2
    plt.plot(range(threshold), acc, color='darkorange',
         lw=lw, label='Semi-ROC curve area '+str(Aire),linewidth=3)
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='--')
    plt.xlabel("Threshold")
    plt.ylabel("True positiv ratio")
    plt.ylim((-0.01,1.01))
    plt.xlim((-1,405))
    plt.title(name)
    plt.legend(loc="lower right")
    fig.savefig(name+".png")
    plt.show()
    
    return 


r1=np.load("eucli4.npy").item()
#Affichage_Accuracy(Benchmark_Fold, Benchmark_SF, r1, 405, "Global")
#Accuracy_according_to_Fold_and_SF(Benchmark_Fold, Benchmark_SF, r1, 405)
#plot_score(Benchmark_Fold, Benchmark_SF, r1)
full_plot(r1, List_Full)