def addGPR2Models(model,cyc):
    '''
    This function uses the function extractGeneAndProteinAssociation to add GPR
    associations to reactions with missing GPR in the model
    Input: 1) cobra model 2) pythonCyc PGDB instance
    Output: cobra model
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    reactions = cyc.reactions
    rxnPresentList = list()
    rxnIDed = dict()
    for CycRxn in reactions.instances:
        CycRxn_id = CycRxn.frameid
        CycRxn_id_adapted = convertCycID2sbmlID(CycRxn_id)
        tempList = list()
        for rxn in model.reactions:
            if CycRxn_id_adapted == rxn.id[0:rxn.id.rindex("_")]:
                tempList.append(rxn)
            elif CycRxn_id_adapted == rxn.id[0:rxn.id.rindex("_")].replace("_NADP","").replace("_NAD",""):
                tempList.append(rxn)
        rxnIDed[CycRxn_id]=tempList
    SoyIgnoreList = ["RXN_9650_p","2_KETO_ADIPATE_DEHYDROG_RXN_m","Phytol_biosynthesis_p" \
                ,"CYSTEINE_AMINOTRANSFERASE_RXN_m","GLYCINE_TRNA_LIGASE_RXN_c" \
                ,"RXN66_1_c","RXN_9648_p","RXN-9651","Plastidial_ATP_Synthase_p" \
                ,"GGPP_biosynthesis_p","RXN_9653_p","lycopene_biosynthesis_p" \
                ,"RXN_2141_p","SUCCINYL_COA_HYDROLASE_RXN_m","PROTON_ATPase_c" \
                ,"MDA_Fd_Ascorbate_p","MercaptoPyruvateSulfurtransferase_m" \
                ,"Phytol_degradation_p","RXN_9652_p","A_B_oxidation_x","unlProtHYPO_c" \
                ,"Mitochondrial_ATP_Synthase_m","IPP_biosynthesis_c","Mehler_Reaction_p" \
                ,"Beta_Oxidation_x","HMBPP_synthesis_p","OROTATE_REDUCTASE_NADH_RXN_p" \
                ,"Ferredoxin_Plastoquinone_Reductase_p","RXN_9651_p","NADPH_Dehydrogenase_p" \
                ,"Plastoquinol_Oxidase_p","SUCCINATE_COA_LIGASE_GDP_FORMING_RXN_m","RXN_1781_v" \
                ,"PREPHENATE_DEHYDROGENASE_NADP_RXN_p","PREPHENATEDEHYDROG_RXN_p" \
                ,"MALEYLACETOACETATE_ISOMERASE_RXN_c","RXN_9654_p","LCYSDESULF_RXN_c","RXN_9958_NAD_m" \
                ,"HEXOKINASE_RXN_MANNOSE_c","PYRUVDEH_RXN_p","PYRUVDEH_RXN_m"] #last 3 lines present in latest version of SoyCyc
    print("--------------\nThis list of metabolic reactions are ignored")
    print(SoyIgnoreList)
    print("--------------")
    IDedlist = set()
    for rxnlist in rxnIDed.values():
        IDedlist = IDedlist.union(set(rxnlist))
    for rxn in set(model.reactions) - IDedlist:
        if not("_tx" in rxn.id or "_pc" in rxn.id or \
               "_mc" in rxn.id or "_xc" in rxn.id or \
               "_im" in rxn.id or "_vc" in rxn.id or \
               "_ec" in rxn.id or "_ep" in rxn.id or \
               "_pr" in rxn.id) \
               and (not "Biomass" in rxn.id) and \
               (not "biomass" in rxn.id) and \
               (not "Protein" in rxn.id) and \
               (not "TRNA_LIGASE" in rxn.id):
               if rxn.id not in SoyIgnoreList:
                   print(rxn.id)
    for k in rxnIDed.keys():
        for v in rxnIDed.get(k):
            rxn = v
            if rxn.gene_reaction_rule == "":
                #print k
                GPR = extractGeneAndProteinAssociation(cyc,k)
                if GPR != "()":
                    GPR = GPR.replace("() or ","")
                    rxn.gene_reaction_rule = GPR
    return model

def extractGeneAndProteinAssociation(cyc,frame_id):
    '''
    This functions adds Gene Associations to cobra model from Pathway Tools via
    PythonCyc
    Input: 1) pythonCyc PGDB instance 2) Frame id of reaction from Pathway Tools
    Output: Gene-Protein-Reaction associations from a PGDB
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    rxn = getFrame(cyc,frame_id)
    if frame_id in cyc.reactions.instances:
        print("Error check if "+frame_id+" is reaction")
        return ""
    else:
        if "enzymatic_reaction" not in dir(rxn):
            return ""
        else:
            enzrxns = cyc.get_frame_objects(rxn.enzymatic_reaction)
            GPR = "(GPR)"
            temp1 = ""
            for enzrxn in enzrxns:
                enz = getFrame(cyc,enzrxn.enzyme)
                if "names" not in dir(enz):
                    continue
                if "gene" not in dir(enz):
                    continue
                if temp1 == "":
                    temp1 = str(enz.frameid)
                else:
                    temp1 = temp1 +" or "+str(enz.frameid)
                temp2 = ""
                for gene in enz.gene:
                    gene = getFrame(cyc,gene)
                    if "accession_1" not in dir(gene):
                        temp1.replace(enz.frameid,"")
                        continue
                    if temp2 == "":
                        temp2 = gene.accession_1
                    else:
                        temp2 = temp2 +" or "+gene.accession_1
                #print temp1
                #print temp2
                temp1 = temp1.replace(enz.frameid,"("+temp2+")")
            GPR = GPR.replace("GPR",temp1)
            return GPR

def getFrame(cyc,frame_id):
    '''
    This function retrieves pythoncyc frame from a PGDB instance
    Input: 1) pythonCyc PGDB instance 2) Frame id from Pathway Tools
    Output: Python instance of a frame
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    frame = cyc.get_frame_objects([frame_id])[0]
    return frame

def convertCycID2sbmlID(id):
    '''
    This function converts Pathway Tools IDs to one that is SBML compliant
    Input: BioCyc IDs
    Output: SBML compliant IDs
    Author: Sanu Shameer (sanushameer@gmail.com)
    '''
    new_id = id.replace(".","_PERIOD_")
    new_id = new_id.replace("%2b","_")
    new_id = new_id.replace("&#039;","_")
    new_id = new_id.replace("&amp;beta;","B")
    new_id = new_id.replace("&beta;","B")
    new_id = new_id.replace("|","")
    new_id = new_id.replace("+-","_")
    new_id = new_id.replace("--","_")
    new_id = new_id.replace("-","_")
    new_id = new_id.replace("+","_")
    new_id = new_id.replace("'","_")
    new_id = new_id.replace("(","_")
    new_id = new_id.replace(")","_")
    new_id = new_id.replace("/","_")
    new_id = new_id.replace("__","_")
    return new_id

def find_average(temp_list):
    '''
    This function calculates the average from a list of numbers
    Input: list
    Output: float
    Author:Sanu Shameer (sanushameer@gmail.com)
    '''
    return sum(temp_list)/len(temp_list)

def adjustObjectiveBasedOnDaylength(diel_leaf,daylength,obj_rxn="diel_biomass",Phloem_day = "X_Phloem_contribution_t1",Phloem_night="X_Phloem_contribution_t2"):
    '''
    This function adjusts the day/night phloem contribution
    in a diel leaf objective.
    Input:  1) diel model, 2) length of day in hours, 3)ID of objective,
            4) Metabolite ID of day-time phloem, 5) Metabolite ID of night-time
            phloem
    Output: 1)diel model
    '''
    ratio = float(3*daylength)/(24-daylength)
    rxn = diel_leaf.reactions.get_by_id(obj_rxn)
    met1 = diel_leaf.metabolites.get_by_id(Phloem_day)
    coeff = abs(rxn.metabolites.get(met1))
    new_coeff = ratio - coeff
    rxn.add_metabolites({met1:-1*new_coeff})
    return diel_leaf

def generateMetaboliteFormula(rxn):
    count = 0
    for met in rxn.metabolites:
        if met.formula=="" or met.formula=="NA" or met.formula == None:
            if met.formula == "NA" or met.formula == None:
                met.formula = ""
            count = count + 1
    if count == 1:
        unb = rxn.check_mass_balance()
        #print(unb.keys())
        for met in rxn.metabolites:
            stoich = rxn.metabolites[met]
            if met.formula == "":
                tempForm = ""
                for a in ["C","H","O"]:
                    if a in unb.keys():
                        if round(unb[a]/stoich,6)==0:
                            continue
                        tempForm = tempForm+a+str(abs(unb[a])/stoich)
                        #print(a)
                        #print(unb[a])
                        #print(stoich)
                        #print(abs(unb[a])/stoich)
                for a in unb.keys():
                    if a in ["C","H","O"]:
                        continue
                    if a=="charge" or round(unb[a]/stoich,6)==0:
                        continue
                    tempForm = tempForm+a+str(abs(unb[a])/stoich)
                met.formula = tempForm
                print(met.id)
                print(tempForm)
    else:
        print("Unable to generate missing metabolite formula")

def removeSpecificMetChargedState(model,metlist):
    for met in metlist:
        met = model.metabolites.get_by_id(met)
        rxn2edit = set(met.reactions)
        defaultForm = model.metabolites.get_by_id(met.id[1:])
        met.remove_from_model()
        for rxn in rxn2edit:
            rxn.add_metabolites({defaultForm:-0.03,
                                 model.metabolites.get_by_id("PROTON_"+defaultForm.compartment):-0.03})
        return model

def updateFAcomposition(model,organ,biomass):

    temp = model.copy()
    temp.reactions.FattyAcid_composition_p.remove_from_model()
    temp.metabolites.Fatty_Acids_p.formula=""
    temp.metabolites.Fatty_Acids_c.formula=""
    temp.metabolites.Long_Chain_Acyl_CoAs_p.formula=""
    from cobra.core import Reaction
    FACP = {"PALMITATE_p":"Palmitoyl_ACPs_p",
            "CPD_9245_p":"Palmitoleoyl_ACP_p",
            "CPD_17412_p":"hexadecadienoate_ACP_p",
            "CPD_17291_p":"hexadecatrienoate_ACP_p",
            "STEARIC_ACID_p":"Stearoyl_ACPs_p",
            "OLEATE_CPD_p":"Oleoyl_ACPs_p",
            "Octadecadienoate_p":"Octadecadienoyl_ACP_p",
            "LINOLENIC_ACID_p":"Octadecatrienoyl_ACP_p",
            "ARACHIDIC_ACID_p":"Arachidoyl_ACPs_p",
            "CPD_16709_p":"Eicosenoyl_ACP_p",
            "DOCOSANOATE_p":"Behenoyl_ACPs_p"}

    PLs = ["ACYL_SN_GLYCEROL_3P_p","L_PHOSPHATIDATE_p","L_PHOSPHATIDATE_m","DIACYLGLYCEROL_p",
           "DIACYLGLYCEROL_r","Triacylglycerols_p","PHOSPHATIDYL_CHOLINE_r",
           "L_1_PHOSPHATIDYL_ETHANOLAMINE_r","L_1_PHOSPHATIDYL_GLYCEROL_p",
           "L_1_PHOSPHATIDYL_GLYCEROL_P_p","L_1_PHOSPHATIDYL_GLYCEROL_P_m",
           "L_1_PHOSPHATIDYL_GLYCEROL_m","2_Lysophosphatidylcholines_r",
           "Lysophosphatidylglycerols_r","CDPDIACYLGLYCEROL_p","CDPDIACYLGLYCEROL_m",
           "D_Galactosyl_12_diacyl_glycerols_p","Galactosyl_galactosyl_diacyl_glycerols_p"]


    for met in PLs:
        met=temp.metabolites.get_by_id(met)
        met.formula=""

    FAdict = dict(biomass[biomass["type"]=="fattyacid"][organ])

    k = organ
    RXN1 = Reaction("Fatty_acid_mix_"+k)
    RXN2 = Reaction("Fatty_acid_ACP_"+k)
    tot = 0
    for met in FAdict.keys():
        RXN1.add_metabolites({temp.metabolites.get_by_id(met):-1*FAdict[met]})
        RXN2.add_metabolites({temp.metabolites.get_by_id(FACP[met]):-1*FAdict[met]})
        tot = tot+FAdict[met]
    print(tot)
    if tot==0:
        RXN1.add_metabolites({temp.metabolites.PALMITATE_p:-1})
        RXN2.add_metabolites({temp.metabolites.Palmitoyl_ACPs_p:-1})
        tot = 1
    RXN1.add_metabolites({temp.metabolites.Fatty_Acids_p:tot})
    RXN1.lower_bound = 1000
    RXN1.upper_bound = 0
    temp.add_reaction(RXN1)

    RXN2.add_metabolites({temp.metabolites.Fatty_acyl_ACP_p:tot})
    RXN2.lower_bound = 1000
    RXN2.upper_bound = 0
    temp.add_reaction(RXN2)

    generateMissingFormula(temp)

    return temp

def generateMissingFormula(model,debug=False):
    loop_counter = 0
    former = 0
    for met in model.metabolites:
        if met.formula == "" or met.formula == "NA":
            former = former +1
    latter = 1
    while True:
        loop_counter = loop_counter+1
        if debug:
            print("Loop = "+str(loop_counter))
        former = latter
        for rxn in model.reactions:
            count = 0
            for met in rxn.metabolites:
                if met.formula=="" or met.formula=="NA" or met.formula == None:
                    if met.formula == "NA" or met.formula == None:
                        met.formula = ""
                    count = count + 1
                    coeff = rxn.metabolites[met]
            if count == 1:
                unb = rxn.check_mass_balance()
                eqn = rxn.reaction
                eqn = " "+eqn+" "
                for met in rxn.metabolites.keys():
                    formula = met.formula
                    if formula == None:
                        formula = "0"
                        NF_list.add(rxn.id)
                    eqn=eqn.replace(" "+met.id+" ","("+formula+")")
                if debug:
                    print(eqn)
                    print(unb)
                for met in rxn.metabolites:
                    if met.formula == "":
                        tempForm = ""
                        for a in sorted(unb.keys()):
                            if a=="charge" or round(unb[a],2)==0:
                                continue
                            num = float(abs(unb[a]))/abs(coeff)
                            if str(round(num))==str(num):
                                tempForm = tempForm+a+str(int(round(num)))
                            else:
                                tempForm = tempForm+a+str(num)
                                #print(a)
                                #print(round(num)==num)
                                #print(round(num))
                                #print(num)
                                #print(type(round(num)))
                                #print(type(num))
                        met.formula = tempForm
                        if debug:
                            print(met.id)
                            print(tempForm)
        latter = 0
        for met in model.metabolites:
            if met.formula == "" or met.formula == "NA":
                latter = latter +1
        if former == latter:
            break

def generateStemModel(model):
	from cobra.core import Reaction
	for met in model.reactions.Phloem_output_tx.metabolites.keys():
		met2 = met.copy()
		if met.id=="sSUCROSE_b":
			met2.id = "SUCROSE_ph"
			met = model.metabolites.get_by_id("SUCROSE_c")
		elif "PROTON" in met.id:
			continue
		else:
			met2.id = met.id.replace("_c","_ph")
		met2.compartment = "ph"
		model.add_metabolites(met2)
		rxn = Reaction(met2.id+"_exchange")
		rxn.add_metabolites({met2:1})
		model.add_reaction(rxn)
		rxn = Reaction(met2.id.replace("_ph","_phloem_uptake"),name=met2.id.replace("_ph","_phloem_uptake"))
		rxn.add_metabolites({met2:-1,model.metabolites.get_by_id("PROTON_e"):-1,
                             met:1,model.metabolites.get_by_id("PROTON_c"):1})
		rxn.lower_bound = 0
		rxn.upper_bound = 1000
		model.add_reaction(rxn)
        #print(rxn.reaction)
	return model

def generateRootModel(model,symbiont=None):
	from cobra.core import Reaction
	for met in model.reactions.Phloem_output_tx.metabolites.keys():
		met2 = met.copy()
		if met.id=="sSUCROSE_b":
			met2.id = "SUCROSE_ph"
			met = model.metabolites.get_by_id("SUCROSE_c")
		elif "PROTON" in met.id:
			continue
		else:
			met2.id = met.id.replace("_c","_ph")
		met2.compartment = "ph"
		model.add_metabolites(met2)
		rxn = Reaction(met2.id+"_exchange")
		rxn.add_metabolites({met2:1})
		model.add_reaction(rxn)

		rxn = Reaction(met2.id.replace("_ph","_phloem_uptake"),name=met2.id.replace("_ph","_phloem_uptake"))
		rxn.add_metabolites({met2:-1,model.metabolites.get_by_id("PROTON_e"):-1,
                             met:1,model.metabolites.get_by_id("PROTON_c"):1})
		rxn.lower_bound = 0
		rxn.upper_bound = 1000
		model.add_reaction(rxn)
	#add xylem reactions
	for met in ["CAII","MGII","KI","NITRATE","SULFATE","AMMONIUM","WATER","GLT","L_ASPARTATE","ASN","GLN"]:
		met2 = model.metabolites.get_by_id(met+"_c").copy()
		met2.id = met+"_xy"
		met2.compartment = "xy"

		rxn = Reaction(met+"_exchange")
		rxn.add_metabolites({met2:-1})
		rxn.lower_bound = 0
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		rxn = Reaction(met+"_xylem_export")
		rxn.add_metabolites({model.metabolites.get_by_id(met+"_c"):-1,met2:1})
		model.add_reaction(rxn)

	if symbiont == None:
		return model
	else:
		model.reactions.Nitrate_tx.upper_bound = 0
		model.reactions.Nitrate_tx.lower_bound = 0

		#adding symbiont compartment
		from cobra import io
		rhizo = io.read_sbml_model(symbiont["path"])
		for met in rhizo.metabolites:
			met.compartment = met.compartment+"_rhizo"
		rhizo.compartments={"c_rhizo":"rhizobe cytosol","e_rhizo":"rhizobe extracellular"}
		model = model+rhizo

		#     rxn = Reaction("Sucrose_exchange_symbiont")
		#     rxn.name = rxn.id.replace("_"," ")
		#     rxn.add_metabolites({model.metabolites.get_by_id("SUCROSE_c"):-1,model.metabolites.get_by_id("cpd00076[e0]"):1})
		#     rxn.lower_bound = -1000
		#     rxn.upper_bound = 1000
		#     model.add_reaction(rxn)

		rxn = Reaction("Alanine_exchange_symbiont")
		rxn.name = rxn.id.replace("_"," ")
		rxn.add_metabolites({model.metabolites.get_by_id("L_ALPHA_ALANINE_c"):-1,
                             model.metabolites.get_by_id("ala__L[e]"):1})
		rxn.lower_bound = -1000
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		rxn = Reaction("Aspartate_exchange_symbiont")
		rxn.name = rxn.id.replace("_"," ")
		rxn.add_metabolites({model.metabolites.get_by_id("L_ASPARTATE_c"):-1,
                             model.metabolites.get_by_id("asp__L[e]"):1})
		rxn.lower_bound = -1000
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		rxn = Reaction("Glutamate_exchange_symbiont")
		rxn.name = rxn.id.replace("_"," ")
		rxn.add_metabolites({model.metabolites.get_by_id("GLT_c"):-1,
                             model.metabolites.get_by_id("glu__L[e]"):1})
		rxn.lower_bound = -1000
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		rxn = Reaction("Malate_exchange_symbiont")
		rxn.name = rxn.id.replace("_"," ")
		rxn.add_metabolites({model.metabolites.get_by_id("MAL_c"):-1,
                             model.metabolites.get_by_id("mal__L[e]"):1})
		rxn.lower_bound = -1000
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		rxn = Reaction("Succinate_exchange_symbiont")
		rxn.name = rxn.id.replace("_"," ")
		rxn.add_metabolites({model.metabolites.get_by_id("SUC_c"):-1,
                             model.metabolites.get_by_id("succ[e]"):1})
		rxn.lower_bound = -1000
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		rxn = Reaction("Ammonium_exchange_symbiont")
		rxn.name = rxn.id.replace("_"," ")
		rxn.add_metabolites({model.metabolites.get_by_id("AMMONIUM_c"):-1,
                             model.metabolites.get_by_id("fixedNH3[e]"):1})
		rxn.lower_bound = -1000
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		return model

def generateSeedModel(model):
	from cobra.core import Reaction
	for met in model.reactions.Phloem_output_tx.metabolites.keys():
		met2 = met.copy()
		if met.id=="sSUCROSE_b":
			met2.id = "SUCROSE_ph"
			met = model.metabolites.get_by_id("SUCROSE_c")
		elif "PROTON" in met.id:
			continue
		else:
			met2.id = met.id.replace("_c","_ph")
		met2.compartment = "ph"
		model.add_metabolites(met2)

		rxn = Reaction(met2.id+"_exchange")
		rxn.add_metabolites({met2:1})
		model.add_reaction(rxn)

		rxn = Reaction(met2.id.replace("_ph","_phloem_uptake"),name=met2.id.replace("_ph","_phloem_uptake"))
		rxn.add_metabolites({met2:-1,model.metabolites.get_by_id("PROTON_e"):-1,
                             met:1,model.metabolites.get_by_id("PROTON_c"):1})
		rxn.lower_bound = 0
		rxn.upper_bound = 1000
		model.add_reaction(rxn)

		#print(rxn.reaction)
	return model

def createEmptyBiomassDataFrame():
	import pandas as pd
	biomass = pd.DataFrame(data={"":["sSUCROSE_b","GLC_c","FRU_c","Starch_b","Cellulose_b","Xylan_b",
	                                 "L_PHOSPHATIDATE_p","PHOSPHATIDYL_CHOLINE_r",
	                                 "L_1_PHOSPHATIDYL_ETHANOLAMINE_r","DIACYLGLYCEROL_p",
	                                 "Galactosyl_galactosyl_diacyl_glycerols_p",
	                                 "D_Galactosyl_12_diacyl_glycerols_p","2_Lysophosphatidylcholines_r",
	                                 "Lysophosphatidylglycerols_r","Triacylglycerols_p",
	                                 "L_1_PHOSPHATIDYL_GLYCEROL_p","L_1_phosphatidyl_inositols_r",
	                                 "SULFOQUINOVOSYLDIACYLGLYCEROL_p","Protein_b",
	                                 "sMAL_b","sCIT_b","sFUM_b","ARG_c","HIS_c","LYS_c","sASP_b",
	                                 "sGLU_b","sSER_b","THR_c","ASN_c","sGLN_b","CYS_c",
	                                 "GLY_c","PRO_c","sALA_b","VAL_c","ILE_c","LEU_c",
	                                 "MET_c","PHE_c","TYR_c","TRP_c","sGABA_b","PALMITATE_p",
	                                 "CPD_9245_p","CPD_17412_p","CPD_17291_p","STEARIC_ACID_p",
	                                 "OLEATE_CPD_p","Octadecadienoate_p","LINOLENIC_ACID_p",
	                                 "ARACHIDIC_ACID_p","CPD_16709_p","DOCOSANOATE_p",
	                                 "SUC_c","FUM_c","MAL_c","CIS_ACONITATE_c","CIT_c","MYO_INOSITOL_c",
	                                 "pHIS_b","pILE_b","pTHR_b","pARG_b","pASN_b","pGLU_b","pPHE_b",
	                                 "pGLN_b","pTYR_b","pMET_b","pASP_b","pVAL_b","pLYS_b","pSER_b",
	                                 "pGLY_b","pALA_b","pLEU_b","pPRO_b","pCYS_b","pTRP_b","COUMARATE_c"],
	                             "type":[""]*81,
	                             "leaf":[0.0]*81,"stem":[0.0]*81,"root":[0.0]*81,"seed":[0.0]*81,},dtype="float64")
	biomass = biomass.set_index("")
	for i in ["pHIS_b","pILE_b","pTHR_b","pARG_b","pASN_b","pGLU_b","pPHE_b","pGLN_b","pTYR_b","pMET_b",
	          "pASP_b","pVAL_b","pLYS_b","pSER_b","pGLY_b","pALA_b","pLEU_b","pPRO_b","pCYS_b","pTRP_b"]:
	    biomass.at[i,"type"]="protein"
	for i in ["PALMITATE_p","CPD_9245_p","CPD_17412_p","CPD_17291_p","STEARIC_ACID_p",
	          "OLEATE_CPD_p","Octadecadienoate_p","LINOLENIC_ACID_p",
	          "ARACHIDIC_ACID_p","CPD_16709_p","DOCOSANOATE_p"]:
	    biomass.at[i,"type"]="fattyacid"
	return biomass
