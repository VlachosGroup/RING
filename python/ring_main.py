import pyring
import pyring.react

#def checkCCRt_COFormationfromOHm0():
#    print("inside the function")
##    substruct = pyring.react.Substructure("O1(-C2)", pyring.react.patternsize("O1(-C2)"))
##    c_substruct = (substruct)._ptr
##    #c_substruct = (<pyring.react.Substructure?>substruct)._ptr
##    pattern = pyring.react.Patternmatch(mol, substruct[0], 1)
##    value = pattern.GetDistinctMatches()
##    print(value)
##    if value == 0:
##        return True
##    else:
##        return False
#
#print("after definition")
##print(pyring.react.GlobalConstraints)
#pyring.react.GlobalConstraints.append(checkCCRt_COFormationfromOHm0())
#print("after definition")
#bo = checkCCRt_COFormationfromOHm0(pyring.react.Molecule("O1(-C2)", 2))
#print(bo)
reactantlist = list()
CompositeAtomsList=[]
CompositeSiteList=[]
Rtypelist = []
reactantlist.append("CCC")
reactantlist.append("O=O")
reactantlist.append("{Pt}")
CompositeAtomsList.append("{Pt}")

Rt_O2Adsorb = pyring.react.ReactionType()
Rt_O2Adsorb.add_reactant_pattern(pyring.react.Substructure("O1(=O2)", pyring.react.patternsize("O1(=O2)")))
Rt_O2Adsorb.add_reactant_pattern(pyring.react.Substructure("{Pt}3[!|~$|>0]", pyring.react.patternsize("{Pt}3[!|~$|>0]")))
Rt_O2Adsorbdupmap2= {}
Rt_O2Adsorbdupmap2[3] = 4
Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap2)
Rt_O2Adsorbdupmap3= {}
Rt_O2Adsorbdupmap3[3] = 5
Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap3)
Rt_O2Adsorbdupmap4= {}
Rt_O2Adsorbdupmap4[3] = 6
Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap4)
Rt_O2Adsorb.disconnect_bond(1,2)
Rt_O2Adsorb.connect_bond(1,3)
Rt_O2Adsorb.connect_bond(1,4)
Rt_O2Adsorb.connect_bond(2,5)
Rt_O2Adsorb.connect_bond(2,6)
Rt_O2Adsorb.set_rule_name("O2Adsorb")
Rtypelist.append(Rt_O2Adsorb) 

Rt_CCScission = pyring.react.ReactionType()
Rt_CCScission.add_reactant_pattern(pyring.react.Substructure("C1(-C2)", pyring.react.patternsize("C1(-C2)")))
Rt_CCScission.add_reactant_pattern(pyring.react.Substructure("{Pt}3[!|~$|>0]", pyring.react.patternsize("{Pt}3[!|~$|>0]")))
Rt_CCScissiondupmap2= {}
Rt_CCScissiondupmap2[3] = 4
Rt_CCScission.AddReactantCopy(1, Rt_CCScissiondupmap2)
Rt_CCScission.disconnect_bond(1,2)
Rt_CCScission.connect_bond(1,3)
Rt_CCScission.connect_bond(2,4)
Rt_CCScission.set_rule_name("CCScission")
Rtypelist.append(Rt_CCScission) 

Rt_CHScission = pyring.react.ReactionType()
Rt_CHScission.add_reactant_pattern(pyring.react.Substructure("C1(-H2)", pyring.react.patternsize("C1(-H2)")))
Rt_CHScission.add_reactant_pattern(pyring.react.Substructure("{Pt}3[!|~$|>0]", pyring.react.patternsize("{Pt}3[!|~$|>0]")))
Rt_CHScissiondupmap2= {}
Rt_CHScissiondupmap2[3] = 4
Rt_CHScission.AddReactantCopy(1, Rt_CHScissiondupmap2)
Rt_CHScission.disconnect_bond(1,2)
Rt_CHScission.connect_bond(1,3)
Rt_CHScission.connect_bond(2,4)
Rt_CHScission.set_rule_name("CHScission")
Rtypelist.append(Rt_CHScission) 

Rt_OHScission = pyring.react.ReactionType()
Rt_OHScission.add_reactant_pattern(pyring.react.Substructure("O1(-H2)", pyring.react.patternsize("O1(-H2)")))
Rt_OHScission.add_reactant_pattern(pyring.react.Substructure("{Pt}3[!|~$|>0]", pyring.react.patternsize("{Pt}3[!|~$|>0]")))
Rt_OHScissiondupmap2= {}
Rt_OHScissiondupmap2[3] = 4
Rt_OHScission.AddReactantCopy(1, Rt_OHScissiondupmap2)
Rt_OHScission.disconnect_bond(1,2)
Rt_OHScission.connect_bond(1,3)
Rt_OHScission.connect_bond(2,4)
Rt_OHScission.set_rule_name("OHScission")
Rtypelist.append(Rt_OHScission) 

Rt_HTransferCtoO = pyring.react.ReactionType()
Rt_HTransferCtoO.add_reactant_pattern(pyring.react.Substructure("C1(-H2)", pyring.react.patternsize("C1(-H2)")))
Rt_HTransferCtoO.add_reactant_pattern(pyring.react.Substructure("{Pt}4(-O3(-{Pt}5))", pyring.react.patternsize("{Pt}4(-O3(-{Pt}5))")))
Rt_HTransferCtoO.disconnect_bond(1,2)
Rt_HTransferCtoO.disconnect_bond(3,4)
Rt_HTransferCtoO.connect_bond(1,4)
Rt_HTransferCtoO.connect_bond(3,2)
Rt_HTransferCtoO.set_rule_name("HTransferCtoO")
Rtypelist.append(Rt_HTransferCtoO) 

Rt_HTransferOHtoO = pyring.react.ReactionType()
Rt_HTransferOHtoO.add_reactant_pattern(pyring.react.Substructure("O1(-H2)", pyring.react.patternsize("O1(-H2)")))
Rt_HTransferOHtoO.add_reactant_pattern(pyring.react.Substructure("{Pt}4(-O3(-{Pt}5))", pyring.react.patternsize("{Pt}4(-O3(-{Pt}5))")))
Rt_HTransferOHtoO.disconnect_bond(1,2)
Rt_HTransferOHtoO.disconnect_bond(3,4)
Rt_HTransferOHtoO.connect_bond(1,4)
Rt_HTransferOHtoO.connect_bond(3,2)
Rt_HTransferOHtoO.set_rule_name("HTransferOHtoO")
Rtypelist.append(Rt_HTransferOHtoO) 

Rt_COFormationfromOH = pyring.react.ReactionType()
Rt_COFormationfromOH.add_reactant_pattern(pyring.react.Substructure("{Pt}2(-C1[!|~C|>1])", pyring.react.patternsize("{Pt}2(-C1[!|~C|>1])")))
Rt_COFormationfromOH.add_reactant_pattern(pyring.react.Substructure("{Pt}5(-O3(-H4))", pyring.react.patternsize("{Pt}5(-O3(-H4))")))
Rt_COFormationfromOH.disconnect_bond(3,5)
Rt_COFormationfromOH.disconnect_bond(1,2)
Rt_COFormationfromOH.connect_bond(1,3)
Rt_COFormationfromOH.set_rule_name("COFormationfromOH")
Rtypelist.append(Rt_COFormationfromOH) 

Rt_COFormationfromO = pyring.react.ReactionType()
Rt_COFormationfromO.add_reactant_pattern(pyring.react.Substructure("{Pt}2(-C1[!|~C|>1])", pyring.react.patternsize("{Pt}2(-C1[!|~C|>1])")))
Rt_COFormationfromO.add_reactant_pattern(pyring.react.Substructure("{Pt}4(-O3(-{Pt}5))", pyring.react.patternsize("{Pt}4(-O3(-{Pt}5))")))
Rt_COFormationfromO.disconnect_bond(3,4)
Rt_COFormationfromO.disconnect_bond(1,2)
Rt_COFormationfromO.connect_bond(1,3)
Rt_COFormationfromO.set_rule_name("COFormationfromO")
Rtypelist.append(Rt_COFormationfromO) 

Rt_H2OwCFormation = pyring.react.ReactionType()
Rt_H2OwCFormation.add_reactant_pattern(pyring.react.Substructure("C1(-H2)", pyring.react.patternsize("C1(-H2)")))
Rt_H2OwCFormation.add_reactant_pattern(pyring.react.Substructure("{Pt}5(-O3(-H4))", pyring.react.patternsize("{Pt}5(-O3(-H4))")))
Rt_H2OwCFormation.disconnect_bond(5,3)
Rt_H2OwCFormation.disconnect_bond(2,1)
Rt_H2OwCFormation.connect_bond(3,2)
Rt_H2OwCFormation.connect_bond(1,5)
Rt_H2OwCFormation.set_rule_name("H2OwCFormation")
Rtypelist.append(Rt_H2OwCFormation) 

Rt_H2OwOFormation = pyring.react.ReactionType()
Rt_H2OwOFormation.add_reactant_pattern(pyring.react.Substructure("O1(-H2)", pyring.react.patternsize("O1(-H2)")))
Rt_H2OwOFormation.add_reactant_pattern(pyring.react.Substructure("{Pt}5(-O3(-H4))", pyring.react.patternsize("{Pt}5(-O3(-H4))")))
Rt_H2OwOFormation.disconnect_bond(5,3)
Rt_H2OwOFormation.disconnect_bond(2,1)
Rt_H2OwOFormation.connect_bond(3,2)
Rt_H2OwOFormation.connect_bond(1,5)
Rt_H2OwOFormation.set_rule_name("H2OwOFormation")
Rtypelist.append(Rt_H2OwOFormation) 

L = pyring.react.LumpingStrategy(False)

test_network = pyring.react.RxnNetGen()
test_network.add_initial_reactants(reactantlist)
test_network.add_reaction_rules(Rtypelist)
test_network.add_lumping_strategy(L)
test_network.add_composite_sites(CompositeSiteList)
test_network.set_calc_thermo(False)
test_network.GenerateNetwork()
test_network.print_rxn_list()
 