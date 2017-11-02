setwd("C:/CUSTOM_R_FILES2")

#SEED_LIST = seq(from = 4870, to = 5000, by=1)
#SEED_LIST = c(5000,4999,4998,4997,4996,4995,4994,4993,4992,4991,4990,4989,4988,4987,4986,4985,4984,4983,4982,4981,4980)
SEED_LIST = c(4991,4984)
TIME_LIST = read.csv("timepoints.csv",header=FALSE)
TIME_MIDPOINT = 600 #455 #360 #450 #300
TIME_IMMPOSTEXT = TIME_MIDPOINT+1
TIME_LIMIT = 1200 #1000 #600

TIME_MIDPOINT_INDEX = which(TIME_LIST==TIME_MIDPOINT)
TIME_LIMIT_INDEX = which(TIME_LIST==TIME_LIMIT)

COLLAPSE_SINGLETON_OPT = "18"
CALC_TREE_INFO_OPT = "31"

SET_SEED_OPT = "11"
ADD_NEW_TRAIT_OPT = "4"
CONTRAIT_STARTVAL = "10.0"
CONTRAIT_LOWLIM = "5.0" # "2.0" #
CONTRAIT_UPLIM = "15.0"
REP_TRUNC_OPT = "r"

RUN_EPOCH_OPT = "6"

TRAIT_EVOL_RULE_OPT = "19"
TRAIT_EVOL_SCHEME_OPT = "4"
TRAIT_EVOL_PUNC_OPT = "y"
TRAIT_EVOL_MEAN_VAL = "0"
TRAIT_EVOL_STD_DEV_VAL = "0.3"

SPECIATION_TYPE_OPT = "6"

TRAIT_EVOL_A_VAL = "5"
TRAIT_EVOL_B_VAL = "-2"
TRAIT_EVOL_C_VAL = "0.1"

BKGD_EXT_OPT = "8"
BKGD_EXT_RATE_VAL = "0.05" #"0.01" 

SAVE_FILE_OPT = "45"
#SAVE_FILENAME = paste("s",SEEDSTR,"_","t",as.character(TIME_LIST[1]),sep="")

NEON_CORE_OPT = "49"

MASSEXT_TYPE_OPT = "13" #"14" #   #
MASSEXT_TYPE_FRAC = "0" #"0.5"     #"0.9" #  #
MASSEXT_TYPE_ANS = "n"  #"y" #   

if ((MASSEXT_TYPE_OPT == "13") && (MASSEXT_TYPE_FRAC == "0"))
{
  MASSEXTTYPE = "no"}   else MASSEXTTYPE = "rand"

if ((MASSEXT_TYPE_OPT == "14") && (MASSEXT_TYPE_ANS == "y")) 
  {MASSEXTTYPE = "hitrt"} else if ((MASSEXT_TYPE_OPT == "14") && (MASSEXT_TYPE_ANS == "n")) MASSEXTTYPE = "lotrt"

SLSIZE = as.numeric(length(SEED_LIST))

#####MAIN LOOP
for (z in 1:SLSIZE)
{
  SEED_VAL = SEED_LIST[z]
  SEEDSTR = as.character(SEED_VAL)
  message(SEEDSTR)
  
  #assemble name of output file
  OUTFILENAME = paste("MESAScript","_s",SEEDSTR,"_w",MASSEXTTYPE,"ext",MASSEXT_TYPE_FRAC,"_BEXT",BKGD_EXT_RATE_VAL,".txt",sep="")
  #open the output file
  OUTFILE = file(OUTFILENAME,"w")
  
  #begin assembling the script
  #create new tree and set random seed
  writeString = "*****begin by creating new tree and set random seed"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "n"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "f"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = SET_SEED_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = SEEDSTR
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "r"
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  writeString = "*****create a new continuous trait and set value"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "d"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = ADD_NEW_TRAIT_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "c"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = CONTRAIT_STARTVAL
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "r"
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  writeString = "*****begin programming queue"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "g"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "p"
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  writeString = "*****force the 1st branching event"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = RUN_EPOCH_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "1"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "y"
  writeLines(writeString,con=OUTFILE,sep="\n")
  #first do left scheme
  writeString = TRAIT_EVOL_RULE_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = TRAIT_EVOL_SCHEME_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = TRAIT_EVOL_MEAN_VAL
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "0"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = TRAIT_EVOL_PUNC_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "n"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "n"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "r"
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  #now do right scheme
  writeString = TRAIT_EVOL_SCHEME_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = TRAIT_EVOL_MEAN_VAL
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = TRAIT_EVOL_STD_DEV_VAL
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = TRAIT_EVOL_PUNC_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "n"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "n"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "r"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "2"  #markovian speciation option
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "1.0"  #forces branching to occur
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "r"  #finish programming branch-forcing epoch
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  writeString = "*****now begin evolving tree with trait-biased evolution, with set limits"
  writeLines(writeString,con=OUTFILE,sep="\n")
  ###begin loop for writing pre-extinction epochs
  for (x in 1:TIME_MIDPOINT_INDEX) {
    writeString = "*****"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = RUN_EPOCH_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    CURR_TIME = as.character(TIME_LIST[x])
    writeString = CURR_TIME
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_RULE_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    #first do left scheme
    writeString = TRAIT_EVOL_SCHEME_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_MEAN_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "0"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_PUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_UPLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_LOWLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = REP_TRUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "r"
    writeLines(writeString,con=OUTFILE,sep="\n")
    
    #now do right scheme
    writeString = TRAIT_EVOL_SCHEME_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_MEAN_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_STD_DEV_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_PUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_UPLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_LOWLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = REP_TRUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "r"
    writeLines(writeString,con=OUTFILE,sep="\n")
    
    writeString = SPECIATION_TYPE_OPT  #pick 3-param trait-biased speciation
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_A_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_B_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_C_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = BKGD_EXT_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = BKGD_EXT_RATE_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "r"  #finish programming epoch
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = COLLAPSE_SINGLETON_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CALC_TREE_INFO_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    for (y in 1:5)
    {
      writeString = "y"
      writeLines(writeString,con=OUTFILE,sep="\n")
    }
    writeString = SAVE_FILE_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    SAVE_FILENAME = paste("s",SEEDSTR,"_","t",as.character(TIME_LIST[x]),"_",sep="")
    writeString = SAVE_FILENAME
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "n"  #make sure file is saved in Nexus format with "n"
    writeLines(writeString,con=OUTFILE,sep="\n")
  }
  
  #make tree neontological before saving
  writeString = "*****make tree neontological before saving"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = NEON_CORE_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = COLLAPSE_SINGLETON_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = CALC_TREE_INFO_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  for (y in 1:5)
  {
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
  }
  writeString = SAVE_FILE_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  SAVE_FILENAME = paste("s",SEEDSTR,"_","t",as.character(TIME_LIST[x]),"NEO_",sep="")
  writeString = SAVE_FILENAME
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "n"  #make sure file is saved in Nexus format with "n"
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  #inflict mass extinction treatment--default is RANDOM
  writeString = paste("*****now prune tree to simulate trait-based extinction, random, fraction = ", MASSEXT_TYPE_FRAC, sep="")
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = MASSEXT_TYPE_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = MASSEXT_TYPE_FRAC
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = MASSEXT_TYPE_ANS
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = NEON_CORE_OPT  #make tree neontological again
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = COLLAPSE_SINGLETON_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = CALC_TREE_INFO_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  for (y in 1:5)
  {
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
  }
  writeString = SAVE_FILE_OPT
  writeLines(writeString,con=OUTFILE,sep="\n")
  SAVE_FILENAME = paste("s",SEEDSTR,"_","t",as.character(TIME_IMMPOSTEXT),"_",sep="")
  writeString = SAVE_FILENAME
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "n"  #make sure file is saved in Nexus format with "n"
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  #now do rest of epochs for post-extinction recovery
  for (x in (TIME_MIDPOINT_INDEX+1):TIME_LIMIT_INDEX) {
    writeString = "*****"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = RUN_EPOCH_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    CURR_TIME = as.character(TIME_LIST[x])
    writeString = CURR_TIME
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_RULE_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    #first do left scheme
    writeString = TRAIT_EVOL_SCHEME_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_MEAN_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "0"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_PUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_UPLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_LOWLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = REP_TRUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "r"
    writeLines(writeString,con=OUTFILE,sep="\n")
    
    #now do right scheme
    writeString = TRAIT_EVOL_SCHEME_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_MEAN_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_STD_DEV_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_PUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_UPLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "y"
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CONTRAIT_LOWLIM
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = REP_TRUNC_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "r"
    writeLines(writeString,con=OUTFILE,sep="\n")
    
    #specify how trait evolution will influence speciation rate
    writeString = SPECIATION_TYPE_OPT  #pick 3-param trait-biased speciation
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_A_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_B_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = TRAIT_EVOL_C_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = BKGD_EXT_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = BKGD_EXT_RATE_VAL
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "r"  #finish programming epoch
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = COLLAPSE_SINGLETON_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = CALC_TREE_INFO_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    for (y in 1:5)
    {
      writeString = "y"
      writeLines(writeString,con=OUTFILE,sep="\n")
    }
    writeString = SAVE_FILE_OPT
    writeLines(writeString,con=OUTFILE,sep="\n")
    SAVE_FILENAME = paste("s",SEEDSTR,"_","t",as.character(TIME_LIST[x]),"_",sep="")
    writeString = SAVE_FILENAME
    writeLines(writeString,con=OUTFILE,sep="\n")
    writeString = "n"  #make sure file is saved in Nexus format with "n"
    writeLines(writeString,con=OUTFILE,sep="\n")
  }
  
  #now run the queue
  writeString = "*****now run the queue"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "r"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "g"
  writeLines(writeString,con=OUTFILE,sep="\n")
  SUMM_FILE_NAME = paste("s",SEEDSTR,"_w",MASSEXTTYPE,"ext",MASSEXT_TYPE_FRAC,"BEXT",BKGD_EXT_RATE_VAL,"smy.txt",sep="")
  writeString = SUMM_FILE_NAME
  writeLines(writeString,con=OUTFILE,sep="\n")
  #make sure to go back to main menu and quit MESA
  writeString = "*****"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "r"
  writeLines(writeString,con=OUTFILE,sep="\n")
  writeString = "q"
  writeLines(writeString,con=OUTFILE,sep="\n")
  
  close(OUTFILE)  

}