
import sys
import argparse
import os

def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python compute_scores.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="compute_scores.py",description="This program can be used to score the predicted interactors based on several parameters.")
    
    parser.add_argument("-fofn","--fofn",help="Enter a list of csv file names for each bait to be scored.",required=True)
    parser.add_argument("-threshold_p_val","--threshold_p_val",default=0.05,help="The maximum cut off for p values for spec score")
    parser.add_argument("-threshold_fold_change","--threshold_fold_change",default=2,help="Minimum fold change for spec score")
    parser.add_argument("-threshold_enrichment_score","--threshold_enrichment_score",default=0.5,help="Minimum enrichment score")
    parser.add_argument("-out_dir","--out_dir",help="Enter the name of the output directory where all the results will be stored",required=True)
    
    return parser.parse_args()

def prepFilesForScoreComputation(options):
    """
    Generates the arguments required for score computation
    """
    # Create the output directory
    os.system("mkdir -p "+options.out_dir)
    
    # Read in the fofn
    fofn=[]
    fhr=open(options.fofn,"r")
    for line in fhr:
        fofn.append(line.strip())
    fhr.close()
    
    salmon_counts=[]
    bait_names_list=[]
    num_of_replicates_list=[]
    out_dir_list=[]
    final_reports_list=[]
    salmon_counts_list=[]
    for eachbait in fofn:
        fhr=open(eachbait,"r")
        for line in fhr:
            if "output_directory" in line:
                out_dir_list.append(line.strip().split(",")[-1])
                bait_names_list.append(line.strip().split(",")[-1].split("/")[-1])
                final_reports_list.append(line.strip().split(",")[-1]+"/"+line.strip().split(",")[-1].split("/")[-1]+"_final_report.csv")
                salmon_counts_list.append(line.strip().split(",")[-1]+"/"+line.strip().split(",")[-1].split("/")[-1]+"_salmon_counts.matrix")
            elif "input_fullpath" in line:
                num_of_replicates_list.append(line.strip().split(",\"")[-1][:-1].count(";")+1)
        fhr.close()
    arguments_for_r_code=[]
    arguments_for_r_code.append(options.threshold_p_val)
    arguments_for_r_code.append(options.threshold_fold_change)
    arguments_for_r_code.append(options.threshold_enrichment_score)
    arguments_for_r_code.append(options.out_dir)
    arguments_for_r_code.append("::".join(salmon_counts_list))
    arguments_for_r_code.append("::".join(final_reports_list))
    arguments_for_r_code.append("::".join(bait_names_list))
    arguments_for_r_code.append("::".join(map(str,num_of_replicates_list)))
    return arguments_for_r_code

def runRCode(arguments_for_r_code):
    os.system("Rscript run_scores_args.R "+" ".join(map(str,arguments_for_r_code)))
    
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    arguments_for_r_code=prepFilesForScoreComputation(options)
    runRCode(arguments_for_r_code)
    
if __name__ == "__main__":
    main()