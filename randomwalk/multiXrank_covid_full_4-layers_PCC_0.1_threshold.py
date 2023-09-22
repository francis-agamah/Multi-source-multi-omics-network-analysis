import multixrank
import sys
#command for quick working example
#multixrank.Example().write(path="airport")

# Docker run arguments
config_argument = sys.argv[1]


#creating an object for the config file
config_file_name = str("/random_walk/Data/config_full_covid_4-layers_PCC_0.1_threshold_" + str(config_argument) + ".yml")
multixrank_obj = multixrank.Multixrank(config=config_file_name, wdir="/random_walk/Data")


print("multixrank_obj created")
print(multixrank_obj)

#perform random walk 
ranking_df = multixrank_obj.random_walk_rank()
print("ranking_df")
print(ranking_df)

##RANKING FOR HYPOTHESIS-DRIVEN SEED SELECTION
#ranking nodes

# Define output directory name
output_dir = str("/random_walk/Data/Output_Files_" + str(config_argument) + "/results_from_4_layers_PCC_0.1_threshold/moderate_hypothesis_driven/")

# Save results to output directory
multixrank_obj.write_ranking(ranking_df, path=output_dir)

multixrank_obj.to_sif(ranking_df, path= str(output_dir + "seed_top3_full_4-layers_PCC_0.1_threshold_moderate_hypothesis_driven.sif"), top=3)


##RANKING FOR DATA-DRIVEN SEED SELECTION
#ranking nodes
#multixrank_obj.write_ranking(ranking_df, path="Output_Files/results_from_4_layers_PCC_0.1_threshold/moderate_data_driven/")

#multixrank_obj.to_sif(ranking_df, path="Output_Files/results_from_4_layers_PCC_0.1_threshold/moderate_data_driven/seed_top3_full_4-layers_PCC_0.1_threshold_moderate_data_driven.sif", top=3)

