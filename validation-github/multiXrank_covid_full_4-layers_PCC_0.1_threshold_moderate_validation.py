import multixrank
#command for quick working example
#multixrank.Example().write(path="airport")

#creating an object for the config file
multixrank_obj = multixrank.Multixrank(config="config_full_covid_4-layers_PCC_0.1_threshold_moderate_validation.yml", wdir="/cbio/users/francis_agamah/SOFT/multiXrank")
print("multixrank_obj created")
print(multixrank_obj)

#perform random walk 
ranking_df = multixrank_obj.random_walk_rank()
print("ranking_df")
print(ranking_df)

##RANKING FOR HYPOTHESIS-DRIVEN SEED SELECTION
#ranking nodes
#multixrank_obj.write_ranking(ranking_df, path="Output_Files/results_from_4_layers_PCC_0.1_threshold_validation/moderate_hypothesis_driven/")

#multixrank_obj.to_sif(ranking_df, path="Output_Files/results_from_4_layers_PCC_0.1_threshold_validation/moderate_hypothesis_driven/seed_top3_full_4-layers_PCC_0.1_threshold_moderate_hypothesis_driven.sif", top=3)


##RANKING FOR DATA-DRIVEN SEED SELECTION
#ranking nodes
multixrank_obj.write_ranking(ranking_df, path="Output_Files/results_from_4_layers_PCC_0.1_threshold_validation/moderate_data_driven/")

multixrank_obj.to_sif(ranking_df, path="Output_Files/results_from_4_layers_PCC_0.1_threshold_validation/moderate_data_driven/seed_top3_full_4-layers_PCC_0.1_threshold_moderate_data_driven.sif", top=3)

