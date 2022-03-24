
getwd() # check directory
# working directory is /Users/m741s365/Desktop/Github/DiNV-Dv1-Genome-Integration
# chimeric alignment position file is in the Mapping directory
# read in that file
chim_pos <- read.csv("Mapping/all_chimeric_positions_positions.csv", header=TRUE, sep=",")
# look at the dataframe
head(chim_pos)
nrow(chim_pos)
# Now I want to find all the rows where one of the chimeric positions (5th column) matchs a main position (2nd column)
matched_chim_pos <- chim_pos[chim_pos$chim_position %in% chim_pos$main_position,]
head(matched_chim_pos)
nrow(matched_chim_pos)
print(matched_chim_pos)

setdiff(chim_pos,matched_chim_pos)
# ok so this 
