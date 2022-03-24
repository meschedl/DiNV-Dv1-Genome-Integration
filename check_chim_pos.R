## Checking to see if all chimeric positions match to a main alignment position
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
# ok so these are the same

# try this same analysis for a different file which had only readnames that I thought could be chimeric

m_chim_pos <- read.csv("Mapping/maybe_chimeric_positions.csv", header=TRUE, sep=",")
head(m_chim_pos)
matched_m_chim_pos <- m_chim_pos[m_chim_pos$chim_position %in% m_chim_pos$main_position,]
# ok so these have the same number of lines so they should be the same files basically 
setdiff(m_chim_pos,matched_m_chim_pos)
# yes they are the same. So all of the chimeric positions match an alignment position 
