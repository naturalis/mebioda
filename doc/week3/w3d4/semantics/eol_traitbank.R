library(traits)
library(taxize)
library(dplyr)

# example...
species <- "Cocos nucifera"

# For EoL traitbank data we need to provide the EoL taxon ID as input parameter. Hence,
# we first need to do a TNRS lookup of these, as follows:
sources <- gnr_datasources() # frame with global names sources
eol_id <- sources[sources$title == "EOL", "id"] # lookup the id of the EOL source
eol_tnrs <- gnr_resolve(species, data_source_ids = c(eol_id), fields = "all") # resolve species
eol_taxon_id <- eol_tnrs[eol_tnrs$matched_name == species,]$local_id # lookup integer id

# Now that we have the taxon id, we query the traitbank
eol_results <- traitbank(eol_taxon_id)
eol_graph <- eol_results[["graph"]] # the interesting bit in the results is the graph

# Here we select the fields with Darwin Core terms
eol_triples <- select( eol_graph, 'dwc:scientificname', 'dwc:measurementtype', 'dwc:measurementvalue', 'units' )

# Write to file
write.csv(eol_triples, file = "eol_triples.csv")
