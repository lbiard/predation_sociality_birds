# predation_sociality_birds

This repository hosts data and R code for Bliard L., Dufour P., Griesser M, Covas R. (2023). Family-living and cooperative breeding in birds are associated with the number of avian predators. *preprint*. 


## GENERAL INFORMATION

1. Title: Data and scripts from "Family-living and cooperative breeding in birds are associated with the number of avian predators".

2. Author Information:
	
        A.  Name: Louis Bliard
		Institution: University of Zurich
		Address: Winterthurerstrasse 190, 8057 Zurich, Switzerland
		Email: louis.bliard@uzh.ch / louis.bliard@evobio.eu
	
        B.  Name: Paul Dufour
		Institution: University of Montpellier, CEFE CNRS
		Address: 1919 route de Mende, 34090 Montpellier, France
	
        C.  Name: Michael Griesser
		Institution: Konstanz University
		Address: Universitätsstrasse 10, 78464 Konstanz, Germany
		Email: mgriesser@ab.mpg.de
    
        D.  Name: Rita Covas
		Institution: CIBIO-InBio, University of Porto
		Address: Campus de Vairão, Rua Padre Armando Quintas, 4485-661 Vairão, Portugal
		Email: rita.covas@cibio.up.pt
    
3. Date of data collection: NA

4. Geographic location of data collection: NA

5. Information about funding sources that supported the collection of the data: RC was funded by FCT fellowship CEECIND/03451/2018. MG was supported by a Heisenberg Grant nr. GR 4650/2-1 by the German Research Foundation DFG.


## SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC0

2. Links to publications that cite or use the data: ADD preprint / paper link

3. Links to other publicly accessible locations of the data: None

4. Links/relationships to ancillary data sets: https://doi.org/10.1371/journal.pbio.2000483

5. Was data derived from another source? Some of the data came from https://doi.org/10.1371/journal.pbio.2000483

6. Recommended citation for this dataset: Bliard L., Dufour P., Griesser M, Covas R. (XXXX). Family-living and cooperative breeding in birds are associated with the number of avian predators [Data set].



## DATA & FILE OVERVIEW

1. File List: 
- `output_ericson.nex`
- `output_hackett.nex`
- `Prum_Jetz_Cooney_9993.tree`
- `distribution_overlap.txt`
- `predation_allometry.txt`
- `cooperation_pred_final.txt`
- `estimate_predator_richness.R`
- `analysis_script.R`

Note that `output_ericson.nex` and `output_hackett.nex` are distributions of phylogenetic trees and are too large to be stored on GitHub, but they will be present in the final repository (Zenodo). These files are only needed to reproduce some results present in the Appendix.

2. Relationship between files, if important: 

The dataset `predation_allometry.txt` and `distribution_overlap.txt` were used to estimate the average number of potential predators in sympatry with each focal species, using the R script `estimate_predator_richness.R`.

The dataset `cooperation_pred_final.txt` was used alongside the phylognetic trees `Prum_Jetz_Cooney_9993.tree`, `output_ericson.nex`, and `output_hackett.nex`, to perform the analyses and produce the figures present in the article, using the R script `analysis_script.R`.

## METHODOLOGICAL INFORMATION
 
1. Methods for processing the data: Raw data

2. Instrument- or software-specific information needed to interpret the data: 
- R v.4.0.5 https://www.r-project.org/
- CmdStanR https://mc-stan.org/cmdstanr/

3. People involved with sample collection, processing, analysis and/or submission: Louis Bliard, Paul Dufour, Michael Griesser, Rita Covas

### DATA-SPECIFIC INFORMATION FOR: `distribution_overlap.txt`

1. Number of variables: 312

2. Number of cases/rows: 2988

3. Variable List: 
- "Match_jetz_birdlife" = matching focal species name btween jetz phylogeny and HBW & BirdLife International checklist.
- "Sp.Scien.jetz" = focal species name based on jetz phylogenetic tree.
- "Sp.Scien.hbw" = focal species name based on HBW & BirdLife International checklist.
- "Total_cells" = total number of grid cells occupied by the focal species (based on its geographical distribution).
- "mass" = mass of the focal species (in grams).
- Columns 6 to 307 are the names of all 302 potential predator species. Each values is the total number of grid cells shared between the geographical distributions of each focal species and each predator species.
- "average_predation_richness" = estimated average richness of potential predators across the distributional range of the species.
- "total_predation_richness" = estimated total richness of potential predators across the distributional range of the species.
- "latitude_med" = median latitude of the grid cells occupied by the focal species.
- "latitude_mean" = mean latitude of the grid cells occupied by the focal species.
- "matching_name" = name of the focal species matching the name used in the `cooperation_pred.txt` dataset.

4. Missing data codes: NA

### DATA-SPECIFIC INFORMATION FOR: `predation_allometry.txt`

1. Number of variables: 10

2. Number of cases/rows: 302

3. Variable List: 
- "jetz.name" = species name of the predator.
- "freq.numbers" = proportion of birds as part of the diet. 0=exclusive bird eater / 1=frequently eats birds / 2=occasionaly eats birds.
- "min.prey" = minimum mass of main preys (in grams).
- "max.prey" = maximum mass of main preys (in grams).
- "mass" = mass of the predator species (in grams).
- "masslog" = log of the predator species mass.
- "min.prey.log" = log of the minimum mass of main preys.
- "max.prey.log" = log of the maximum mass of main preys.
- "estimated_min_prey" = estimated mass of the smallest prey the predator species can eat (based on the predator-prey allometry).
- "estinated_max_prey" = estimated mass of the largest prey the predator species can eat (based on the predator-prey allometry).

4. Missing data codes: NA

### DATA-SPECIFIC INFORMATION FOR: `cooperation_pred_final.txt`

1. Number of variables: 22

2. Number of cases/rows: 2988

3. Variable List: 
- "tip_label" = species name.
- "fam_sys_known50" = social system ("coop_families"=cooperative breeding species / "family"=family-living species / "no_fam"=non family living species).
- "devo_mode" = developmental mode (altricial / precocial).
- "region" = biogeographic realm.
- "Prcp.Mean" = mean yearly precipitation (in mm).
- "Temp.Mean" = mean yearly temperature (in Celsius).
- "Prcp.P" = precipitation predictability during the whole year.
- "Temp.P" = temperautre predictability during the whole year.
- "Prcp.Var" = precipitation variance.
- "Temp.Var" = temperature variance.
- "Foraging" = foraging strategy.
- "mass" =
- "habitat" = habitat openess (0 = fully vegetated habitat / 100 = fully open habitat).
- "average_predation_richness" = estimated average richness of potential predators across the distributional range of the species.
- "total_predation_richness" = estimated total richness of potential predators across the distributional range of the species.
- "latitude_mean" = mean latitude of the grid cells occupied by the focal species.
- "latitude_med" = median latitude of the grid cells occupied by the focal species.
- "fam_sys_known50" = social system ("coop_families"=cooperative breeding species / "family"=family-living species / "no_fam"=non family living species).
- "cooperative_breeding" = whether the species is a cooperative breeder (0=no / 1=yes).
- "phylo" = species name, matching the names in the phylogenetic tree files.
- "comment.CB" = comment about the social system classification.
- "mov_min" = migration status based on the minimum movement of the species (i.e, partial migrants would be considered migrants).

4. Missing data codes: NA
