# -*- coding: utf-8 -*-
"""This script takes a reference sequence of and
some interest sequences, and computes the distances
between those sequences in NCBI assembled genomes.
Was used to show that some species have the reference
and interest sequences close to one another, while they
are distant in other species."""
import pandas as pd
import seaborn as sns
import numpy as np

# CONSTANTS
REFERENCE_PROTEIN = "UB2H"

INTEREST_PROTEINS = ["AlllpoP",
                     "lpoB3",
                     "RecA",
                     "SecY"]

PROTEINS = [REFERENCE_PROTEIN] + INTEREST_PROTEINS

FILES = ["tblastnUB2Hconsensus.txt",
         "tblastnAlllpoP.txt",
         "tblastn3lpoB.txt",
         "tblastnRecAPAO1.txt",
         "tblastnSecYPAO1.txt"]

HEADER = {"subject_id": np.str,
          "evalue": np.str,
          "subject_tax_ids": np.float,
          "subject_sci_names": np.str,
          "%%_query_coverage_per_subject": np.float,
          "s_start": np.float,
          "s_end": np.float}

DATAFRAMES = {k: pd.read_csv(v,
                             names=HEADER.keys(),
                             dtype=HEADER,
                             comment="#",
                             sep="\t") for k,v in dict(zip(PROTEINS, FILES)).items()}

# Add 'center' column to dataframes
for k,v in DATAFRAMES.items():
    v["center"] = (v["s_start"] + v["s_end"]) / 2

# make empty dataframe containers
tt_common_id = DATAFRAMES[REFERENCE_PROTEIN][["subject_id",
                         "subject_tax_ids", "subject_sci_names"]].copy()
tt_has_both = DATAFRAMES[REFERENCE_PROTEIN][["subject_tax_ids", "subject_sci_names"]].copy()
tt_different_chromosome = tt_has_both.copy()
centers = DATAFRAMES[REFERENCE_PROTEIN][["subject_id",
                    "subject_tax_ids", "subject_sci_names", "center"]].copy()

for p in INTEREST_PROTEINS:
    # does the chromosome/plasmid have both the reference protein
    # and an analog of the BLASTed protein?
    tt_common_id[p] = tt_common_id.subject_id.isin(DATAFRAMES[p].subject_id)
    
    # does the taxon have both the reference protein and an analog
    # of the BLASTed protein, either on same or different chromosomes/plasmid?
    tt_has_both[p] = tt_has_both.subject_tax_ids.isin(DATAFRAMES[p].subject_tax_ids)
    
    # group central positions of all analog proteins in the same dataset
    centers = pd.merge(centers,
                       DATAFRAMES[p][["subject_id", "center"]],
                       suffixes=("", p), how="outer", on="subject_id")

tt_same_chromosome = DATAFRAMES[REFERENCE_PROTEIN].groupby("subject_tax_ids")[["subject_sci_names"]].apply(lambda x: x.sum())  # group by taxa, keep names
distances = centers[["subject_id", "subject_tax_ids", "subject_sci_names"]].copy()
for p in INTEREST_PROTEINS:
    # does the taxon have both the reference protein and an analog
    # of the BLASTed protein on the same chromosome/plasmid?
    tt_same_chromosome[p] = tt_common_id.groupby("subject_tax_ids")[p].max()
    
    # does the taxon have both the reference protein and an analog
    # of the BLASTed protein on distinct chromosomes/plasmids?
    tt_different_chromosome[p] = ~tt_has_both.subject_tax_ids.isin(tt_same_chromosome.index)
    
    # what are the distances between the center of the reference protein
    # and that of each protein analog we want to test?
    distances[p] = abs(centers["center{}".format(p)].sub(centers["center"]))

tt_different_chromosome.drop_duplicates(subset=["subject_tax_ids"], inplace=True)
        
# testing harness
"""
for n in ["tt_common_id",
          "tt_same_chromosome",
          "tt_has_both",
          "tt_different_chromosome",
          "distances"]:
    print("{}: {}".format(n,
          tt_common_id.equals(pd.read_pickle("tt_common_id"))))
"""
