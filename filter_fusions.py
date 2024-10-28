import os
import subprocess
import pathlib
import glob
import pandas as pd
import argparse

# Notes
"""
Generates a tsv file containing tiered fusion candidates that have been filtered based on criteria set by Malachi Griffith, Obi Griffith, and Kelsy Cotto.

Author: Kelsy Cotto
Date: August 2024
"""


# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
# Future impovements: require the user to enter either the WB OR the list of files
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(
        description="Get FDA qc stats from various files and determine if they pass or fail."
    )

    parser.add_argument(
        "-WB",
        help="the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt",
    )

    # The name of the final results folder
    parser.add_argument(
        "-f", "--fin_results", help="Name of the final results folder in gcp immuno"
    )

    return parser.parse_args()


# home_dir = os.path.expanduser("~")
# sample = "JLF-100-021"

args = parse_arguments()
print(args)
final_result = args.fin_results
fusion_candidates_df = pd.read_csv(
    args.WB
    + final_result
    + "/rnaseq/fusioninspector_evidence/finspector.FusionInspector.fusions.tsv",
    delimiter="\t",
)

fusion_dir = f"{args.WB}{final_result}/fusion_review"
pathlib.Path(fusion_dir).mkdir(parents=True, exist_ok=True)

subprocess.run(
    f"gsutil cp gs://jlf-immuno-outputs-archive/annotation_files_for_review/CancerGeneCensus-Mar2023.tsv {fusion_dir}",
    shell=True,
    check=True,
)

##filter based on support

print(fusion_candidates_df)

filt_condition = (
    fusion_candidates_df["JunctionReadCount"]
    + fusion_candidates_df["SpanningFragCount"]
) > 5
fusion_candidates_df_filt = fusion_candidates_df.loc[filt_condition]


filt_condition = fusion_candidates_df["JunctionReadCount"] > 0
fusion_candidates_df_filt = fusion_candidates_df_filt.loc[filt_condition]
##filter readthroughs
# NOT a readthrough. Defined as:
# [Left Chr and Right Chr are different] OR
# [chromosome are the same BUT Left Strand and Right Strand are different] OR
# [chromosome and strand are the same BUT ABS(Left Pos - Right Pos) < 1,000,000] OR
# [Fusion GeneA Name OR Fusion GeneB Name matches a known fusion driver gene]

cancer_genes = pd.read_csv(f"{fusion_dir}/CancerGeneCensus-Mar2023.tsv", delimiter="\t")
fusion_genes_df = cancer_genes[
    cancer_genes["Role in Cancer"].str.contains("fusion", case=False, na=False)
]
fusion_genes = fusion_genes_df["Gene Symbol"].tolist()

# Create the OR condition
condition = (
    (fusion_candidates_df_filt["#FusionName"].str.split("--").str[0].isin(fusion_genes))
    | (
        fusion_candidates_df_filt["#FusionName"]
        .str.split("--")
        .str[1]
        .isin(fusion_genes)
    )
    | (
        fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[0]
        != fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[0]
    )
    | (
        (
            fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[0]
            == fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[0]
        )
        & (
            (
                fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[2]
                != fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[2]
            )
            | (
                abs(
                    fusion_candidates_df_filt["LeftBreakpoint"]
                    .str.split(":")
                    .str[1]
                    .astype(int)
                    - fusion_candidates_df_filt["RightBreakpoint"]
                    .str.split(":")
                    .str[1]
                    .astype(int)
                )
                > 1000000
            )
        )
    )
)

# Filter the DataFrame for Class I
filtered_df = fusion_candidates_df_filt.loc[condition]
passed_candidates = filtered_df["#FusionName"].tolist()
print(passed_candidates)

pvacfuse_file_mhc_i = glob.glob(
    f"{args.WB}{final_result}/pVACfuse/mhc_i/*.all_epitopes.aggregated.tsv"
)[0]
pvacfuse_df_mhc_i = fusion_candidates_df = pd.read_csv(
    pvacfuse_file_mhc_i,
    delimiter="\t",
)

condition = (
    pvacfuse_df_mhc_i["Gene"].str.replace("_", "--").isin(passed_candidates)
) & (pvacfuse_df_mhc_i["Num Passing Peptides"] > 0)
pvacfuse_df_mhc_i["Tier"] = pvacfuse_df_mhc_i["Tier"].astype("string")

# Set 'Status' column to 'Fusion Detected' where condition is True
pvacfuse_df_mhc_i.loc[condition, "Tier"] = "Review"

# Set 'Status' column to 'No Fusion' where condition is False
pvacfuse_df_mhc_i.loc[~condition, "Tier"] = "Poor"

pvacfuse_df_mhc_i.to_csv(
    f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_i.review.csv", index=False
)


# Filter the DataFrame for Class II
pvacfuse_file_mhc_ii = glob.glob(
    f"{args.WB}{final_result}/pVACfuse/mhc_ii/*.all_epitopes.aggregated.tsv"
)[0]
pvacfuse_df_mhc_ii = fusion_candidates_df = pd.read_csv(
    pvacfuse_file_mhc_ii,
    delimiter="\t",
)

condition = (
    pvacfuse_df_mhc_ii["Gene"].str.replace("_", "--").isin(passed_candidates)
) & (pvacfuse_df_mhc_ii["Num Passing Peptides"] > 0)
pvacfuse_df_mhc_ii["Tier"] = pvacfuse_df_mhc_ii["Tier"].astype("string")

# Set 'Status' column to 'Fusion Detected' where condition is True
pvacfuse_df_mhc_ii.loc[condition, "Tier"] = "Review"

# Set 'Status' column to 'No Fusion' where condition is False
pvacfuse_df_mhc_ii.loc[~condition, "Tier"] = "Poor"

pvacfuse_df_mhc_ii.to_csv(
    f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_ii.review.csv", index=False
)
