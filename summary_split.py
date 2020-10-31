"""This script splits the variants in summary statistics into two list:
1. Genome-wide significant variant list.
2. The rest of the variants (PRS) list.
Variants in the genome-wide significant list are LD clumped.
Then the variants in the PRS list are clumped with the LD free genome-wide list to
remove correlation between variants in different lists.
The script takes two inputs:
1) summary statistics file in .giant2 format
2) Reference genetic data for LD clumping in plink format (.bed, .bim, .fam file trio).
"""

import argparse
import pandas as pd
import subprocess
import os


def summary_ld_clumping(summary_file, genotype_data):
    #configuring paths
    working_dir = os.path.dirname(os.path.realpath(summary_file))
    temp_dir = os.path.join(working_dir, "summary_split_temp_folder")
    summary_name = os.path.basename(os.path.realpath(summary_file))
    name_split = summary_name.split(".")
    init_genome_wide_snps_name = ".".join(name_split[:-1]) + ".gw_init." + name_split[-1]
    init_genome_wide_snps_path = os.path.join(temp_dir, init_genome_wide_snps_name)
    genome_wide_snp_clumping_path = os.path.join(temp_dir, "genome_wide_snp_clumping")
    genome_wide_prs_clumping_path = os.path.join(temp_dir, "genome_wide_prs_clumping")
    merged_list_path = os.path.join(temp_dir, "merged_list")
    final_prs_snps_name = ".".join(name_split[:-1]) + ".PRS." + name_split[-1]
    final_genome_wide_snps_name = ".".join(name_split[:-1]) + ".genome_wide." + name_split[-1]
    final_prs_snps_path = os.path.join(working_dir, final_prs_snps_name)
    final_genome_wide_snps_path = os.path.join(working_dir, final_genome_wide_snps_name)
    #subprocess commands
    cmd1 = ["mkdir", temp_dir]
    cmd2 = [
        "plink",
        "--noweb",
        "--bfile", genotype_data,
        "--clump", init_genome_wide_snps_path,
        "--clump-field", "P.2gc",
        "--clump-p1", "1",
        "--clump-p2", "1",
        "--clump-r2", "0.1",
        "--clump-kb", "500",
        "--clump-snp-field", "MarkerName",
        "--out", genome_wide_snp_clumping_path
    ]
    cmd3 = [
        "plink",
        "--noweb",
        "--bfile", genotype_data,
        "--clump", merged_list_path,
        "--clump-field", "P.2gc",
        "--clump-p1", "1e-8",
        "--clump-p2", "1",
        "--clump-r2", "0.1",
        "--clump-kb", "500",
        "--clump-snp-field", "MarkerName",
        "--out", genome_wide_prs_clumping_path
    ]
    cmd4 = ["rm", "-r", temp_dir]
    #creating folder for temporary files
    p = subprocess.run(cmd1)
    #splitting summary file into genome wide significant and other (PRS) SNPs
    summary = pd.read_csv(summary_file, delim_whitespace=True)
    genome_wide_snps = summary[summary["P.2gc"] <= 1*10**(-8)]
    prs_snps = summary[summary["P.2gc"] > 1*10**(-8)]
    genome_wide_snps.to_csv(init_genome_wide_snps_path, sep=" ", index=False, float_format="%g")
    #removing SNPs in LD in genome-wide set
    p = subprocess.run(cmd2)
    genome_wide_snp_clumping = pd.read_csv(genome_wide_snp_clumping_path + ".clumped", delim_whitespace=True)
    genome_wide_snps_list = genome_wide_snp_clumping["SNP"].values
    genome_wide_snps = genome_wide_snps[genome_wide_snps["MarkerName"].isin(genome_wide_snps_list)]
    #removing SNPs from the PRS list in LD with genome-wide SNPs
    merged_list = genome_wide_snps.append(prs_snps)
    merged_list.to_csv(merged_list_path, sep=" ", index=False, float_format="%g")
    p = subprocess.run(cmd3)
    genome_wide_prs_clumping = pd.read_csv(genome_wide_prs_clumping_path + ".clumped", delim_whitespace=True)
    removed_snps = genome_wide_prs_clumping["SP2"].values
    removed_snps_list = []
    for snps in removed_snps:
        snps = snps.replace("(1)", "")
        snps = snps.split(",")
        for snp in snps:
            removed_snps_list.append(snp)
    prs_snps = prs_snps[~prs_snps["MarkerName"].isin(removed_snps_list)]
    prs_snps.to_csv(final_prs_snps_path, sep=" ", index=False, float_format="%g")
    genome_wide_snps.to_csv(final_genome_wide_snps_path, sep=" ", index=False, float_format="%g")
    #removing intermediate files
    p = subprocess.run(cmd4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--summary", required=True, help="Summary statistics file")
    parser.add_argument("-g", "--genotype_data", required=True,
                        help="Name of the plink reference genotype file trio - .bed, .bim, .fam")
    args = parser.parse_args()
    summary_ld_clumping(summary_file=args.summary,genotype_data=args.genotype_data)