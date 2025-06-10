__version__ = "0.0.1"

import time
from datetime import timedelta
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd

pops = [
    "AAFA",
    "AFA",
    "AFB",
    "AINDI",
    "AISC",
    "ALANAM",
    "AMIND",
    "API",
    "CARB",
    "CARHIS",
    "CARIBI",
    "CAU",
    "EURCAU",
    "FILII",
    "HAWI",
    "HIS",
    "JAPI",
    "KORI",
    "MENAFC",
    "MSWHIS",
    "NAM",
    "NCHI",
    "SCAHIS",
    "SCAMB",
    "SCSEAI",
    "VIET",
]


loci_combinations = [
    "A",
    "C",
    "B",
    "DRB1",
    "DRBX",
    "DQA1",
    "DQB1",
    "DPA1",
    "DPB1",
    "A~C~B",
    "A~B~DRB1",
    "A~C~B~DRB1",
    "A~B~DRB1~DQB1",
    "C~B~DRB1~DQB1",
    "A~C~B~DRB1~DQB1",
    "A~C~B~DRBX~DRB1~DQB1",
    "DPB1~DPA1",
    "DQA1~DQB1",
    "DRB1~DQB1",
    "DRBX~DRB1",
    "DRBX~DQB1",
    "DRBX~DRB1~DQA1~DQB1",
    "A~C~B~DRBX~DRB1~DQB1~DPB1",
    "DRBX~DRB1~DQA1~DQB1~DPB1~DPA1",
    "A~C~B~DRBX~DRB1~DQA1~DQB1~DPB1~DPA1",
]


def determine_loci_order():
    # Sample haplotype from the frequency file
    full_loci_haplotype = "A*03:01~C*07:02~B*07:02~DRB5*01:01~DRB1*15:01~DQB1*06:02~DQA1*01:02~DPA1*01:03~DPB1*04:01"  # noqa: E501
    loci = list(map(lambda a: a.split("*")[0], full_loci_haplotype.split("~")))
    drbx_index = loci.index("DRB5")
    loci_order = "~".join(loci[:drbx_index] + ["DRBX"] + loci[drbx_index + 1 :])
    return loci_order


def df_from_freq_files(freqs_folder="freqs/"):
    freq_dir = Path(freqs_folder)
    # Glob for all .freqs.gz files
    csv_gz_files = sorted(freq_dir.glob("*.freqs.gz"))
    freq_dfs = dict()
    print("Total files:", len(csv_gz_files))
    for csv_gz_file in csv_gz_files:
        population = csv_gz_file.name.removesuffix(".freqs.gz")
        print(population, csv_gz_file)
        freq_df = pd.read_csv(
            csv_gz_file,
            usecols=["Haplo", "Freq"],
            dtype={"Haplo": str, "Freq": np.float64},
            compression="gzip",
            engine="pyarrow",
        )
        # Rename the Freq column to be population
        freq_df = freq_df.rename(
            columns={"Haplo": "Haplotype", "Freq": f"{population}-Freq"}
        )
        freq_dfs[population] = freq_df

    return freq_dfs


def merge_all_with_total(freq_dfs):
    print("Merging...")
    merged_df = reduce(
        lambda left, right: pd.merge(
            left, right, on="Haplotype", how="outer", copy=False
        ),
        freq_dfs.values(),
    ).fillna(0.0)
    # Create a TotalFreq column that sums all the freq for all population
    merged_df["TotalFreq"] = merged_df.filter(like="-Freq").sum(axis=1)
    return merged_df


def keep_top_million(merged_df):
    # Keep the top million
    final_df = merged_df.sort_values(by="TotalFreq", ascending=False).head(1_000_000)
    # Rename frequency column to only be population
    final_df.columns = final_df.columns.str.replace(r"-Freq$", "", regex=True)

    return final_df


def generate_excel_names(loci_set):
    filename = loci_set + ".xlsx"
    # Max length of Excel sheet is 31
    sheet_name = loci_set[: 31 - 3] + "..." if len(loci_set) > 31 else loci_set
    return filename, sheet_name


def get_output_dir():
    output_dir = Path("./Excel")
    if not output_dir.exists():
        output_dir.mkdir()
    return output_dir


def save_to_excel(final_df, filename, sheet_name):
    # Save to an Excel File
    file_path = get_output_dir() / filename
    print("Frequencies: Saving to file:", file_path)
    final_df.to_excel(file_path, index=False, sheet_name=sheet_name)


def create_loci_combo_freq_file(final_df, loci_combo):
    loci_columns = loci_combo.split("~")
    loci_combo_df = (
        final_df.groupby(loci_columns)[pops + ["TotalFreq"]].sum().reset_index()
    )
    loci_combo_df = loci_combo_df.sort_values(by="TotalFreq", ascending=False).head(
        1_000_000
    )  # Keep top million

    if len(loci_columns) == 1:
        loci_combo_df = loci_combo_df.rename(columns={loci_combo: "Haplotype"})
    else:
        loci_combo_df.insert(
            0, "Haplotype", loci_combo_df[loci_columns].agg("~".join, axis=1)
        )
        loci_combo_df = loci_combo_df.drop(loci_columns, axis=1)

    print(loci_combo_df.head())
    filename, sheet_name = generate_excel_names(loci_combo)
    save_to_excel(loci_combo_df, filename, sheet_name)
    save_summary(loci_combo_df, loci_combo)


def save_summary(df, loci_combo):
    # Basic summary stats
    summary = df.describe()
    summary.loc["sum"] = df.sum()
    print(summary)

    filename = f"{loci_combo}-summary.xlsx"
    file_path = get_output_dir() / filename
    print("Summary: Saving to file:", file_path)
    summary.to_excel(file_path)


def create_full_locus_all_pops_df(freqs_dir):
    freq_dfs = df_from_freq_files(freqs_dir)
    merged_df = merge_all_with_total(freq_dfs)
    final_df = keep_top_million(merged_df)

    return final_df


def make_freqs(freqs_dir):
    print("Hello from nmdp-bioinformatics!")

    loci_set = determine_loci_order()

    # full_locus_freqs_df = create_full_locus_all_pops_df(freqs_dir)
    # # Save the 9-locus combined as a parquet file
    # full_locus_freqs_df.to_parquet(f"{loci_set}-haplotype.parquet", index=False)

    full_locus_freqs_df = pd.read_parquet(
        "A~C~B~DRBX~DRB1~DQA1~DQB1~DPA1~DPB1-haplotype.parquet"
    )

    full_locus_freqs_df[loci_set.split("~")] = full_locus_freqs_df[
        "Haplotype"
    ].str.split("~", expand=True)

    #
    # For all locus combo, generate Excel files
    for i, loci_combo in enumerate(loci_combinations, start=1):
        start_time = time.perf_counter()
        print(f"{i}/{len(loci_combinations)}:", loci_combo, "=>", loci_combo.split("~"))
        create_loci_combo_freq_file(full_locus_freqs_df, loci_combo)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        # Format the elapsed time
        formatted_time = str(timedelta(seconds=elapsed_time))
        print(f"Time to process {loci_combo}: {formatted_time}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog="make-freqs",
        usage="make-freqs -d <frequency directory>",
        description="Make frequency Files",
    )
    parser.add_argument("-d", "--freqs-dir", default="./freqs")
    args = parser.parse_args()
    make_freqs(args.freqs_dir)


if __name__ == "__main__":
    main()
