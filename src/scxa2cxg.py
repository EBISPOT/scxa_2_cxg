import json
import os
import re
from ftplib import FTP

import anndata as ad
import httpx
import pandas as pd
import curies
import numpy as np

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/{study}/"
H5AD_EXT_FILE = ".project.h5ad"
MODIFIED_H5AD_EXT_FILE = "_modified.project.h5ad"
METADATA_EXT_FILE = ".idf.txt"
SDRF_EXT_FILE = ".sdrf.txt"
ASSAY_MAPPING = {
    "smart-like": "EFO:0010184",
    "smart-seq": "EFO:0008930",
    "smart-seq2": "EFO:0008931",
    "10xV1": "EFO:0010183",
    "10xV1a": "EFO:0010183",
    "10xV1i": "EFO:0010183",
    "10xV2": "EFO:0009899",
    "10x 5' v1": "EFO:0011025",
    "10xV3": "EFO:0009922",
    "10x Ig enrichment": "EFO:0010715",
    "10x feature barcode (cell surface protein profiling)": "EFO:0030011",
    "drop-seq": "EFO:0008722",
    "seq-well": "EFO:0008919",
    "SCRB-seq": "EFO:0010004",
    "MARS-seq": "EFO:0008796",
    "CEL-seq": "EFO:0008679",
    "CEL-seq2": "EFO:0010010",
    "STRT-seq": "EFO:0008953"
}
CURIE_CONVERTER = curies.get_obo_converter()


def download_files(study: str):
    """
    Download H5AD and metadata files for a study
    """
    download_dir = os.path.join("downloads")
    study_dir = os.path.join(download_dir, study)
    h5ad_path = f"{study_dir}/{study}{H5AD_EXT_FILE}"
    metadata_path = f"{study_dir}/{study}{METADATA_EXT_FILE}"
    sdrf_path = f"{study_dir}/{study}{SDRF_EXT_FILE}"

    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    if not os.path.isdir(study_dir):
        os.mkdir(study_dir)

    if not os.path.exists(h5ad_path):
        res = httpx.get(
            f"{BASE_URL.format(study=study)}{study}{H5AD_EXT_FILE}"
        )
        with open(h5ad_path, 'wb') as f:
            f.write(res.content)

    if not os.path.exists(metadata_path):
        res = httpx.get(
            f"{BASE_URL.format(study=study)}{study}{METADATA_EXT_FILE}"
        )
        with open(metadata_path, 'wb') as f:
            f.write(res.content)

    if not os.path.exists(sdrf_path):
        res = httpx.get(
            f"{BASE_URL.format(study=study)}{study}{SDRF_EXT_FILE}"
        )
        with open(sdrf_path, 'wb') as f:
            f.write(res.content)


def compress_url(url: str):
    """
    Compress URL using OBO Converter
    """
    if str(url) == str(np.nan):
        return np.nan

    # EFO namespace is not in OBO Foundry context
    if "EFO_" in str(url):
        return url.replace("http://www.ebi.ac.uk/efo/EFO_", "EFO:")

    curie = CURIE_CONVERTER.compress(url)
    if not curie:
        return np.nan
    return curie


def convert_and_save(study: str):
    """
    Change column names to be conform to CxG schema v4.0.0
    """
    ann_data = ad.read_h5ad(
        f"downloads/{study}/{study}{H5AD_EXT_FILE}",
        backed="r+"
    )
    new_obs = ann_data.obs.copy()
    new_obs.rename(columns={
        "developmental_stage": "development_stage",
        "developmental_stage_ontology": "development_stage_ontology_term_id",
        "inferred_sex": "sex",
        "inferred_sex_ontology": "sex_ontology_term_id",
        "sex_ontology": "sex_ontology_term_id",
        "cell_type_ontology": "cell_type_ontology_term_id",
        "inferred_cell_type_-_ontology_labels": "cell_type",
        "inferred_cell_type_-_ontology_labels_ontology": "cell_type_ontology_term_id",
        "authors_cell_type_-_ontology_labels": "cell_type",
        "authors_cell_type_-_ontology_labels_ontology": "cell_type_ontology_term_id",
        "organism_ontology": "organism_ontology_term_id",
        "organism_part": "tissue",
        "organism_part_ontology": "tissue_ontology_term_id",
        "disease_ontology": "disease_ontology_term_id"
        },
        inplace=True
    )
    duplicated_cols = [
        col for col in new_obs.columns if re.search(r"[a-z]\.[1-9]", col)
    ]
    new_obs.drop(labels=duplicated_cols, axis="columns", inplace=True)

    cols_terms = [col for col in new_obs.columns if col.endswith("term_id")]
    for ont_term_col in cols_terms:
        new_obs[ont_term_col] = new_obs[ont_term_col].apply(compress_url)

    meta_sdrf = pd.read_csv(
        f"downloads/{study}/{study}{SDRF_EXT_FILE}", sep="\t"
    )
    meta_sdrf.set_index(meta_sdrf["Source Name"], inplace=True)
    new_obs["assay_ontology_term_id"] = ASSAY_MAPPING[meta_sdrf["Comment[library construction]"].iloc[0]]

    new_var = ann_data.var.copy()
    new_var["feature_is_filtered"] = False

    metadata = pd.read_csv(
        f"downloads/{study}/{study}{METADATA_EXT_FILE}", sep="\t"
    )
    metadata.set_index(metadata["MAGE-TAB Version"], inplace=True)

    new_uns = ann_data.uns.copy()
    new_uns["title"] = metadata.loc[["Investigation Title"]]["1.1"].values[0]
    new_uns["default_embedding"] = "X_umap_neighbors_n_neighbors_20"
    new_uns["citation"] = f"Publication: https://doi.org/{metadata.loc[['Publication DOI']]['1.1'].values[0]}"

    ann_data.obs = new_obs
    ann_data.var = new_var
    ann_data.uns = new_uns
    ann_data.write_h5ad(
        f"downloads/{study}/{study}_modified{H5AD_EXT_FILE}",
        compression="gzip"
    )

    return ann_data


def get_studies(filter_study: str):
    """
    Get studies in the FTP with filter
    """
    ftp = FTP("ftp.ebi.ac.uk")
    ftp.login()
    ftp.cwd("pub/databases/microarray/data/atlas/sc_experiments/")
    studies = []

    for filename in ftp.nlst():
        if filename.startswith(filter_study):
            if "E-CURD-100" not in filename:
                studies.append(filename)

    ftp.quit()
    return studies


def get_studies_downloaded(root_path, filter_study):
    """
    Get the studies located in the downloads folder
    """
    for dir_name in os.listdir(root_path):
        if dir_name.startswith(filter_study):
            yield dir_name


def check_modified(study: str):
    """
    Check if conversion worked
    """
    ann_data = ad.read_h5ad(
        f"downloads/{study}/{study}_modified{H5AD_EXT_FILE}",
        backed="r+"
    )
    if ann_data:
        print("good")
