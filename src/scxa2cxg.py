import os
import re
from ftplib import FTP

import anndata as ad
import curies
import httpx
import numpy as np
import pandas as pd

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/{study}/"
H5AD_EXT_FILE = ".project.h5ad"
MODIFIED_H5AD_EXT_FILE = "_modified.project.h5ad"
METADATA_EXT_FILE = ".idf.txt"
SDRF_EXT_FILE = ".sdrf.txt"
CLUSTER_EXT_FILE = ".clusters.tsv"
ASSAY_MAPPING = {
    "smart-like": "EFO:0010184",
    "smart-seq": "EFO:0008930",
    "smart-seq2": "EFO:0008931",
    "10xv1": "EFO:0010183",
    "10xv1a": "EFO:0010183",
    "10xv1i": "EFO:0010183",
    "10xv2": "EFO:0009899",
    "10x5prime": "EFO:0030004",
    "10x 5' v1": "EFO:0011025",
    "10xv3": "EFO:0009922",
    "10x ig enrichment": "EFO:0010715",
    "10x feature barcode (cell surface protein profiling)": "EFO:0030011",
    "drop-seq": "EFO:0008722",
    "seq-well": "EFO:0008919",
    "scrb-seq": "EFO:0010004",
    "mars-seq": "EFO:0008796",
    "cel-seq": "EFO:0008679",
    "cel-seq2": "EFO:0010010",
    "strt-seq": "EFO:0008953"
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
    cluster_path = f"{study_dir}/{study}{CLUSTER_EXT_FILE}"

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

    if not os.path.exists(cluster_path):
        res = httpx.get(
            f"{BASE_URL.format(study=study)}{study}{CLUSTER_EXT_FILE}"
        )
        with open(cluster_path, 'wb') as f:
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


def get_cluster(cell_id: str, cluster_df: pd.DataFrame):
    """
    Find cluster for a cell
    """
    sel = cluster_df[cell_id]
    if sel.any():
        return f"Cluster {sel.values[0]}"
    return ""


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
        "organism_ontology": "organism_ontology_term_id",
        "organism_part": "tissue",
        "organism_part_ontology": "tissue_ontology_term_id",
        "disease_ontology": "disease_ontology_term_id",
        "individual": "donor_id"
        },
        inplace=True
    )
    duplicated_cols = [
        col for col in new_obs.columns if re.search(r"[a-z]\.[1-9]", col)
    ]
    new_obs.drop(labels=duplicated_cols, axis="columns", inplace=True)
    if (
        "inferred_cell_type_-_ontology_labels" in new_obs.columns and
        "cell_type" not in new_obs.columns
    ):
        new_obs.rename(columns={
            "inferred_cell_type_-_ontology_labels": "cell_type",
            "inferred_cell_type_-_ontology_labels_ontology":
                "cell_type_ontology_term_id",
            },
            inplace=True
        )
    else:
        new_obs.rename(columns={
            "inferred_cell_type_-_ontology_labels": "inferred_cell_type",
            "inferred_cell_type_-_ontology_labels_ontology":
                "inferred_cell_type_ontology_term_id",
            },
            inplace=True
        )
    if (
        "authors_cell_type_-_ontology_labels" in new_obs.columns and
        "cell_type" not in new_obs.columns and
        "inferred_cell_type_-_ontology_labels" not in new_obs.columns
    ):
        new_obs.rename(columns={
            "authors_cell_type_-_ontology_labels": "cell_type",
            "authors_cell_type_-_ontology_labels_ontology":
                "cell_type_ontology_term_id",
            },
            inplace=True
        )
    else:
        new_obs.rename(columns={
            "authors_cell_type_-_ontology_labels":
                "authors_cell_type_ontology_label"
            },
            inplace=True
        )
    cols_terms = [col for col in new_obs.columns if col.endswith("term_id")]
    for ont_term_col in cols_terms:
        new_obs[ont_term_col] = new_obs[ont_term_col].apply(compress_url)

    meta_sdrf = pd.read_csv(
        f"downloads/{study}/{study}{SDRF_EXT_FILE}", sep="\t"
    )
    meta_sdrf.columns = meta_sdrf.columns.str.replace(' [', '[')
    if "donor_id" in new_obs.columns:
        # This code changes the values in assay_ontology_term of new_obs to the
        # corresponding values in Comment[library construction] of meta_sdrf
        # where Characteristics[individual] of meta_sdrf matches donor of new_obs.
        new_obs["assay_ontology_term_id"] = new_obs["donor_id"].map(
            meta_sdrf.drop_duplicates(subset="Characteristics[individual]")
            .set_index("Characteristics[individual]")["Comment[library construction]"].map(
                lambda x: ASSAY_MAPPING[x.lower()]
            )
        )
    else:
        new_obs["assay_ontology_term_id"] = ASSAY_MAPPING[
            meta_sdrf["Comment[library construction]"].iloc[0].lower()
        ]
    new_obs["suspension_type"] = meta_sdrf["Material Type"].iloc[0].lower()

    cluster_df = pd.read_csv(
        f"downloads/{study}/{study}{CLUSTER_EXT_FILE}", sep="\t"
    )
    default_clusters = cluster_df[cluster_df["sel.K"] == True]
    new_obs["cluster_nb"] = new_obs.index.map(lambda x: get_cluster(x, default_clusters))

    if "cell_type" not in new_obs.columns:
        new_obs["cell_type"] = "cell"
        new_obs["cell_type_ontology_term_id"] = "CL:0000000"

    new_var = ann_data.var.copy()
    new_var["feature_is_filtered"] = False

    metadata = pd.read_csv(
        f"downloads/{study}/{study}{METADATA_EXT_FILE}", sep="\t", header=None
    )
    new_uns = ann_data.uns.copy()
    new_uns["title"] = metadata[metadata[0].fillna('').str.startswith("Investigation Title")].values[0][0]
    new_uns["default_embedding"] = "X_umap_neighbors_n_neighbors_20"
    if not metadata[metadata[0].fillna('').str.startswith('Publication DOI')].empty:
        new_uns["citation"] = f"Publication: https://doi.org/{metadata[metadata[0].fillna('').str.startswith('Publication DOI')].values[0][0]}"
    new_uns["dataset_curie"] = f"SCXA:{study}"
    new_uns["schema_reference"] = "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md"

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
        return ann_data
    return None
