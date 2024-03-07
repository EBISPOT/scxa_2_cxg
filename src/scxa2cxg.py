import os
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd
import curies

import urllib.request as request

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/{study}/"
H5AD_EXT_FILE = ".project.h5ad"
METADATA_EXT_FILE = ".cell_metadata.tsv"

def download_files(study: str):
    download_dir = os.path.join("downloads")
    study_dir = os.path.join(download_dir, study)
    h5ad_path = f"{study_dir}/{study}{H5AD_EXT_FILE}"
    metadata_path = f"{study_dir}/{study}{METADATA_EXT_FILE}"

    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    if not os.path.isdir(study_dir):
        os.mkdir(study_dir)

    if not os.path.exists(h5ad_path):
        request.urlretrieve(
            f"{BASE_URL.format(study=study)}{study}{H5AD_EXT_FILE}",
            h5ad_path
        )

    if not os.path.exists(metadata_path):
        request.urlretrieve(
            f"{BASE_URL.format(study=study)}{study}{METADATA_EXT_FILE}",
            metadata_path
        )


def convert(study):
    download_files(study)

    ann_data = ad.read_h5ad(f"downloads/{study}/{study}{H5AD_EXT_FILE}")
    new_obs = pd.DataFrame(index=ann_data.obs.index.copy())
    new_obs.rename(columns={
        "developmental_stage": "development_stage",
        "developmental_stage_ontology": "developmental_stage_ontology_term_id",
        "inferred_sex": "sex",
        "inferred_sex_ontology": "sex_ontology_term_id",
        "inferred_cell_type_-_ontology_labels": "cell_type",
        "inferred_cell_type_-_ontology_labels_ontology": "cell_type_ontology_term_id",
        "organism_ontology": "organism_ontology_term_id",
        "organism_part": "tissue",
        "organism_part_ontology": "tissue_ontology_term_id"
    }
    )
    ann_data.obs = new_obs
    ad.AnnData.write_h5ad(ann_data, f"downloads/{study}/{study}_modified{H5AD_EXT_FILE}", compression="gzip")


convert("E-CURD-134")
