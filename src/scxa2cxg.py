import json
import os
from ftplib import FTP
from urllib import request

import anndata as ad
import pandas as pd

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/{study}/"
H5AD_EXT_FILE = ".project.h5ad"
METADATA_EXT_FILE = ".cell_metadata.tsv"


def download_files(study: str):
    """
    Download H5AD and metadata files for a study
    """
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


def convert(study: str):
    """
    Change column names to be conform to CxG schema v4.0.0
    """
    ann_data = ad.read_h5ad(f"downloads/{study}/{study}{H5AD_EXT_FILE}")
    new_obs = pd.DataFrame(data=ann_data.obs.copy(), index=ann_data.obs.index.copy())
    obs_cols = ann_data.obs_keys()

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
        },
        inplace=True
    )

    ann_data.obs = new_obs
    ad.AnnData.write_h5ad(
        ann_data,
        f"downloads/{study}/{study}_modified{H5AD_EXT_FILE}",
        compression="gzip"
    )

    return obs_cols


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
            studies.append(filename)

    ftp.quit()
    return studies


def main(chunk: int = 10):
    """
    Main function to process studies and convert them into CxG schema
    """
    studies = get_studies("E-GEO")

    obs_by_experiment = {}
    for study in studies[:chunk]:
        print(f"Downloading files for {study}")
        download_files(study)
        print("Download completed")
        obs_columns = convert(study)
        obs_by_experiment[study] = obs_columns

    with open("obs_columns_by_experiment.json", 'w', encoding='utf-8') as f:
        json.dump(obs_by_experiment, f, ensure_ascii=False, indent=4)


if __name__ == '__main__':
    main()
