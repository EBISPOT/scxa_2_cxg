import json
import os
from ftplib import FTP

import anndata as ad
import httpx

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

    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    if not os.path.isdir(study_dir):
        os.mkdir(study_dir)

    if not os.path.exists(h5ad_path):
        res = httpx.get(f"{BASE_URL.format(study=study)}{study}{H5AD_EXT_FILE}")
        with open(h5ad_path, 'wb') as f:
            f.write(res.content)


def convert(study: str):
    """
    Change column names to be conform to CxG schema v4.0.0
    """
    ann_data = ad.read_h5ad(
        f"downloads/{study}/{study}{H5AD_EXT_FILE}",
        backed="r+"
    )
    new_obs = ann_data.obs.copy()
    new_var = ann_data.var.copy()
    obs_cols = ann_data.obs_keys()
    print(ann_data.var_keys())

    new_obs.rename(columns={
        "developmental_stage": "development_stage",
        "developmental_stage_ontology": "developmental_stage_ontology_term_id",
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
    duplicated_cols = [col for col in new_obs.columns if ".1" in col]
    new_obs.drop(labels=duplicated_cols, axis="columns", inplace=True)

    new_var.rename(columns={
        "type": "feature_biotype",
        "gene_name": "feature_name",
    })

    print(new_var.gene_name.drop_duplicates().tolist())
    # print(new_var.gene_version.drop_duplicates().tolist()) # not available in all experiments
    print(new_var.gene_biotype.drop_duplicates().tolist())
    print(new_var.gene_source.drop_duplicates().tolist())
    print(new_var.type.drop_duplicates().tolist())
    print(new_var.source.drop_duplicates().tolist())
    print(new_var.mito.drop_duplicates().tolist())

    ann_data.obs = new_obs
    ann_data.var = new_var
    ann_data.write_h5ad(
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
            if "E-CURD-100" not in filename:
                studies.append(filename)

    ftp.quit()
    return studies


def get_studies_downloaded(root_path):
    """
    Get the studies located in the downloads folder
    """
    for dir_name in os.listdir(root_path):
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


def main(chunk: int = 25):
    """
    Main function to process studies and convert them into CxG schema
    """
    # studies = get_studies("E-HCAD")

    obs_by_experiment = {}
    # for study in studies[:chunk]:
    for study in get_studies_downloaded("downloads"):
        print(f"Downloading files for {study}")
        download_files(study)
        print("Download completed")
        obs_columns = convert(study)
        # print(obs_columns)
        check_modified(study)
        obs_by_experiment[study] = obs_columns

    with open("obs_columns_by_experiment.json", 'w', encoding='utf-8') as f:
        json.dump(obs_by_experiment, f, ensure_ascii=False, indent=4)


if __name__ == '__main__':
    main()
