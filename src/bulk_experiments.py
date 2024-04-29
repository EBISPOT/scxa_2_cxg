"""
Script to bulk process studies and convert them into RDF after Enrichment
"""
from scxa_enrichment import generate_rdf
from scxa2cxg import get_studies_downloaded, H5AD_EXT_FILE


# studies = get_studies("E-HCAD")

# cols_by_experiment = {}
# for study in studies[:chunk]:
for study in get_studies_downloaded("downloads", "E-GEOD-103771"):
    study_path = f"downloads/{study}/{study}_modified{H5AD_EXT_FILE}"
    output_path = f"downloads/{study}/{study}.owl"
    generate_rdf(study_path, output_path)
    # experiment = {}
    # print(f"Downloading files for {study}")
    # download_files(study)
    # print("Download completed")
    # conv_anndata = convert(study)
    # print(obs_columns)
    # check_modified(study)
    # experiment["obs_columns"] = conv_anndata.obs_keys()
    # experiment["uns_columns"] = conv_anndata.uns_keys()
    # experiment["var_columns"] = conv_anndata.var_keys()
    # experiment["obsm_columns"] = conv_anndata.obsm_keys()
    # cols_by_experiment[study] = experiment

# with open("obs_columns_by_experiment.json", 'w', encoding='utf-8') as f:
#     json.dump(cols_by_experiment, f, ensure_ascii=False, indent=4)