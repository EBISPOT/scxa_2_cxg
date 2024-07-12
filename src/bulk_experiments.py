"""
Script to bulk process studies and convert them into RDF after Enrichment
"""
import argparse
import itertools
import logging
from typing import List

from pandasaurus_cxg.enrichment_analysis import AnndataEnrichmentAnalyzer
from pandasaurus_cxg.graph_generator.graph_generator import GraphGenerator

from scxa2cxg import (MODIFIED_H5AD_EXT_FILE, convert_and_save, download_files,
                      get_studies, get_studies_downloaded)

logging.basicConfig(level=logging.INFO)

def generate_rdf(anndata_file_path: str, author_cell_type_list: List[str], output_rdf_path: str):
    """
    Generate RDF graph based on the given Anndata file and author cell type list.

    Args:
        anndata_file_path (str): The file path of the Anndata file.
        author_cell_type_list (List[str]): A list of author cell types.
        output_rdf_path (str): The file path to save the RDF graph.

    Returns:
        None
    """
    aea = AnndataEnrichmentAnalyzer(anndata_file_path, author_cell_type_list)
    aea.analyzer_manager.co_annotation_report()
    gg = GraphGenerator(aea)
    gg.generate_rdf_graph()
    gg.set_label_adding_priority(author_cell_type_list)
    gg.add_label_to_terms()

    metadata_field_list = ["tissue", "disease", "development_stage", "organism"]
    for field_name in metadata_field_list:
        if field_name in aea.enricher_manager.anndata.obs.columns:
            continue
        metadata_field_list.remove(field_name)
    gg.add_metadata_nodes(metadata_fields=metadata_field_list)
    gg.save_rdf_graph(file_name=output_rdf_path)

def bulk_process(study_filter: str, chunk: int, download: bool, modified: bool, exp_list: List[str]):
    """
    Process bulk experiments.

    Args:
        study_filter (str): Filter to select specific studies.
        chunk (int): Number of studies to process at a time.
        download_files (bool): Flag indicating whether to download files.

    Returns:
        None
    """
    # if download:
    #     # Get the studies in FTP containing the filter to download
    #     studies = get_studies(study_filter)
    # else:
    #     # Get the studies located in the downloads folder
    #     studies = get_studies_downloaded("downloads", study_filter)

    author_cell_type_list = [
        "cluster_nb",
        "inferred_cell_type",
        "authors_cell_type",
        "authors_cell_type_ontology_label"
    ]
    studies_iter = iter(exp_list)
    while True:
        chunk_of_studies = list(itertools.islice(studies_iter, chunk))
        if not chunk_of_studies:
            break

        for study in chunk_of_studies:
            logging.info("Processing study: %s", study)
            study_path = f"downloads/{study}/{study}{MODIFIED_H5AD_EXT_FILE}"
            output_path = f"downloads/{study}/{study}"
            if download:
                logging.info("Downloading files...")
                download_files(study)
                logging.info("Files downloaded")
            if not modified:
                logging.info("Converting and saving...")
                convert_and_save(study)
            logging.info("Generating RDF...")
            generate_rdf(study_path, author_cell_type_list, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to bulk process studies and convert them into RDF after Enrichment"
    )
    parser.add_argument("--study_filter", type=str, help="Name of the study to filter")
    parser.add_argument("--chunk_size", type=int, help="Chunk size of studies to process")
    parser.add_argument(
        "--download", action="store_true", help="Flag to indicate whether to download files or not"
    )
    parser.add_argument("--modified", action="store_true", help="Flag to indicate whether to use modified files")
    args = parser.parse_args()

    exp_list = ["E-MTAB-8698",
                "E-MTAB-9444",
                "E-MTAB-10519",
                "E-MTAB-10628",
                "E-MTAB-7194",
                "E-MTAB-7195",
                "E-GEOD-141273",
                "E-GEOD-100058",
                "E-GEOD-134722",
                "E-GEOD-136162",
                "E-GEOD-146040",
                "E-GEOD-172231",
                "E-GEOD-152495",
                "E-GEOD-141807",
                "E-GEOD-157775",
                "E-GEOD-103771",
                "E-GEOD-126139",
                "E-GEOD-157202",
                "E-GEOD-125948",
                "E-GEOD-125948",
                "E-GEOD-147601",
                "E-CURD-21",
                "E-CURD-134",
                "E-CURD-91",
                "E-CURD-90",
                "E-CURD-92",
                "E-CURD-56",
                "E-CURD-124",
                "E-CURD-87"]
    bulk_process(args.study_filter, args.chunk_size, args.download, args.modified, exp_list)
