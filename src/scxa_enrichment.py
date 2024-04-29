from pandasaurus_cxg.enrichment_analysis import AnndataEnrichmentAnalyzer
from pandasaurus_cxg.graph_generator.graph_generator import GraphGenerator


# def generate_rdf(anndata_file_path: str, author_cell_type_list: List[str], output_rdf_path: str):
def generate_rdf(anndata_file_path: str, output_rdf_path: str):
    # aea = AnndataEnrichmentAnalyzer(anndata_file_path, author_cell_type_list)
    aea = AnndataEnrichmentAnalyzer(anndata_file_path)
    aea.analyzer_manager.co_annotation_report()
    gg = GraphGenerator(aea)
    gg.generate_rdf_graph()
    # gg.set_label_adding_priority(author_cell_type_list)
    gg.add_label_to_terms()
    gg.save_rdf_graph(file_name=output_rdf_path)


if __name__ == "__main__":
    pass
