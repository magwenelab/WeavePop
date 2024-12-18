from shiny import reactive
from shiny.express import input, render, ui
import query_database as qdb
import io
import datetime
from Bio import SeqIO
import pandas as pd

mydb='database.db'

with ui.navset_pill(id="Database"):
    with ui.nav_panel("Home"):
        ui.h1(ui.markdown("FungalPop Database"), style="padding-top: 20px;padding-bottom: 20px;")
        ui.markdown(
            """
            This database is the result of FungalPop, a workflow for the discovery of genomic diversity in a population.   
            In summary, the analyses consisted of the following steps:
            1. Processing of the annotation of the reference genome. 
            2. Map the short reads of each sample to the corresponding reference genome using [Snippy](https://github.com/tseemann/snippy).  
                * Reference-based assemblies.
                * Called variants.
                * Depth of coverage.
            3. Annotate the assembly of each sample by lifting over the annotation of the corresponding reference genome to extract the DNA and protein sequences of each isoform of each gene using [AGAT](https://agat.readthedocs.io/en/latest/index.html).
            4. Annotate the predicted effects of the called variants using [SnpEff](https://pcingola.github.io/SnpEff/).
            5. Analyze the depth of coverage to identify copy number variants (deletions and duplications).
            
            """
        )
        ui.h1("Available data",style="padding-top: 20px;padding-bottom: 10px;")
        ui.h4("Metadata",style="padding-top: 10px;padding-bottom: 10px;")
        "Metadata of the samples, including the strain, sample ID, lineage, etc."
        ui.h4(" Reference Genomes’ Annotations",style="padding-top: 10px;padding-bottom: 10px;")
        "Table with the description of the genes in the reference genome of each lineage. Including the nested features of the genes."
        ui.h4("Sequences",style="padding-top: 10px;padding-bottom: 10px;")
        "DNA and protein sequences of each isoform of each gene in all samples."
        ui.h4("Variants",style="padding-top: 10px;padding-bottom: 10px;")
        ui.markdown("""
                    Variants (SNPs, INDELs, and MNPs) and their predicted effects. Visit [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/) to see the description of the effects and impacts.  
                    Since the variants were called comparing each sample to the reference genome of its lineage, 
                    the variants have an ID that corresponds to the lineage, and they are associated with the samples 
                    they were found in.    
                    The table contains one row per combination of variant, effect, and samples where the variant it is present in.                     
                    """)
        ui.h4("Copy Number Variants",style="padding-top: 10px;padding-bottom: 10px;")
        "Table with predicted duplicated and deleted regions in the samples."
        ui.h4("Glossary",style="padding-top: 10px;padding-bottom: 10px;")
        "Definitions of the terms used in the tables."
        ui.h1("",style="padding-top: 20px;padding-bottom: 10px;")
        
    with ui.nav_panel("Metadata"):
        ui.h1("Metadata", style="padding-top: 20px;padding-bottom: 20px;")

        with ui.card():
            ui.card_header("Download metadata table") 
            @render.download(
                label="Download",
                filename=lambda: f"metadata-{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.tsv")
            def down_metadata():
                df = qdb.get_metadata(db=mydb)
                with io.BytesIO() as buf:
                    df.to_csv(buf, index=False, sep="\t")
                    yield buf.getvalue()
                
        with ui.card():
            ui.card_header("Preview table") 
            df = reactive.value(qdb.get_metadata(db=mydb))                  
            @render.data_frame
            def show_m():
                return render.DataGrid(
                df(),
                editable=True,
                selection_mode="rows",
                filters=True,
            )

    with ui.nav_panel(" Reference Genomes’ Annotations"):
        ui.h1(" Reference Genomes’ Annotations", style="padding-top: 20px;padding-bottom: 20px;")
        with ui.layout_columns(col_widths=(6,3,3), min_height="200px"):
            with ui.navset_card_pill():
                with ui.nav_panel("Gene name"):
                    ui.input_selectize(
                        "gene_name_g",
                        "Gene names",
                        choices=qdb.list_gene_names(db = mydb),
                        multiple=True,
                        width="100%",
                    )
                with ui.nav_panel("Gene IDs"):
                    ui.input_selectize(
                        "gene_id_g",
                        "Gene IDs",
                        choices=qdb.list_gene_ids(db = mydb),  
                        multiple=True,
                        width="100%",
                    )
                with ui.nav_panel("Description"):
                    ui.input_selectize(
                        "description_g",
                        "Description",
                        choices= qdb.list_descriptions(db = mydb),
                        width="100%",
                        multiple=True,
                    )
                with ui.nav_panel("Location"):
                    ui.input_selectize(
                        "chromosomes_g",
                        "Chromosome",
                        choices= qdb.list_chromosomes(db = mydb),
                        width="100%",
                        multiple=True,
                    )
                    with ui.layout_columns(col_widths=(6,6)):
                        ui.input_numeric(
                            "start_g",
                            "Start position",
                            value = None,
                            width="100%",
                        )
                        ui.input_numeric(
                            "end_g",
                            "End position",
                            value = None,
                            width="100%",
                        )
            with ui.card():
                ui.card_header("Lineage")
                ui.input_selectize(
                        "lineage_g",
                        "Lineages",
                        choices=qdb.list_lineages(db = mydb),
                        multiple=True,
                        width="100%",
                    )
            with ui.card():
                ui.card_header("Feature") 
                ui.input_checkbox_group(
                    "feature_type_g",
                    "Feature type",
                    choices= qdb.list_feature_types(db = mydb),
                    selected= qdb.list_feature_types(db = mydb),
                    width="100%",
                )
        with ui.navset_card_pill(): 
            with ui.nav_panel("Preview table"):
                ui.input_action_button("preview_g", "Preview")
                @render.data_frame
                @reactive.event(input.preview_g)
                def show_g():
                    df = qdb.genes(db=mydb, gene_name=input.gene_name_g(), gene_id=input.gene_id_g(),
                                   chromosome=input.chromosomes_g(), start=input.start_g(),end=input.end_g(),
                                   feature_type=input.feature_type_g(), description=input.description_g(), lineage=input.lineage_g())
                    if df.shape[0] > 500:
                        return df.head(500)
                    else:
                        return df
            with ui.nav_panel("Download table"):
                @render.text()
                def check_genes():
                    try:
                        df = qdb.genes(db=mydb, gene_name=input.gene_name_g(), gene_id=input.gene_id_g(),
                                        chromosome=input.chromosomes_g(), start=input.start_g(),end=input.end_g(),
                                        feature_type=input.feature_type_g(), description=input.description_g(), lineage=input.lineage_g())
                    except Exception as e:
                        return f"Error: {e}"
                @render.download(
                    label="Download",
                    filename=lambda: f"genes-{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.tsv")
                def down_genes():
                    df = qdb.genes(db=mydb, gene_name=input.gene_name_g(), gene_id=input.gene_id_g(),
                                        chromosome=input.chromosomes_g(), start=input.start_g(),end=input.end_g(),
                                        feature_type=input.feature_type_g(), description=input.description_g(),
                                        lineage=input.lineage_g())                         
                    with io.BytesIO() as buf:
                        df.to_csv(buf, index=False, sep="\t")
                        yield buf.getvalue()
                                   
    with ui.nav_panel("Sequences"):
        ui.h1("Sequences", style="padding-top: 20px;padding-bottom: 20px;")
        list_datasets = qdb.list_datasets(db = mydb)
        if len(list_datasets) > 1:
            with ui.card():
                ui.input_checkbox_group(
                    "dataset_sq",
                    "Dataset",
                    qdb.list_datasets(db = mydb))
        with ui.layout_columns(col_widths=(6,6), min_height="200px"):
            with ui.navset_card_pill(): 
                with ui.nav_panel("Strains"):
                    @render.ui
                    def show_strains_sq():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_sq())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "strain_sq",
                                "Strains",
                                choices=qdb.list_strains(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
                with ui.nav_panel("Sample IDs"):
                    @render.ui
                    def show_samples_sq():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_sq())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "sample_sq",
                                "Sample IDs",
                                choices=qdb.list_samples(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
                with ui.nav_panel("Lineage"):
                    @render.ui
                    def show_lineage_sq():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_sq())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "lineage_sq",
                                "Lineages",
                                choices=qdb.list_lineages(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
            with ui.navset_card_pill(): 
                with ui.nav_panel("Gene names"):
                    ui.input_selectize(
                        "gene_name_sq",
                        "Gene names",
                        choices=qdb.list_gene_names(db = mydb),
                        multiple=True,
                        width="100%",
                    )
                with ui.nav_panel("Gene IDs"):
                    ui.input_selectize(
                        "gene_id_sq",
                        "Gene IDs",
                        choices=qdb.list_gene_ids(db = mydb),
                        multiple=True,
                        width="100%",
                    )   
        with ui.navset_card_pill(): 
            with ui.nav_panel("Preview sequence counts"):
                ui.input_action_button("count_seqs", "Show number of sequences")
                @render.data_frame
                @reactive.event(input.count_seqs)
                def seq_counts():
                    available_input = input.__dict__.get('_map', {}).keys()
                    d = input.dataset_sq() if 'dataset_sq' in available_input else None
                    s = input.sample_sq() if 'sample_sq' in available_input else None
                    st = input.strain_sq() if 'strain_sq' in available_input else None
                    l = input.lineage_sq() if 'lineage_sq' in available_input else None
                    
                    seqs_df = qdb.sequences(db = mydb, dataset = d, seq_type = "PROTEIN",
                        gene_id = input.gene_id_sq(), gene_name=input.gene_name_sq(),
                        sample = s, strain = st, lineage= l)

                    df_counts = pd.DataFrame({"Total samples": [seqs_df["sample"].nunique()],
                                                "Total genes": [seqs_df["gene_id"].nunique()],
                                                "Total sequences": [seqs_df.shape[0]]})
                    return df_counts

            with ui.nav_panel("Download FASTA files"):
                with ui.layout_columns(col_widths=(6,6)):
                    @render.download(
                        label="Download protein sequences",
                        filename=lambda: f"sequences-{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.faa")
                    def down_prots():
                        try:
                            available_input = input.__dict__.get('_map', {}).keys()
                            d = input.dataset_sq() if 'dataset_sq' in available_input else None
                            s = input.sample_sq() if 'sample_sq' in available_input else None
                            st = input.strain_sq() if 'strain_sq' in available_input else None
                            l = input.lineage_sq() if 'lineage_sq' in available_input else None
                            
                            seqs_df = qdb.sequences(db = mydb, dataset = d, seq_type = "PROTEIN",
                                gene_id = input.gene_id_sq(), gene_name=input.gene_name_sq(),
                                sample = s, strain = st, lineage= l)

                            seqs_records = qdb.df_to_seqrecord(seqs_df)
                            seqs_text = qdb.seqrecord_to_text(seqs_records)
                        except Exception as e:
                            with io.BytesIO() as buf:
                                buf.write(f"Error: {e}".encode())
                                yield buf.getvalue()
                            return f"Error: {e}"
                        with io.BytesIO() as buf:
                            buf.write(seqs_text.encode())
                            yield buf.getvalue()  
                    @render.download(
                        label="Download DNA sequences",
                        filename=lambda: f"sequences-{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.fna")
                    def down_dna():
                        try:
                            available_input = input.__dict__.get('_map', {}).keys()
                            d = input.dataset_sq() if 'dataset_sq' in available_input else None
                            s = input.sample_sq() if 'sample_sq' in available_input else None
                            st = input.strain_sq() if 'strain_sq' in available_input else None
                            l = input.lineage_sq() if 'lineage_sq' in available_input else None
                            
                            seqs_df = qdb.sequences(db = mydb, dataset = d, seq_type = "DNA",
                                gene_id = input.gene_id_sq(), gene_name=input.gene_name_sq(),
                                sample = s, strain = st, lineage= l)

                            seqs_records = qdb.df_to_seqrecord(seqs_df)
                            seqs_text = qdb.seqrecord_to_text(seqs_records)
                        except Exception as e:
                            with io.BytesIO() as buf:
                                buf.write(f"Error: {e}".encode())
                                yield buf.getvalue()
                            return f"Error: {e}"
                        with io.BytesIO() as buf:
                            buf.write(seqs_text.encode())
                            yield buf.getvalue() 
                
    with ui.nav_panel("Variants"):
        ui.h1("Variants and their predicted effects", style="padding-top: 20px;padding-bottom: 20px;")
        ui.markdown(
            """
            The variants of each sample were predicted relative to the reference genome 
            of the corresponding lineage so **the results are biased towards the 
            strain that was used as reference!**
            """)
        list_datasets = qdb.list_datasets(db = mydb)
        if len(list_datasets) > 1:
            with ui.card():
                ui.input_checkbox_group(
                    "dataset_v",
                    "Dataset",
                    qdb.list_datasets(db = mydb))
        
        with ui.layout_columns(col_widths=(3,4,5), min_height="300px"):
            with ui.navset_card_pill(): 
                with ui.nav_panel("Strains"):
                    @render.ui
                    def show_strains():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_v())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "strain",
                                "Strains",
                                choices=qdb.list_strains(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
                with ui.nav_panel("Sample IDs"):
                    @render.ui
                    def show_samples():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_v())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "sample",
                                "Sample IDs",
                                choices=qdb.list_samples(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
                with ui.nav_panel("Lineage"):
                    @render.ui
                    def show_lineage():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_v())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "lineage",
                                "Lineages",
                                choices=qdb.list_lineages(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
            with ui.navset_card_pill(): 
                with ui.nav_panel("Gene names"):
                    ui.input_selectize(
                        "gene_name",
                        "Gene names",
                        choices=qdb.list_gene_names(db = mydb),
                        multiple=True,
                        width="100%",
                    )
                with ui.nav_panel("Gene IDs"):
                    ui.input_selectize(
                        "gene_id",
                        "Gene IDs",
                        choices=qdb.list_gene_ids(db = mydb),
                        multiple=True,
                        width="100%",
                    )
                with ui.nav_panel("Location"):
                    ui.input_selectize(
                        "chromosomes_v",
                        "Chromosome",
                        choices= qdb.list_chromosomes(db = mydb),
                        width="100%",
                        multiple=True,
                    )
                    with ui.layout_columns(col_widths=(6,6)):
                        ui.input_numeric(
                            "start_v",
                            "Start position",
                            value = None,
                            width="100%",
                        )
                        ui.input_numeric(
                            "end_v",
                            "End position",
                            value = None,
                            width="100%",
                        )
                
            with ui.navset_card_pill(): 
                with ui.nav_panel("Impact"):
                    ui.input_checkbox_group(
                            "impact",
                            "Impact",
                            qdb.list_impacts(db = mydb),
                        )
                with ui.nav_panel("Effect types"):
                    ui.input_selectize(
                            "effect_type",
                            "Effect type",
                            choices= qdb.list_effect_types(db = mydb),
                            multiple=True,
                            width="100%",
                    )
        
        with ui.navset_card_pill(): 
            with ui.nav_panel("Preview table"):
                ui.input_action_button("show_effects_g", "Preview")
                "For large tables only the first 500 rows will be shown."
                @render.data_frame
                @reactive.event(input.show_effects_g)
                def show_input():
                    available_input = input.__dict__.get('_map', {}).keys()
                    d = input.dataset_v() if 'dataset_v' in available_input else None
                    s = input.sample() if 'sample' in available_input else None
                    st = input.strain() if 'strain' in available_input else None
                    l = input.lineage() if 'lineage' in available_input else None
                    
                    df = qdb.effects(db = mydb, dataset = d, 
                                    sample = s, strain = st, lineage = l,
                                    gene_name = input.gene_name(), gene_id = input.gene_id(),
                                    impact = input.impact(), effect_type=input.effect_type(),
                                    chromosome = input.chromosomes_v(), start = input.start_v(), end = input.end_v())
                    if df.shape[0] > 500:
                        return df.head(500)
                    else:
                        return df
            with ui.nav_panel("Download table"):
                @render.download(
                    label="Download",
                    filename=lambda: f"variants-{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.tsv")
                def down_df():
                    try:
                        available_input = input.__dict__.get('_map', {}).keys()
                        d = input.dataset_v() if 'dataset_v' in available_input else None
                        s = input.sample() if 'sample' in available_input else None
                        st = input.strain() if 'strain' in available_input else None
                        l = input.lineage() if 'lineage' in available_input else None
                        
                        df = qdb.effects(db = mydb, dataset = d, 
                                        sample = s, strain = st, lineage = l,
                                        gene_name = input.gene_name(), gene_id = input.gene_id(),
                                        impact = input.impact(), effect_type=input.effect_type(),
                                        chromosome = input.chromosomes_v(), start = input.start_v(), end = input.end_v())
                    except Exception as e:
                        with io.BytesIO() as buf:
                            buf.write(f"Error: {e}".encode())
                            yield buf.getvalue()
                        return f"Error: {e}"
                    with io.BytesIO() as buf:
                        df.to_csv(buf, index=False, sep="\t")
                        yield buf.getvalue()
                    
                
    with ui.nav_panel("Copy Number Variants"):
        ui.h1("Copy Number Variants", style="padding-top: 20px;padding-bottom: 20px;")
        list_datasets = qdb.list_datasets(db = mydb)
        if len(list_datasets) > 1:
            with ui.card():
                ui.input_checkbox_group(
                    "dataset_cnv",
                    "Dataset",
                    qdb.list_datasets(db = mydb))
        with ui.layout_columns(col_widths=(4,4,4), min_height="300px"):
            with ui.navset_card_pill(): 
                with ui.nav_panel("Strains"):
                    @render.ui
                    def show_strains_cnv():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_cnv())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "strain_cnv",
                                "Strains",
                                choices=qdb.list_strains(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
                with ui.nav_panel("Sample IDs"):
                    @render.ui
                    def show_samples_cnv():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_cnv())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "sample_cnv",
                                "Sample IDs",
                                choices=qdb.list_samples(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
                with ui.nav_panel("Lineage"):
                    @render.ui
                    def show_lineage_cnv():
                        if len(list_datasets) > 1:
                            mydataset = tuple(input.dataset_cnv())
                        else:
                            mydataset = tuple(list_datasets)
                        return ui.TagList(
                            ui.input_selectize(
                                "lineage_cnv",
                                "Lineages",
                                choices=qdb.list_lineages(db = mydb, dataset=mydataset),
                                multiple=True,
                                width="100%"))
            with ui.card():
                ui.card_header("Location") 
                ui.input_selectize(
                    "chromosomes_cnv",
                    "Chromosome",
                    choices= qdb.list_chromosomes(db = mydb),
                    width="100%",
                    multiple=True,
                )
                with ui.layout_columns(col_widths=(6,6)):
                    ui.input_numeric(
                        "start_cnv",
                        "Start position",
                        value = None,
                        width="100%",
                    )
                    ui.input_numeric(
                        "end_cnv",
                        "End position",
                        value = None,
                        width="100%",
                    )
                ui.input_slider(
                        "size_cnv",
                        "Size range",
                        min = 0,
                        max = qdb.get_cnv_max_length(db = mydb),
                        step = 500,
                        value = [0, qdb.get_cnv_max_length(db = mydb)],
                        width="100%",
                )                
            with ui.card(): 
                ui.card_header("CNV")
                ui.input_select(
                        "cnv_cnv",
                        "CNV type",
                        choices=qdb.list_cnv_types(db = mydb),
                    )
                ui.input_slider(
                        "repeats_threshold_cnv",
                        "Repeats threshold",
                        min = 0,
                        max = 1,
                        step = 0.05,
                        value = 0.5,
                        width="100%",
                )
                "The repeats threshold is the proportion of the feature that is allowed to be covered by repetitive sequences."
        with ui.navset_card_pill(): 
            with ui.nav_panel("Preview table"):
                ui.input_action_button("preview_cnv", "Preview")
                "For large tables only the first 500 rows will be shown."
                @render.data_frame
                @reactive.event(input.preview_cnv)
                def show_cnv():
                    available_input = input.__dict__.get('_map', {}).keys()
                    d = input.dataset_cnv() if 'dataset_cnv' in available_input else None
                    s = input.sample_cnv() if 'sample_cnv' in available_input else None
                    st = input.strain_cnv() if 'strain_cnv' in available_input else None
                    l = input.lineage_cnv() if 'lineage_cnv' in available_input else None
                    df = qdb.get_cnv(db=mydb, dataset= d,
                                     sample= s, strain= st, lineage= l,
                                     chromosome=input.chromosomes_cnv(), start=input.start_cnv(),end=input.end_cnv(), 
                                     min_size=input.size_cnv()[0], max_size=input.size_cnv()[1], 
                                     cnv=input.cnv_cnv(),repeat_fraction=input.repeats_threshold_cnv())
                    if df.shape[0] > 500:
                        return df.head(500)
                    else:
                        return df
            with ui.nav_panel("Download table"):
                @render.download(
                    label="Download",
                    filename=lambda: f"cnv-{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.tsv")
                def down_cnv():
                    try:
                        available_input = input.__dict__.get('_map', {}).keys()
                        d = input.dataset_cnv() if 'dataset_cnv' in available_input else None
                        s = input.sample_cnv() if 'sample_cnv' in available_input else None
                        st = input.strain_cnv() if 'strain_cnv' in available_input else None
                        l = input.lineage_cnv() if 'lineage_cnv' in available_input else None
                        df = qdb.get_cnv(db=mydb,dataset= d,
                            sample= s, strain= st, lineage= l,
                            chromosome=input.chromosomes_cnv(), start=input.start_cnv(),end=input.end_cnv(), 
                            min_size=input.size_cnv()[0], max_size=input.size_cnv()[1], 
                            cnv=input.cnv_cnv(),repeat_fraction=input.repeats_threshold_cnv())
                    except Exception as e:
                        with io.BytesIO() as buf:
                            buf.write(f"Error: {e}".encode())
                            yield buf.getvalue()
                        return f"Error: {e}"
                    with io.BytesIO() as buf:
                        df.to_csv(buf, index=False, sep="\t")
                        yield buf.getvalue()
                        
                        
    with ui.nav_panel("Glossary"):
        ui.h1("Glossary", style="padding-top: 20px;padding-bottom: 20px;")
        ui.h3("Common")
        ui.markdown(
            """
            **Strain / strain**: The strain name of the samples.  
            **Sample ID / sample**: The unique identifier of the sample. It's the same as the SRA sample accession.  
            **Lineage/ lineage**: The lineage of the sample.  
            **Gene name / gene_name**: The name of the gene if available, if not, the gene ID is used.  
            **Gene ID / gene_id**: The unique identifier of the gene.  
            **Location**: Coordinates of a feature (gene features, variants or copy-number variants) in a chromosome.  
            **Chromosome / chromosome**: The chromosome where the feature is located.    
            **Start position / start**: The start position of the feature in the chromosome (first position in each chromosome is 1).    
            **End position / end**: The end position of the feature in the chromosome.    
            """
        )
        ui.h3("Genes")
        ui.markdown(
            """
            **Description / description**: The description of the gene product.  
            **Feature type / primary_tag**: The type of annotated gene feature, (i.e. gene, mRNA, exon, CDS, three_prime_UTR, five_prime_UTR). Genes with more than one isoform will have multiple mRNA features.   
            **feature_id**: The unique identifier of the feature, it includes the gene ID, the type of feature and number within the gene.  
            **parent**: The gene ID of the gene that the feature is part of, empty for the feature type gene.  
            **identical_to_main_ref**: In the mRNA feature type, if the isoform sequence of the lineage specific reference genome is identical to the main reference sequence the value is Yes, otherwise No.  
            **start_top_mutation**: In the mRNA feature type, if the isoform sequence of the reference genome has a missing stop or start codon or a inframe stop codon, compared to the sequence fo the main reference the value will be one of missing_stop_codon, missing_start_codon or inframe_stop_codon.  
            """
        )
        ui.h3("Variants")
        ui.markdown(
            """
            **Impact / impact**: The SnpEff putative impact of a variant, one of HIGH, MODERATE, LOW, MODIFIER. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **Effect type / effect_type**: The type of effect the variant has, each effect type is associated with an impact. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **var_id**: The unique identifier of the variant. It includes the lineage and a number.   
            **position**: The position of the variant in the chromosome.  
            **reference**: The allele in the reference genome.   
            **alternative**: The alternative allele.  
            **transcript_id**: The unique identifier of the transcript. It is the same as the feature_id of mRNA features in the Genes table.   
            **effect**: The functional class of the effect type, one of SILENT, MISSENSE, NONSENSE or none. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **codon_change**: Old/new codon sequences or distance to transcript. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **amino_acid_change**: Old amino acid, position, and new amino acid. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **amino_acid_length**: Length of the protein in amino acids. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **transcript_biotype**: Transcript bioType, if available. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **gene_coding**: CODING or NON_CODING classification. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **exon_rank**: Exon number in the gene. Check [SnpEff](https://pcingola.github.io/SnpEff/snpeff/inputoutput/).  
            **mean_depth_normalized**: The mean depth of coverage of the gene normalized by the genome-wide depth of coverage. It is a proxy for the copy number of the gene.  
            **mean_mapq**: The mean mapping quality of the reads that cover the gene.  
            """)
        ui.h3("Copy number variants")
        ui.markdown(
            """
            **Feature type / CNV**: The type of copy number variant (deletion or duplication).  
            **Size range**: The range of sizes of the features to include in the table.  
            **region_size**: The size of the feature in base pairs.  
            **Repeats threshold**: Parameter to filter out features with a repeat_fraction larger than this threshold.  
            **repeat_fraction**: The proportion of the feature that is covered by repetitive sequences.  
            """)
