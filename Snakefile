import os

# -----------------------------------
# Configuration
# -----------------------------------
DIRECTORIES   = [d.strip().rstrip("/") for d in config["directories"].split(",")]
FOODB         = config.get("foodb", False)
HOST          = config.get("host", False)
METABOLOME    = config.get("metabolome", False)
E_WEIGHTS     = config.get("e_weights", False)
N_WEIGHTS     = config.get("n_weights", False)
INCLUDE_ORGS  = config.get("include_orgs", False)
ABUNDANCE_COL = config.get("abundance_col", " ") 
ALL_FOOD      = config.get("all_food", False)

with open("VERSION") as f:
    PIPELINE_VERSION = f.read().strip()

print("Running with config:")
print(f"  Directories:    {DIRECTORIES}")
print(f"  FooDB:          {FOODB}")
print(f"  Host:           {HOST}")
print(f"  Metabolome:     {METABOLOME}")
print(f"  E Weights:      {E_WEIGHTS}")
print(f"  N Weights:      {N_WEIGHTS}")
print(f"  Include Orgs:   {INCLUDE_ORGS}")
print(f"  Abundance Col:  {ABUNDANCE_COL}")
print(f"  All Food:       {ALL_FOOD}")

# -------------------------------------------------------------
# Select correct food_meta file depending on ALL_FOOD
# -------------------------------------------------------------
def select_meta_file(wildcards):
    """
    If ALL_FOOD == True, use precomputed metadata from Data/AllFood.
    Otherwise use standard output_fdborhost/food_meta.csv and require the
    CreateFoodMetadata_fdb rule to generate it.
    """
    if ALL_FOOD:
        return f"Data/AllFood/food_meta.csv"
    else:
        return f"{wildcards.dir}/output_fdb/food_meta.csv"

def select_meta_file_host(wildcards):
    if ALL_FOOD:
        return "Data/AllFood/food_meta.csv"
    else:
        return f"{wildcards.dir}/output_host/food_meta.csv"


# -----------------------------------
# Rule all – gather outputs across dirs
# -----------------------------------
rule all:
    input:
        # run info
        [f"{d}/run_info.txt" for d in DIRECTORIES],

        # FooDB outputs
        *([f"{d}/output_fdb/food_compound_report.html" for d in DIRECTORIES]
          if FOODB and not ALL_FOOD else []),
        *([f"{d}/output_fdb/microbe_compound_report.html" for d in DIRECTORIES]
          if FOODB and INCLUDE_ORGS and N_WEIGHTS else []),
        *([f"{d}/output_fdb/graph/graph_results_report.html" for d in DIRECTORIES]
          if FOODB else []),
          
        # host outputs
        *([f"{d}/output_host/food_compound_report.html" for d in DIRECTORIES]
          if HOST and not ALL_FOOD else []),
        *([f"{d}/output_host/microbe_compound_report.html" for d in DIRECTORIES]
          if HOST and INCLUDE_ORGS and N_WEIGHTS else []),
        *([f"{d}/output_host/graph/graph_results_report.html" for d in DIRECTORIES]
          if HOST else [])

# ---------------------------
# Run metadata generation
# ---------------------------
rule write_run_info:
    output:
        "{dir}/run_info.txt"
    run:
        from src.SupplementalFunctions.run_info import write_run_info
        write_run_info(PIPELINE_VERSION, config, output[0])

# ---------------------------
# FOODB rules
# ---------------------------
if FOODB:

    rule all_fdb:
        input: 
            "{dir}/output_fdb/food_compound_report.html",
            "{dir}/output_fdb/AMON_output/rn_dict.json",
            "{dir}/output_fdb/graph/M_nodes_df.csv",
            "{dir}/output_fdb/graph/M_edges_df.csv",
            "{dir}/output_fdb/graph/M_AbundanceDistribution.png",
            "{dir}/output_fdb/graph/M_FoodFrequencyDistribution.png", 
            "{dir}/output_fdb/graph/network_summary.txt",
            "{dir}/output_fdb/microbe_compound_report.html" if INCLUDE_ORGS and N_WEIGHTS else [],
            "{dir}/output_fdb/graph/graph_results.csv",
            "{dir}/output_fdb/graph/graph_results_report.html",
            "{dir}/output_fdb/MetabolomeComparison_Report.html" if METABOLOME else []


    rule CreateFoodMetadata_fdb:
        input: 
            f_file = "{dir}/foodb_foods_dataframe.csv"
        output: 
            f_meta = "{dir}/output_fdb/food_meta.csv"
        conda: "environment.yaml"
        shell:
            """
            Rscript src/foodb_proc/comp_FoodDB.R \
                --diet_file {input.f_file} \
                --content Data/Content.csv \
                --ExDes_file Data/CompoundExternalDescriptor.csv \
                --meta_o_file {output.f_meta}
            """
    if not ALL_FOOD:
        rule CreateCompoundReport_fdb:
            input: 
                f_meta = select_meta_file,
                graphs = "{dir}/output_fdb/graph/M_nodes_df.csv" 
            output: 
                report = "{dir}/output_fdb/food_compound_report.html"
            conda: "environment.yaml"
            shell:
                """
                python {workflow.basedir}/src/foodb_proc/RenderCompoundAnalysis.py \
                    --food_file {input.f_meta} \
                    --output {output.report}
                """
    
    rule PrepareAMONOutput_fdb:
        input: 
            dir="{dir}"
        output:
            touch("{dir}/output_fdb/AMON_output/.prepared")
        run:
            import os, shutil
            outdir = os.path.join(input.dir, "output_fdb", "AMON_output")
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            os.makedirs(outdir, exist_ok=True)
            # create a dummy file so Snakemake sees this as complete
            open(output[0], 'w').close()


    rule RunAMON_fdb:
        input: 
            prep="{dir}/output_fdb/AMON_output/.prepared",
            kos = "{dir}/noquote_ko.txt"
        output: 
            rn_json = "{dir}/output_fdb/AMON_output/rn_dict.json"
        conda: "environment.yaml"
        shell:
            """
            rm -rf {wildcards.dir}/output_fdb/AMON_output
            amon \
                -i {input.kos} \
                -o {wildcards.dir}/output_fdb/AMON_output \
                --save_entries
            """

    rule GraphCreation_fdb:
        input: 
            f_meta = select_meta_file,
            rn_json = "{dir}/output_fdb/AMON_output/rn_dict.json",
            m_meta = "{dir}/ko_taxonomy_abundance.csv"
        params:
            flags = " ".join(filter(None, [
                "--n_weights" if N_WEIGHTS else "",
                "--e_weights" if E_WEIGHTS else "",
                "--org" if INCLUDE_ORGS else ""
            ])),
            abundance = ABUNDANCE_COL,
            graph_dir = "{dir}/output_fdb/graph/"
        output:
            nodes = "{dir}/output_fdb/graph/M_nodes_df.csv",
            edges = "{dir}/output_fdb/graph/M_edges_df.csv",
            summary = "{dir}/output_fdb/graph/network_summary.txt"
        conda: "environment.yaml"
        shell:
            """
            mkdir -p {params.graph_dir}
            python src/foodb_proc/main_foodb.py \
                --f {input.f_meta} \
                --r {input.rn_json} \
                --m_meta {input.m_meta} \
                {params.flags} \
                --a {params.abundance} \
                --o {params.graph_dir}
            """

    if INCLUDE_ORGS and N_WEIGHTS: 
        rule MicrobeCompoundReport_fdb:
            input:
                nodes = "{dir}/output_fdb/graph/M_nodes_df.csv",
                edges = "{dir}/output_fdb/graph/M_edges_df.csv"
            output:
                report = "{dir}/output_fdb/microbe_compound_report.html"
            conda: "environment.yaml"
            shell:
                """
                python {workflow.basedir}/src/RenderCompoundAnalysis_Microbe.py \
                    --node_file {input.nodes} \
                    --edge_file {input.edges} \
                    --output {output.report}
                """
    
    rule RunGraph_fdb: 
        input: 
            nodes = "{dir}/output_fdb/graph/M_nodes_df.csv",
            edges = "{dir}/output_fdb/graph/M_edges_df.csv"
        output: 
            output = "{dir}/output_fdb/graph/graph_results.csv"
        conda: "environment.yaml"
        shell:
            """
            python src/run_graph.py \
                --n {input.nodes} \
                --e {input.edges} \
                --o {output.output}
            """
        
    rule PatternReport_fdb: 
        input: 
            graph_res = "{dir}/output_fdb/graph/graph_results.csv",
            rxn_json = "{dir}/output_fdb/AMON_output/rn_dict.json"
        output: 
            output = "{dir}/output_fdb/graph/graph_results_report.html"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/RenderGraphResults_Report.py \
                --patterns {input.graph_res} \
                --rxn_json {input.rxn_json} \
                --output {output.output}
            """
    
    rule MetabolomeComparisonReport_fdb: 
        input: 
            graph_res = "{dir}/output_fdb/graph/graph_results.csv",
            metabolome = "{dir}/metabolome.csv"
        output: 
            output = "{dir}/output_fdb/"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/RenderMetabolomeComparison.py \
                --patterns {input.graph_res} \
                --metabolome {input.metabolome} \
                --output {output.output}
            """

if HOST: 

    rule all_host:
        input: 
            "{dir}/output_host/food_compound_report.html",
            "{dir}/output_host/AMON_output/rn_dict.json",
            "{dir}/output_host/graph/nodes_df.csv",
            "{dir}/output_host/graph/edges_df.csv",
            "{dir}/output_host/graph/AbundanceDistribution.png",
            "{dir}/output_host/graph/FoodFrequencyDistribution.png", 
            "{dir}/output_host/graph/network_summary.txt",
            "{dir}/output_host/microbe_compound_report.html" if INCLUDE_ORGS and N_WEIGHTS else [],
            "{dir}/output_host/graph/graph_results.csv",
            "{dir}/output_host/graph/graph_results_report.html", 
            "{dir}/output_host/MetabolomeComparison_Report.html" if METABOLOME else []

    rule CreateFoodMetadata_host:
        input: 
            f_file = "{dir}/foodb_foods_dataframe.csv"
        output: 
            f_meta = "{dir}/output_host/food_meta.csv"
        conda: "environment.yaml"
        shell:
            """
            Rscript src/foodb_proc/comp_FoodDB.R \
                --diet_file {input.f_file} \
                --content Data/Content.csv \
                --ExDes_file Data/CompoundExternalDescriptor.csv \
                --meta_o_file {output.f_meta}
            """
    if not ALL_FOOD:
        rule CreateCompoundReport_host:
            input: 
                f_meta = select_meta_file_host,
                graphs = "{dir}/output_host/graph/nodes_df.csv" 
            output: 
                report = "{dir}/output_host/food_compound_report.html"
            conda: "environment.yaml"
            shell:
                """
                python {workflow.basedir}/src/foodb_proc/RenderCompoundAnalysis.py \
                    --food_file {input.f_meta} \
                    --output {output.report}
                """
    
    rule PrepareAMONoutput_host:
        input: 
            dir="{dir}"
        output:
            touch("{dir}/output_host/AMON_output/.prepared")
        run:
            import os, shutil
            outdir = os.path.join(input.dir, "output_host", "AMON_output")
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            os.makedirs(outdir, exist_ok=True)
            # create a dummy file so Snakemake sees this as complete
            open(output[0], 'w').close()


    rule RunAMON_host:
        input: 
            prep="{dir}/output_host/AMON_output/.prepared",
            mkos = "{dir}/noquote_ko.txt",
            hkos = "{dir}/noquote_ko_host.txt"
        output: 
            rn_json = "{dir}/output_host/AMON_output/rn_dict.json"
        conda: "environment.yaml"
        shell:
            """
            rm -rf {wildcards.dir}/output_host/AMON_output
            amon \
            --gene_set {input.mkos} \
            --output_dir {wildcards.dir}/output_host/AMON_output \
            --other_gene_set {input.hkos} \
            --gene_set_name "microbe" \
            --other_gene_set_name "host" \
            --save_entries
            """

    rule GraphCreation_host:
        input: 
            f_meta = select_meta_file_host,
            rn_json = "{dir}/output_host/AMON_output/rn_dict.json",
            m_meta = "{dir}/ko_taxonomy_abundance.csv",
            h_meta = "{dir}/host_ko_abundance.csv" 
        params:
            flags = " ".join(filter(None, [
                "--n_weights" if N_WEIGHTS else "",
                "--e_weights" if E_WEIGHTS else "",
                "--org" if INCLUDE_ORGS else ""
            ])),
            microbe_abundance = ABUNDANCE_COL, 
            graph_dir = "{dir}/output_host/graph/"
        output:
            nodes = "{dir}/output_host/graph/nodes_df.csv",
            edges = "{dir}/output_host/graph/edges_df.csv",
            summary = "{dir}/output_host/graph/network_summary.txt"
        conda: "environment.yaml"
        shell:
            """
            mkdir -p {params.graph_dir}
            python src/Host/main_host.py \
                --f {input.f_meta} \
                --r {input.rn_json} \
                --m_meta {input.m_meta} \
                --h_meta {input.h_meta} \
                {params.flags} \
                --microbe_a {params.microbe_abundance} \
                --o {params.graph_dir}
            """

    if INCLUDE_ORGS and N_WEIGHTS: 
        rule MicrobeCompoundReport_host:
            input:
                nodes = "{dir}/output_host/graph/nodes_df.csv",
                edges = "{dir}/output_host/graph/edges_df.csv"
            output:
                report = "{dir}/output_host/microbe_compound_report.html"
            conda: "environment.yaml"
            shell:
                """
                python {workflow.basedir}/src/RenderCompoundAnalysis_Microbe.py \
                    --node_file {input.nodes} \
                    --edge_file {input.edges} \
                    --output {output.report}
                """
    
    rule RunGraph_host: 
        input: 
            nodes = "{dir}/output_host/graph/nodes_df.csv",
            edges = "{dir}/output_host/graph/edges_df.csv"
        output: 
            output = "{dir}/output_host/graph/graph_results.csv"
        conda: "environment.yaml"
        shell:
            """
            python src/Host/host_run_graph.py \
                --n {input.nodes} \
                --e {input.edges} \
                --o {output.output}
            """

    rule PatternReport_host: 
        input: 
            graph_res = "{dir}/output_host/graph/graph_results.csv",
            rxn_json = "{dir}/output_host/AMON_output/rn_dict.json"
        output: 
            output = "{dir}/output_host/graph/graph_results_report.html"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/Host/RenderHostGraphResults_Report.py \
                --patterns {input.graph_res} \
                --rxn_json {input.rxn_json} \
                --output {output.output}
            """

    rule MetabolomeComparisonReport_host: 
        input: 
            graph_res = "{dir}/output_host/graph/graph_results.csv",
            metabolome = "{dir}/metabolome.csv"
        output: 
            output = "{dir}/output_host/"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/RenderMetabolomeComparison.py \
                --patterns {input.graph_res} \
                --metabolome {input.metabolome} \
                --output {output.output}
            """