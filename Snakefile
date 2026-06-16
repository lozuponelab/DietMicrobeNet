import os

# -----------------------------------
# Configuration
# -----------------------------------
DIRECTORIES     = [d.strip().rstrip("/") for d in config["directories"].split(",")]
DIETMICROBE     = config.get("dietmicrobe", False)
DIETMICROBEHOST = config.get("dietmicrobehost", False)
METABOLOME      = config.get("metabolome", False)
E_WEIGHTS       = config.get("e_weights", False)
N_WEIGHTS       = config.get("n_weights", False)
INCLUDE_ORGS    = config.get("include_orgs", False)
ABUNDANCE_COL   = config.get("abundance_col", " ") 
ALL_FOOD        = config.get("all_food", False)

with open("VERSION") as f:
    PIPELINE_VERSION = f.read().strip()

print("Running with config:")
print(f"  Directories:     {DIRECTORIES}")
print(f"  DietMicrobe:     {DIETMICROBE}")
print(f"  DietMicrobeHost: {DIETMICROBEHOST}")
print(f"  Metabolome:      {METABOLOME}")
print(f"  E Weights:       {E_WEIGHTS}")
print(f"  N Weights:       {N_WEIGHTS}")
print(f"  Include Orgs:    {INCLUDE_ORGS}")
print(f"  Abundance Col:   {ABUNDANCE_COL}")
print(f"  All Food:        {ALL_FOOD}")

# -------------------------------------------------------------
# Select correct food_meta file depending on ALL_FOOD
# -------------------------------------------------------------
def select_meta_file(wildcards):
    """
    If ALL_FOOD == True, use precomputed metadata from Data/AllFood.
    Otherwise use standard output_dmORdmh/food_meta.csv and require the
    CreateFoodMetadata rule to generate it.
    """
    if ALL_FOOD:
        return f"Data/AllFood/food_meta.csv"
    else:
        return f"{wildcards.dir}/output_dm/food_meta.csv"

def select_meta_file_host(wildcards):
    if ALL_FOOD:
        return "Data/AllFood/food_meta.csv"
    else:
        return f"{wildcards.dir}/output_dmh/food_meta.csv"


# -----------------------------------
# Rule all – gather outputs across dirs
# -----------------------------------
rule all:
    input:
        # run info
        [f"{d}/run_info.txt" for d in DIRECTORIES],

        # DIETMICROBE outputs
        *([f"{d}/output_dm/food_compound_report.html" for d in DIRECTORIES]
          if DIETMICROBE and not ALL_FOOD else []),
        *([f"{d}/output_dm/microbe_compound_report.html" for d in DIRECTORIES]
          if DIETMICROBE and INCLUDE_ORGS and N_WEIGHTS else []),
        *([f"{d}/output_dm/graph/graph_results_report.html" for d in DIRECTORIES]
          if DIETMICROBE else []),
          
        # DIETMICROBEHOST outputs
        *([f"{d}/output_dmh/food_compound_report.html" for d in DIRECTORIES]
          if DIETMICROBEHOST and not ALL_FOOD else []),
        *([f"{d}/output_dmh/microbe_compound_report.html" for d in DIRECTORIES]
          if DIETMICROBEHOST and INCLUDE_ORGS and N_WEIGHTS else []),
        *([f"{d}/output_dmh/graph/graph_results_report.html" for d in DIRECTORIES]
          if DIETMICROBEHOST else [])

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
# DIETMICROBE rules
# ---------------------------
if DIETMICROBE:

    rule all_dm:
        input: 
            "{dir}/output_dm/food_compound_report.html",
            "{dir}/output_dm/AMON_output/rn_dict.json",
            "{dir}/output_dm/graph/M_nodes_df.csv",
            "{dir}/output_dm/graph/M_edges_df.csv",
            "{dir}/output_dm/graph/M_AbundanceDistribution.png",
            "{dir}/output_dm/graph/M_FoodFrequencyDistribution.png", 
            "{dir}/output_dm/graph/network_summary.txt",
            "{dir}/output_dm/microbe_compound_report.html" if INCLUDE_ORGS and N_WEIGHTS else [],
            "{dir}/output_dm/graph/graph_results.csv",
            "{dir}/output_dm/graph/graph_results_report.html",
            "{dir}/output_dm/MetabolomeComparison_Report.html" if METABOLOME else []


    rule CreateFoodMetadata_dm:
        input: 
            f_file = "{dir}/foodb_foods_dataframe.csv"
        output: 
            f_meta = "{dir}/output_dm/food_meta.csv"
        conda: "environment.yaml"
        shell:
            """
            Rscript src/dietmicrobe/comp_FoodDB.R \
                --diet_file {input.f_file} \
                --content Data/Content.csv \
                --ExDes_file Data/CompoundExternalDescriptor.csv \
                --meta_o_file {output.f_meta}
            """

    if not ALL_FOOD:
        rule CreateCompoundReport_dm:
            input: 
                f_meta = select_meta_file,
                graphs = "{dir}/output_dm/graph/M_nodes_df.csv" 
            output: 
                report = "{dir}/output_dm/food_compound_report.html"
            conda: "environment.yaml"
            shell:
                """
                python {workflow.basedir}/src/dietmicrobe/RenderCompoundAnalysis.py \
                    --food_file {input.f_meta} \
                    --output {output.report}
                """
    
    rule PrepareAMONOutput_dm:
        input: 
            dir="{dir}"
        output:
            touch("{dir}/output_dm/AMON_output/.prepared")
        run:
            import os, shutil
            outdir = os.path.join(input.dir, "output_dm", "AMON_output")
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            os.makedirs(outdir, exist_ok=True)
            # create a dummy file so Snakemake sees this as complete
            open(output[0], 'w').close()


    rule RunAMON_dm:
        input: 
            prep="{dir}/output_dm/AMON_output/.prepared",
            kos = "{dir}/noquote_ko.txt"
        output: 
            rn_json = "{dir}/output_dm/AMON_output/rn_dict.json"
        conda: "environment.yaml"
        shell:
            """
            rm -rf {wildcards.dir}/output_dm/AMON_output
            amon \
                -i {input.kos} \
                -o {wildcards.dir}/output_dm/AMON_output \
                --save_entries
            """

    rule GraphCreation_dm:
        input: 
            f_meta = select_meta_file,
            rn_json = "{dir}/output_dm/AMON_output/rn_dict.json",
            m_meta = "{dir}/ko_taxonomy_abundance.csv"
        params:
            flags = " ".join(filter(None, [
                "--n_weights" if N_WEIGHTS else "",
                "--e_weights" if E_WEIGHTS else "",
                "--org" if INCLUDE_ORGS else ""
            ])),
            abundance = ABUNDANCE_COL,
            graph_dir = "{dir}/output_dm/graph/"
        output:
            nodes = "{dir}/output_dm/graph/M_nodes_df.csv",
            edges = "{dir}/output_dm/graph/M_edges_df.csv",
            summary = "{dir}/output_dm/graph/network_summary.txt"
        conda: "environment.yaml"
        shell:
            """
            mkdir -p {params.graph_dir}
            python src/dietmicrobe/main_foodb.py \
                --f {input.f_meta} \
                --r {input.rn_json} \
                --m_meta {input.m_meta} \
                {params.flags} \
                --a {params.abundance} \
                --o {params.graph_dir}
            """

    if INCLUDE_ORGS and N_WEIGHTS: 
        rule MicrobeCompoundReport_dm:
            input:
                nodes = "{dir}/output_dm/graph/M_nodes_df.csv",
                edges = "{dir}/output_dm/graph/M_edges_df.csv"
            output:
                report = "{dir}/output_dm/microbe_compound_report.html"
            conda: "environment.yaml"
            shell:
                """
                python {workflow.basedir}/src/RenderCompoundAnalysis_Microbe.py \
                    --node_file {input.nodes} \
                    --edge_file {input.edges} \
                    --output {output.report}
                """
    
    rule RunGraph_dm: 
        input: 
            nodes = "{dir}/output_dm/graph/M_nodes_df.csv",
            edges = "{dir}/output_dm/graph/M_edges_df.csv"
        output: 
            output = "{dir}/output_dm/graph/graph_results.csv"
        conda: "environment.yaml"
        shell:
            """
            python src/dietmicrobe/run_graph.py \
                --n {input.nodes} \
                --e {input.edges} \
                --o {output.output}
            """
        
    rule PatternReport_dm: 
        input: 
            graph_res = "{dir}/output_dm/graph/graph_results.csv",
            rxn_json = "{dir}/output_dm/AMON_output/rn_dict.json"
        output: 
            output = "{dir}/output_dm/graph/graph_results_report.html"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/dietmicrobe/RenderGraphResults_Report.py \
                --patterns {input.graph_res} \
                --rxn_json {input.rxn_json} \
                --output {output.output}
            """
    
    rule MetabolomeComparisonReport_dm: 
        input: 
            graph_res = "{dir}/output_dm/graph/graph_results.csv",
            metabolome = "{dir}/metabolome.csv"
        output: 
            output = "{dir}/output_dm/"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/RenderMetabolomeComparison.py \
                --patterns {input.graph_res} \
                --metabolome {input.metabolome} \
                --output {output.output}
            """

if DIETMICROBEHOST: 

    rule all_host:
        input: 
            "{dir}/output_dmh/food_compound_report.html",
            "{dir}/output_dmh/AMON_output/rn_dict.json",
            "{dir}/output_dmh/graph/nodes_df.csv",
            "{dir}/output_dmh/graph/edges_df.csv",
            "{dir}/output_dmh/graph/AbundanceDistribution.png",
            "{dir}/output_dmh/graph/FoodFrequencyDistribution.png", 
            "{dir}/output_dmh/graph/network_summary.txt",
            "{dir}/output_dmh/microbe_compound_report.html" if INCLUDE_ORGS and N_WEIGHTS else [],
            "{dir}/output_dmh/graph/graph_results.csv",
            "{dir}/output_dmh/graph/graph_results_report.html", 
            "{dir}/output_dmh/MetabolomeComparison_Report.html" if METABOLOME else []

    rule CreateFoodMetadata_dmh:
        input: 
            f_file = "{dir}/foodb_foods_dataframe.csv"
        output: 
            f_meta = "{dir}/output_dmh/food_meta.csv"
        conda: "environment.yaml"
        shell:
            """
            Rscript src/dietmicrobe/comp_FoodDB.R \
                --diet_file {input.f_file} \
                --content Data/Content.csv \
                --ExDes_file Data/CompoundExternalDescriptor.csv \
                --meta_o_file {output.f_meta}
            """
    if not ALL_FOOD:
        rule CreateCompoundReport_dmh:
            input: 
                f_meta = select_meta_file_host,
                graphs = "{dir}/output_dmh/graph/nodes_df.csv" 
            output: 
                report = "{dir}/output_dmh/food_compound_report.html"
            conda: "environment.yaml"
            shell:
                """
                python {workflow.basedir}/src/dietmicrobe/RenderCompoundAnalysis.py \
                    --food_file {input.f_meta} \
                    --output {output.report}
                """
    
    rule PrepareAMONoutput_dmh:
        input: 
            dir="{dir}"
        output:
            touch("{dir}/output_dmh/AMON_output/.prepared")
        run:
            import os, shutil
            outdir = os.path.join(input.dir, "output_dmh", "AMON_output")
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            os.makedirs(outdir, exist_ok=True)
            # create a dummy file so Snakemake sees this as complete
            open(output[0], 'w').close()


    rule RunAMON_host:
        input: 
            prep="{dir}/output_dmh/AMON_output/.prepared",
            mkos = "{dir}/noquote_ko.txt",
            hkos = "{dir}/noquote_ko_host.txt"
        output: 
            rn_json = "{dir}/output_dmh/AMON_output/rn_dict.json"
        conda: "environment.yaml"
        shell:
            """
            rm -rf {wildcards.dir}/output_dmh/AMON_output
            amon \
            --gene_set {input.mkos} \
            --output_dir {wildcards.dir}/output_dmh/AMON_output \
            --other_gene_set {input.hkos} \
            --gene_set_name "microbe" \
            --other_gene_set_name "host" \
            --save_entries
            """

    rule GraphCreation_host:
        input: 
            f_meta = select_meta_file_host,
            rn_json = "{dir}/output_dmh/AMON_output/rn_dict.json",
            m_meta = "{dir}/ko_taxonomy_abundance.csv",
            h_meta = "{dir}/host_ko_abundance.csv" 
        params:
            flags = " ".join(filter(None, [
                "--n_weights" if N_WEIGHTS else "",
                "--e_weights" if E_WEIGHTS else "",
                "--org" if INCLUDE_ORGS else ""
            ])),
            microbe_abundance = ABUNDANCE_COL, 
            graph_dir = "{dir}/output_dmh/graph/"
        output:
            nodes = "{dir}/output_dmh/graph/nodes_df.csv",
            edges = "{dir}/output_dmh/graph/edges_df.csv",
            summary = "{dir}/output_dmh/graph/network_summary.txt"
        conda: "environment.yaml"
        shell:
            """
            mkdir -p {params.graph_dir}
            python src/dietmicrobehost/main_host.py \
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
                nodes = "{dir}/output_dmh/graph/nodes_df.csv",
                edges = "{dir}/output_dmh/graph/edges_df.csv"
            output:
                report = "{dir}/output_dmh/microbe_compound_report.html"
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
            nodes = "{dir}/output_dmh/graph/nodes_df.csv",
            edges = "{dir}/output_dmh/graph/edges_df.csv"
        output: 
            output = "{dir}/output_dmh/graph/graph_results.csv"
        conda: "environment.yaml"
        shell:
            """
            python src/dietmicrobehost/host_run_graph.py \
                --n {input.nodes} \
                --e {input.edges} \
                --o {output.output}
            """

    rule PatternReport_host: 
        input: 
            graph_res = "{dir}/output_dmh/graph/graph_results.csv",
            rxn_json = "{dir}/output_dmh/AMON_output/rn_dict.json"
        output: 
            output = "{dir}/output_dmh/graph/graph_results_report.html"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/dietmicrobehost/RenderHostGraphResults_Report.py \
                --patterns {input.graph_res} \
                --rxn_json {input.rxn_json} \
                --output {output.output}
            """

    rule MetabolomeComparisonReport_host: 
        input: 
            graph_res = "{dir}/output_dmh/graph/graph_results.csv",
            metabolome = "{dir}/metabolome.csv"
        output: 
            output = "{dir}/output_dmh/"
        conda: "environment.yaml"
        shell: 
            """
            python {workflow.basedir}/src/RenderMetabolomeComparison.py \
                --patterns {input.graph_res} \
                --metabolome {input.metabolome} \
                --output {output.output}
            """