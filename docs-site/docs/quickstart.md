# Quickstart 

## Prerequisites 

- python>=3.10
- r-base>=4.4.2
- conda 

This program has been tested on Mac M1 and Ubuntu/linux 

## Installation 

```
git clone https://github.com/lozuponelab/DietMicrobeNet.        # clone repo
cd DietMicrobeNet                                               # move into this project directory
conda env create -f DMnet_env.yaml                              # create environment
conda activate DietMicrobeNet                                   # activate environment 
pip install -e .                                                # set up directory structure 

wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/ESXx7vpypQFOt4iVv6x-ErkBykpAVS1fppQjYZkrxkDnAA?download=1' -O Data/CompoundExternalDescriptor.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EYJUYQWmY9VDlYZIAXpzpvEBzhrnViFZQjrikXIla_aPPg?download=1' -O Data/Content.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EXyRAlYs1htNlcwz5T67BxQBGO7HfOjmfIBlkOydM0BIAw?download=1' -O Data/Food.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EbY2fD3JTcNLomKFqQhY5jABAXN-60A80PmkngRynazocg?download=1' -O Data/hmdb.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EZ1pyHd616RFkR9zG6kenuoBhZDroHYTbaGmEfwpxFOHLg?download=1' -O Data/AllFood/food_meta.csv
```

## Inputs

Your sample directory must contain the following files:

| File | Description | Required |
|------|-------------| -------- |
| `ko_taxonomy_abundance.csv` | Microbiome KO abundances and taxonomy | ✅ |
| `noquote_ko.txt` | KO list without quotes | ✅ |
| `foodb_foods_dataframe.csv` | Diet data for FooDB-based analysis | ❌  only if **foodb** used and not **all foods**|
| `host_ko_abundance.csv` | Host KO abundances | ❌  only if **host** is used|
| `noquote_ko_host.txt` | KO list without quotes from host | ❌  only if **host** is used|
| `metabolome.csv` | List of KEGG compounds | ❌  only if **metabolome** is used |

## Example

The easiest way to run DMnet is through the included `run_workflow.py` wrapper script.

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--directories` | ✅ | One or more **absolute paths** to sample directories, space-separated and quoted |
| `--dietmicrobe` | ❌ | Enable FooDB-based analysis |
| `--dietmicrobehost` | ❌ | Enable host-based analysis |
| `--metabolome` | ❌ | Include if comparison to known metabolome is wanted |
| `--e-weights` | ❌ | Weight edges by read abundance |
| `--n-weights` | ❌ | Weight nodes by food frequency |
| `--include-orgs` | ❌ | Include organism-level information |
| `--abundance-col` | ❌ | Column name for abundance values (default: `Abundance_RPKs`) |
| `--all-food` | ❌ | Use all foods from FooDB instead of sample-specific diet file |
| `--cores` | ❌ | Number of cores (default: `1`) |
| `--profile` | ❌ | Snakemake profile to use |
| `--dry-run` / `-n` | ❌ | Preview jobs without executing |

!!! tip
    `--cores`, `--profile`, and `--dry-run` are Snakemake-specific arguments.
    See [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for details.

### Run with test data

Test data is included in `Data/test_sample/`. To run the full pipeline on it:

```bash
python run_workflow.py \
    --directories "/absolute/path/to/Data/test_sample" \
    --dietmicrobe \
    --dietmicrobehost \
    --metabolome\
    --e-weights \
    --n-weights \
    --include-orgs \
    --abundance-col "Abundance_RPKs"
```

!!! note
    Always use absolute paths for `--directories`. We recommend doing a dry-run 
    first with `--dry-run` to verify the pipeline is configured correctly before executing.

## Outputs

All outputs are written into your sample directory under `output_fdb/` and/or `output_gen/`,
depending on which analysis modes were enabled.

### Directory Structure

```
my_directory/
├── ko_taxonomy_abundance.csv
├── noquote_ko.txt
├── foodb_foods_dataframe.csv
├── kegg_organisms_dataframe.csv
├── run_info.txt                         # pipeline metadata for this run
├── output_dm/                           # dietmicrobe-based analysis outputs
│   ├── food_meta.csv
│   ├── food_compound_report.html
│   ├── microbe_compound_report.html     # only if --include-orgs and --n-weights
│   ├── MetabolomeComparison_Report.html # only if --metabolome
│   ├── AMON_output/
│   │   ├── AMON_log.txt
│   │   ├── gene_set_1_enrichment.tsv
│   │   ├── kegg_mapper.tsv
│   │   ├── origin_table.tsv
│   │   ├── enrichment_heatmap.png
│   │   ├── co_dict.json
│   │   ├── ko_dict.json
│   │   └── rn_dict.json
│   └── graph/
│       ├── M_nodes_df.csv
│       ├── M_edges_df.csv
│       ├── M_AbundanceDistribution.png
│       ├── M_FoodFrequencyDistribution.png
│       ├── network_summary.txt
│       ├── graph_results.csv
│       └── graph_results_report.html
└── output_dmh/                          # Host-based analysis outputs
    ├── food_meta.csv
    ├── food_compound_report.html
    ├── microbe_compound_report.html     # only if --include-orgs and --n-weights
    ├── MetabolomeComparison_Report.html # only if --metabolome
    ├── AMON_output/
    │   ├── AMON_log.txt
    │   ├── gene_set_1_enrichment.tsv
    │   ├── kegg_mapper.tsv
    │   ├── origin_table.tsv
    │   ├── enrichment_heatmap.png
    │   ├── co_dict.json
    │   ├── ko_dict.json
    │   └── rn_dict.json
    └── graph/
        ├── nodes_df.csv
        ├── edges_df.csv
        ├── AbundanceDistribution.png
        ├── FoodFrequencyDistribution.png
        ├── network_summary.txt
        ├── graph_results.csv
        └── graph_results_report.html
```

### Key Output Files

| File | Description |
|------|-------------|
| `run_info.txt` | Pipeline version, run date, and full config used for this run |
| `food_compound_report.html` | Compounds identified in each food item |
| `microbe_compound_report.html` | Compounds predicted to be produced by microbes |
| `graph_results_report.html` | Network pattern analysis report with Neo4j query results |
| `network_summary.txt` | Summary statistics of the constructed network |
| `graph_results.csv` | Raw graph analysis results |
| `*_nodes_df.csv` | Node dataframe with optional frequency weights |
| `*_edges_df.csv` | Edge dataframe with optional abundance weights |
| `*_AbundanceDistribution.png` | Histogram of edge weights (requires `--e-weights`) |
| `*_FoodFrequencyDistribution.png` | Histogram of node weights (requires `--n-weights`) |

!!! note
    `microbe_compound_report.html` is only generated when both `--include-orgs` and
    `--n-weights` are specified.

## Next Steps 

Once networks and patterns have been found for each sample you can continue to do:

1. [inter-sample comparison](intersample_comp.md)
2. Run your own analyses 😃!