# Output Files

All outputs are written into your sample directory under `output_fdb/` and/or `output_gen/`,
depending on which analysis modes were enabled.

## Directory Structure

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

## Key Output Files

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