# DietMicrobe[Host]Net: Network Modeling of Dietary Effect on Microbial Metabolism

The purpose of this code will be to create a metabolic network where nodes represent compounds and edges represent reactions. Compounds will originate from either food items (e.g. apple), microbes, host, a combination thereof or neither (cannot be predicted). Graphs are created in memory, where patterns of origin will be found and stored (for details regarding patterns see [docs](https://lozuponelab.github.io/DietMicrobeNet/latest/usage/#step-5-build-graph-and-extract-patterns) ). Information such as food frequency, read abundance, and taxonomy can also be conserved within the graph depending on user preferences. This will allow a user to find potential instances of dietary metabolism by microbes and host.

## General Workflow 

Every user will need to set up and install this program the same way. To do this go to [Installation](https://lozuponelab.github.io/DietMicrobeNet/latest/installation/) in the repo's [Docs](https://lozuponelab.github.io/DietMicrobeNet/latest/). 

After getting set up, you will have decide if you want to run each step manually or use the provided Snakemake Workflow. 

* For **manual approach** follow instructions in the usage section: [Running the pipeline](https://lozuponelab.github.io/DietMicrobeNet/latest/usage/).
* For **Snakemake Workflow** (RECOMMENDED) follow [Quick Start Example](https://lozuponelab.github.io/DietMicrobeNet/latest/quickstart/) instructions
* [Inter-Sample Comparison](https://lozuponelab.github.io/DietMicrobeNet/latest/intersample_comp/) can be performed after manual or snakemake approach

## Outputs

The main outputs of this program will be: 

1. HTML report on what we know about the compounds originating from the **food** (IF using specific FFQs, if all foods are used no report is created) 
2. HTML report on what we know about the compounds originating from the **microbes** (IF taxonomic and abundance information is provided)
3. HTML report on compounds and reactions involved in dietary metabolism 
4. Figures describing distribution of node and edge weights (representing food frequency and read abundance respectively)
5. Summary text file for each graph 


If you include a metabolome: 

1. HTML report on which compounds found in the patterns were also found in a metabolomics experiment that you have provided. 

If you run the graph comparison script you will generate: 

1.  Summary files on the results of three different patterns that were searched for in the graph 
2. Figures describing clustered similarity of graphs 

## Contact

For all programmatic questions contact: 

- Angela Sofia Burkhart Colorado (main developer) at angelasofia.burkhartcolorado@cuanschutz.edu

This work was done under the auspices of the [Lozupone Lab](https://www.lozuponelab.com/). If you have any questions regarding future publications and general management of this tool contact:

- Dr. Catherine Lozupone (PI) at catherine.lozupone@cuanschutz.edu