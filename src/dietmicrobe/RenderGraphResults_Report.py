import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Generate HTML report of graph query results.")
parser.add_argument("--patterns", default="default.csv", help="Path to the node data file.")
parser.add_argument("--rxn_json", default="rxn_default.json", help="Path to JSON file created from AMON output.")
parser.add_argument("--output", default="filepath.html", help="Path to HTML report file")
args = parser.parse_args()

this_dir = os.path.dirname(__file__)
rmd_path = os.path.join(this_dir, "GraphResults_Report.Rmd")
# Ensure output path is absolute
output_path = os.path.abspath(args.output)

# Construct the R command with parameters
r_command = f'''
Rscript -e "rmarkdown::render(
  input = '{rmd_path}',
  params = list(patterns = '{args.patterns}', reactions =  '{args.rxn_json}'),
  output_file = '{output_path}')"
'''

# Execute the R command
subprocess.run(r_command, shell=True)