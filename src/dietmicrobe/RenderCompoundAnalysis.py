import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Generate R Markdown report with parameters.")
parser.add_argument("--food_file", default="default.csv", help="Path to the data file.")
parser.add_argument("--output", default="filepath.html", help="Path to HTML report file")
args = parser.parse_args()

this_dir = os.path.dirname(__file__)
rmd_path = os.path.join(this_dir, "CompoundAnalysis_FooDB.Rmd")
# Ensure output path is absolute
output_path = os.path.abspath(args.output)

# Construct the R command with parameters
r_command = f'''
Rscript -e "rmarkdown::render(
  input = '{rmd_path}',
  params = list(food_file = '{args.food_file}'),
  output_file = '{output_path}')"
'''

# Execute the R command
subprocess.run(r_command, shell=True)
