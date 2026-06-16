import streamlit as st
import requests
import os
import signal
import pandas as pd

# Streamlit app
st.title("Whole Genome and/or Metabolome Representation of Diet")

st.write("""
This app allows you to customize a set of food items using FooDB (food metabolomes) that can represent a patients diet. 
Dataframes will be created based on user selections which can be downloaded as a CSV file. 
         
Disclaimer: Version 1.0 of FooDB is used from the CSV file added April7, 2020.
""")

# start the FooDB selection
st.header("FooDB Food Item Selection")

# Load the food CSV (adjust the path if needed) ####################Local Directory 
#food_df = pd.read_csv('/Users/burkhang/Code_Projs/DietMicrobeNet/Data/food.csv')

# Get the current directory of the script
script_dir = os.path.dirname(__file__)
file_path = os.path.join(script_dir, '..', 'data', 'Food.csv')

# Load CSV
food_df = pd.read_csv(file_path)

# Subset to relevant columns
food_df = food_df[['id', 'name', 'name_scientific']]

# Allow the user to select multiple foods
selected_foods = st.multiselect("Select Food Item", food_df['name'].tolist())

# Filter the DataFrame to include only selected foods
foodb_df = food_df[food_df["name"].isin(selected_foods)]

# Display the DataFrame with user input column
if not foodb_df.empty:
    st.write("### Enter values (between 1 and 100) for each selected food")

    food_values = []
    for index, row in foodb_df.iterrows():
        value = st.number_input(
            label=f"Enter value for {row['name']} ({row['id']})",
            min_value=1,
            max_value=100,
            step=1,
            key=f"food_input_{index}"
        )
        food_values.append(value)

    foodb_df["food_frequency"] = food_values

    st.write("### Selected Food Item DataFrame with Input", foodb_df)

    csv = foodb_df.to_csv(index=False)

    st.download_button(
        label="Download Selected Food Items DataFrame as CSV",
        data=csv,
        file_name="foodb_foods_dataframe.csv",
        mime="text/csv",
        key="download_foodb"
    )
else:
    st.write("No foods selected.")

# put a kill button so user can close application easily 
st.write("""
**Once you're done hit the *Stop Server* button wait 5 seconds close the tab**
""")
if st.button("Stop Server"):
    st.warning("Stopping the Streamlit server...")
    os.kill(os.getpid(), signal.SIGTERM)