import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pyarrow as pa
    import hail as hl

    hl.init()
    import numpy as np
    return hl, np, pd, plt, sns


@app.cell
def _(hl):
    mt = hl.read_matrix_table('kawasaki_filtered.mt')
    return (mt,)


@app.cell
def _(mt):
    mt.describe()
    return


@app.cell
def _(hl, mt):
    # --- Sex Estimation Function ---
    def estimate_sex_from_x_chromosome_alternative_plot(
        mt: hl.MatrixTable,
        female_threshold: float = 0.4,  
        male_threshold: float = 0.4     
    ):

        print("Starting sex estimation...")
        sex_table = hl.impute_sex(
            mt.GT,
            female_threshold=female_threshold,
            male_threshold=male_threshold
        )
        print("Sex estimation completed.")

        sex_df = sex_table.to_pandas()

        sex_df['is_female_str'] = sex_df['is_female'].map({True: 'Female', False: 'Male'})
        sex_df['is_female_str'].fillna('Unknown', inplace=True)  # Fill NaNs

        sex_df['is_female_str'] = sex_df['is_female_str'].astype('category')

        print("\n--- First 5 Rows of Sex Estimation Results ---")
        print(sex_df.head())

        print("\n--- Predicted Sex Distribution ---")
        print(sex_df['is_female_str'].value_counts())

        return sex_df


    sex_df = estimate_sex_from_x_chromosome_alternative_plot(mt)

    print("\n--- Sex Estimation Results (Hail Table Summary) ---")
    return (sex_df,)


@app.cell
def _(plt, sex_df, sns):
    plt.figure(figsize=(10, 6))

    sns.histplot(
        sex_df,
        x='f_stat',
        hue='is_female_str', 
        kde=True,
        bins=50,
        palette='viridis',
        alpha=0.7
    )

    plt.savefig('sex_estimation.svg', format='svg', bbox_inches='tight')
    plt.show()
    return


@app.cell
def _(mt):
    mt.rows().show(200)
    return


@app.cell
def _(mt):
    mt.describe()
    return


@app.cell
def _(mt):
    mt.pheno.Sex.show()
    return


@app.cell
def _(mt, pd, sex_df):

    cols_ht = mt.cols()

    reported_sex_ht = cols_ht.select(
        reported_sex = cols_ht.pheno.Sex,
        status = cols_ht.case_control_status
    )

    reported_sex_df = reported_sex_ht.to_pandas().reset_index()

    print("--- Reported Gender Information ---")
    print(reported_sex_df.head())
    print(f"Total {len(reported_sex_df)} Reported gender information was obtained for the sample.")

    cols_ht = mt.cols()


    if 's_corrected' not in sex_df.columns:
        sex_df_with_key = sex_df.reset_index()
    else:
        sex_df_with_key = sex_df.copy()


    merged_df = pd.merge(
        sex_df_with_key,  
        reported_sex_df,  
        on='s_corrected', 
        how='left'        
    )

    print("Merge completed")
    print(merged_df.head())
    return (merged_df,)


@app.cell
def _(merged_df, np):
    print("\n--- Only 'Case' Examples Are Filtered ---")

    cases_only_df = merged_df[merged_df['status'].str.lower() == 'case'].copy()

    if cases_only_df.empty:
        print("WARNING: No instance found with status 'case'. Please check the value in the 'status' column.")
    else:
        print(f"Total {len(cases_only_df)} 'case' samples were found.")
        print("First 5 rows of the filtered 'case' data:")
        print(cases_only_df.head())

    print("\n---'Case' Group Data Cleaning and Comparison ---")

    sex_map = {
        'Male': 'Male',
        'Female': 'Female'}

    cases_only_df['reported_sex_clean'] = cases_only_df['reported_sex'].map(sex_map).fillna('Unknown')

    conditions = [
        # 1. If calculation is 'Unknown' (f_stat fell between threshold values)
        (cases_only_df['is_female_str'] == 'Unknown'),

        # 2. If reported is 'Unknown' (Should not occur since we filtered for 'case')
        (cases_only_df['reported_sex_clean'] == 'Unknown'),

        # 3. If Calculated and Reported MATCH
        (cases_only_df['is_female_str'] == cases_only_df['reported_sex_clean']),

        # 4. If Calculated and Reported DO NOT MATCH
        (cases_only_df['is_female_str'] != cases_only_df['reported_sex_clean'])
    ]

    choices = [
        'Inferred as Unknown', 
        'Reported as Unknown',  
        'Match',                
        'Mismatch'              
    ]

    cases_only_df['validation_status'] = np.select(conditions, choices, default='Error')

    print("\n--- 'Case' Group Comparison Status Summary ---")
    print(cases_only_df['validation_status'].value_counts())
    return (cases_only_df,)


@app.cell
def _(cases_only_df):
    print("\n--- 'Case' Group SEX MISMATCHES ---")


    mismatches_cases_df = cases_only_df[cases_only_df['validation_status'] == 'Mismatch']

    if mismatches_cases_df.empty:
        print("No mismatch found between reported and inferred sex in the 'Case' group.")
    else:
        print(f"A total of {len(mismatches_cases_df)} mismatches were found in 'case' samples:")
        print(mismatches_cases_df[[
            's_corrected', 
            'f_stat', 
            'is_female_str',    
            'reported_sex',      
            'status',            
            'validation_status'
        ]])
    return


@app.cell
def _(cases_only_df, plt, sns):
    print("\n--- 'Case' Group Comparison Plot ---")

    plt.figure(figsize=(10, 8)) 
    sns.stripplot(
        data=cases_only_df,
        x='reported_sex_clean',  
        y='f_stat',              
        hue='is_female_str',     
        jitter=True,             
        dodge=True,              
        alpha=0.7,               
        palette='Set1'           
    )

    plt.title("Sex Validation ('Case' Group Only) vs. F-Stat")
    plt.xlabel("Reported Sex")
    plt.ylabel('F-Stat (Calculated)')
    plt.legend(title='Calculated Sex')
    plt.grid(axis='y', linestyle='--', alpha=0.7) 

    plt.savefig('sex_validation_case_group.svg', format='svg', bbox_inches='tight')
    plt.show()

    return


if __name__ == "__main__":
    app.run()
