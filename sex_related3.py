import marimo

__generated_with = "0.14.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pyarrow as pa

    return mo, pd, plt, sns


@app.cell
def _():
    import hail as hl

    hl.init()
    return (hl,)


@app.cell
def _(hl):
    mt = hl.read_matrix_table('/home/ubuntu/kawasaki/kawasaki_filtered.mt')
    return (mt,)


@app.cell
def _(mt):
    print("--- MatrixTable Schema ---")
    mt.describe()

    return


@app.cell
def _(hl, mt):

    # --- Sex Estimation Function ---
    def estimate_sex_from_x_chromosome_alternative_plot(
        mt: hl.MatrixTable,
        female_threshold: float = -0.25,  # If F-stat is below this threshold, call female
        male_threshold: float = -0.25     # If F-stat is above this threshold, call male
    ):

        print("Starting sex estimation...")
        # 1. Call hl.impute_sex function.
        # This function operates on mt.GT and returns a Hail Table.
        sex_table = hl.impute_sex(
            mt.GT,
            female_threshold=female_threshold,
            male_threshold=male_threshold
        )
        print("Sex estimation completed.")

        # Convert to Pandas DataFrame
        sex_df = sex_table.to_pandas()

        # First convert True/False to text and then add unknowns
        sex_df['is_female_str'] = sex_df['is_female'].map({True: 'Women', False: 'Men'})
        sex_df['is_female_str'].fillna('Unknown', inplace=True)  # Fill NaNs

        # Convert 'is_female_str' column to categorical type
        # This helps Seaborn handle categorical data better
        sex_df['is_female_str'] = sex_df['is_female_str'].astype('category')
        # >>> END <<<

        print("\n--- First 5 Rows of Sex Estimation Results ---")
        print(sex_df.head())

        print("\n--- Predicted Sex Distribution ---")
        # Use 'is_female_str' column here
        print(sex_df['is_female_str'].value_counts())
    
        return sex_df

    # 2. Call estimate_sex_from_x_chromosome_alternative_plot function.
    # The function will both visualize (if implemented) and return the Hail Table results.
    sex_df = estimate_sex_from_x_chromosome_alternative_plot(mt)

    # 3. Display the sex estimation table (Hail Table).
    # This shows the schema and first few rows of the table.
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

    return


@app.cell
def _(mo, pd):
    # CSV read
    meta_data_df = pd.read_csv("/home/ubuntu/kawasaki/meta.csv", sep='\t')

    # Marimo'da göster
    mo.ui.table(meta_data_df)
    return (meta_data_df,)


@app.cell
def _():
    # --- Set Column Names ---
    # These column names are set according to the header row of your meta.csv file.
    id_column_name = 'Case Number'
    gender_column_name = 'Sex'

    return gender_column_name, id_column_name


@app.cell
def _(id_column_name, meta_data_df):
    # Format the Case Numbers in the metadata with the 'M_' prefix
    # Example: 15 -> M_15, 256 -> M_256
    meta_data_df['formatted_sample_id'] = 'M_' + meta_data_df[id_column_name].astype(str)
    print("\n--- Metadata IDs formatted with 'M_' prefix (First 5 rows) ---")
    print(meta_data_df[['Case Number', 'formatted_sample_id']].head())

    return


@app.cell
def _(meta_data_df):
    # --- Filter sex_df (Hail predictions) to include only patients present in the metadata ---
    # Get all FORMATTED Case Numbers (patient IDs) from meta_data_df
    patient_formatted_ids = meta_data_df['formatted_sample_id'].tolist()

    return (patient_formatted_ids,)


@app.cell
def _(patient_formatted_ids, sex_df):
    # Select rows in sex_df where the 's' column (sample IDs) matches these formatted IDs
    sex_df_patient = sex_df[sex_df['s'].astype(str).isin(patient_formatted_ids)].copy()

    print(f"\nsex_df (Only Patients in Metadata) size: {len(sex_df_patient)}")
    print("sex_df (Patients) first 5 rows:")
    print(sex_df_patient.head())

    return (sex_df_patient,)


@app.cell
def _(gender_column_name, meta_data_df, pd, sex_df_patient):
    # --- Merge Filtered Data and Compare ---
    # We will now merge using the 'formatted_sample_id' column.
    comparison_df = pd.merge(
        sex_df_patient,
        meta_data_df[['formatted_sample_id', gender_column_name]],
        left_on='s',
        right_on='formatted_sample_id',
        how='left'
    )


    return (comparison_df,)


@app.cell
def _(comparison_df, gender_column_name):
    # The 'Sex' column in the metadata is already in 'Male'/'Female' format, so we can use it directly.
    comparison_df.rename(columns={gender_column_name: 'true_sex_str'}, inplace=True)

    return


@app.cell
def _(comparison_df):
    # Check the accuracy of the predictions
    comparison_df['prediction_correct'] = (comparison_df['is_female_str'] == comparison_df['true_sex_str'])

    print("\n--- Accuracy of Sex Predictions (Only Kawasaki Patients) ---")
    print(comparison_df['prediction_correct'].value_counts())

    return


@app.cell
def _(comparison_df):
    # By the way, a note: during quality control, the sample with ID M_36 had a sex that did not match the one in the phenotype file. You should check this during QC.

    # View incorrect predictions (if any)
    incorrect_predictions = comparison_df[comparison_df['prediction_correct'] == False]
    if not incorrect_predictions.empty:
        print("\n--- Kawasaki Patient Samples with Incorrect Sex Predictions (Details) ---")
        print(incorrect_predictions[['s', 'is_female_str', 'true_sex_str', 'f_stat']])
    else:
        print("\nSex prediction is correct for all Kawasaki patient samples!")

    return


@app.cell
def _(comparison_df, mo):
    # Overall accuracy percentage
    accuracy_percentage = comparison_df['prediction_correct'].mean() * 100
    print(f"\nOverall Prediction Accuracy for Kawasaki Patient Samples: {accuracy_percentage:.2f}%")

    # Display the resulting DataFrame in Marimo
    mo.ui.table(comparison_df)

    return


if __name__ == "__main__":
    app.run()
