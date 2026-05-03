import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from scipy import stats

st.set_page_config(layout="wide")
st.title("🧬 Gene Expression Cancer RNA-Seq Agent")
st.write("Compare two cancer types using RNA-seq gene expression data.")

st.subheader("Step 1: Load files")

data_file = st.file_uploader("Upload data.csv", type="csv", key="data_file")
label_file = st.file_uploader("Upload labels.csv or label.csv", type="csv", key="label_file")

if data_file is not None and label_file is not None:
    data_df = pd.read_csv(data_file)
    labels_df = pd.read_csv(label_file)

    st.success("Both files uploaded successfully.")
    st.write("**data.csv shape:**", data_df.shape)
    st.write("**labels file shape:**", labels_df.shape)

    st.subheader("Step 2: Preview data")
    st.write("**Gene expression preview**")
    st.dataframe(data_df.head())

    st.write("**Labels preview**")
    st.dataframe(labels_df.head())

    st.subheader("Step 3: Detect label column")

    possible_label_cols = [col for col in labels_df.columns if labels_df[col].dtype == "object"]

    if len(possible_label_cols) == 0:
        st.error("No text label column found in the labels file.")
    else:
        label_col = st.selectbox("Choose the cancer label column", possible_label_cols)

        labels = labels_df[label_col].dropna().astype(str)
        unique_groups = sorted(labels.unique())

        st.write("**Cancer groups found:**", unique_groups)

        if len(unique_groups) < 2:
            st.error("Need at least two cancer groups to compare.")
        else:
            st.subheader("Step 4: Choose two cancer types")

            col1, col2 = st.columns(2)
            with col1:
                group1 = st.selectbox("Group 1", unique_groups, index=0)
            with col2:
                group2 = st.selectbox("Group 2", unique_groups, index=1 if len(unique_groups) > 1 else 0)

            st.subheader("Step 5: Analysis settings")

            max_genes = min(500, data_df.shape[1])
            num_genes = st.slider("How many genes to analyze (faster = fewer)", 50, max_genes, 200, step=50)

            if st.button("🚀 Run RNA-seq Comparison", type="primary"):
                if group1 == group2:
                    st.error("Please choose two different cancer groups.")
                else:
                    if len(data_df) != len(labels_df):
                        st.error("data.csv and labels file do not have the same number of rows.")
                    else:
                        analysis_df = data_df.copy()
                        analysis_df["CancerType"] = labels.values

                        gene_cols = analysis_df.select_dtypes(include=[np.number]).columns.tolist()
                        gene_cols = gene_cols[:num_genes]

                        g1_df = analysis_df[analysis_df["CancerType"] == group1]
                        g2_df = analysis_df[analysis_df["CancerType"] == group2]

                        st.write(f"Samples in {group1}: {len(g1_df)}")
                        st.write(f"Samples in {group2}: {len(g2_df)}")

                        results = []

                        for gene in gene_cols:
                            g1_vals = g1_df[gene].dropna()
                            g2_vals = g2_df[gene].dropna()

                            if len(g1_vals) > 1 and len(g2_vals) > 1:
                                _, pval = stats.ttest_ind(g1_vals, g2_vals, equal_var=False)

                                g1_mean = g1_vals.mean()
                                g2_mean = g2_vals.mean()

                                if g1_mean != 0:
                                    fold_change = g2_mean / g1_mean
                                else:
                                    fold_change = np.nan

                                results.append({
                                    "Gene": gene,
                                    f"{group1}_mean": round(g1_mean, 4),
                                    f"{group2}_mean": round(g2_mean, 4),
                                    "Fold_Change": round(fold_change, 4) if pd.notnull(fold_change) else np.nan,
                                    "P_Value": pval
                                })

                        if len(results) == 0:
                            st.warning("No valid results were produced.")
                        else:
                            results_df = pd.DataFrame(results)
                            results_df["P_Value"] = pd.to_numeric(results_df["P_Value"], errors="coerce")
                            results_df["Fold_Change"] = pd.to_numeric(results_df["Fold_Change"], errors="coerce")

                            results_df = results_df.dropna(subset=["P_Value", "Fold_Change"])
                            results_df = results_df.sort_values("P_Value").reset_index(drop=True)

                            results_df["neg_log10_p"] = -np.log10(results_df["P_Value"].replace(0, 1e-300))

                            st.subheader("Top results")
                            st.dataframe(results_df.head(20))

                            fig = px.scatter(
                                results_df.head(200),
                                x="Fold_Change",
                                y="neg_log10_p",
                                hover_name="Gene",
                                color="P_Value",
                                title=f"Volcano-style Plot: {group1} vs {group2}",
                                labels={"neg_log10_p": "-log10(P_Value)"}
                            )
                            st.plotly_chart(fig, use_container_width=True)

                            sig_hits = len(results_df[results_df["P_Value"] < 0.05])
                            top_gene = results_df.iloc[0]["Gene"]

                            st.markdown(f"""
### 🎯 Research Summary

**Compared:** {group1} vs {group2}  
**Genes analyzed:** {len(gene_cols)}  
**Significant hits (p < 0.05):** {sig_hits}  
**Top gene:** {top_gene}  
**Interpretation:** This app screens RNA-seq gene-expression features to identify genes with different average expression between two cancer types.

**Note:** This is a prototype screening workflow for demo purposes, not a publication-ready differential expression pipeline.
""")

                            csv_data = results_df.to_csv(index=False)
                            st.download_button(
                                "💾 Download results CSV",
                                csv_data,
                                "rna_seq_results.csv",
                                "text/csv"
                            )
else:
    st.info("Upload both data.csv and labels.csv to begin.")
