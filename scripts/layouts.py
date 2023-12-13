import streamlit as st
from scripts.datasets import define_color_scheme

ALPHABET_COLORS, APP_COLORS, ALL_COLORS = define_color_scheme()


def pca_layout(cds):
    with st.expander('Show PCA'):

        st.write('### PCA Options')
        c1, c2, c3, c4 = st.columns(4)
        num_components = c1.number_input("Number of PCs", min_value=2, max_value=50, value=10)
        num_genes = c2.number_input("Number of genes to use", min_value=int(num_components),
                                    value=int(max(100, cds.count_data.shape[0] * 0.1)),
                                    max_value=int(cds.count_data.shape[0]),
                                    help='By default, uses top 10% most variable barcodes')
        choose_by = 'variance'
        num_genes = int(num_genes)
        num_components = int(num_components)
        pc_df, percent_variance = cds.get_principal_components(num_components, num_genes, choose_by)
        missing_meta = " ,".join(list(pc_df[pc_df.isna().any(axis=1)].index))
        if missing_meta:
            st.write(f"The following samples have missing_metadata and will not be shown: {missing_meta}")
        pc_df = pc_df[~pc_df.isna().any(axis=1)]  # todo this should be included in the function
        pc_labels = [f'PC{i}' for i in range(1, num_components + 1)]
        experiment_vars = [c for c in pc_df.columns if c not in pc_labels]
        desired_width, desired_height, font_size = None, None, 24
        if st.checkbox('Modify PCA plot dimensions'):
            wc, hc, font = st.columns(3)
            desired_width = wc.number_input('Width', min_value=200, max_value=1000, value=800)
            desired_height = hc.number_input('Height', min_value=200, max_value=1000, value=400)
            font_size = font.number_input("Font size", min_value=8, max_value=40, value=24)
        pc_x = c1.selectbox('X-axis component', pc_labels)
        pc_y = c2.selectbox('Y-axis component', [pc for pc in pc_labels if pc != pc_x])
        pc_to_highlight = c3.radio('Variable to highlight', experiment_vars)
        pc_to_symbol = c4.radio('Variable to show as symbol', [None] + experiment_vars)
        pc_df = pc_df.sort_values(pc_to_highlight)
        fig1, fig2, fig3 = cds.pca_figure(pc_df, pc_x, pc_y, pc_to_highlight, percent_variance, pc_to_symbol,
                                          experiment_vars,
                                          ALL_COLORS, desired_width, desired_height, font_size=font_size)
        st.plotly_chart(fig1, use_container_width=False)
        c5, c6 = st.columns(2)
        c5.write('### Scree Plot')
        c5.plotly_chart(fig2)
        c6.write(f'### PCs summarized by {pc_to_highlight}')
        c6.plotly_chart(fig3, use_container_width=True)


def barcode_abundance_layout(cds):
    st.write('## Barcode Abundance')
    with st.expander('Show Barcode Abundance'):
        c1, c2 = st.columns(2)
        compare_condition = c1.selectbox('Which conditions to compare?', sorted(cds.sample_data.columns))
        condition_categories = c1.multiselect(f'Categories of {compare_condition} to display',
                                              ['All'] + list(cds.sample_data[compare_condition].unique()),
                                              default='All')
        filter_condition = c2.selectbox("Filter by", ['No filter'] + list(cds.sample_data.columns),
                                        index=0)
        if filter_condition == 'No filter':
            filter_categories = []
        else:
            filter_categories = c2.multiselect(f'Which category(ies) of {filter_condition} to keep?',
                                               list(cds.sample_data[filter_condition].unique()))
        if 'All' in condition_categories:
            condition_categories = list(cds.sample_data[compare_condition].unique())
        default_genes =  [ex for ex in  ['dcuS', 'dcuR'] if ex in cds.count_data[cds.gene_name_col].unique()]
        genes = st.multiselect("Choose gene(s) of interest", cds.count_data[cds.gene_name_col].unique(), default=default_genes)
        if len(genes) * len(condition_categories) > 40:
            st.write('Too many genes/categories to display, consider choosing fewer genes.')
        else:
            gene_df = cds.count_data[cds.count_data[cds.gene_name_col].isin(genes)]
            ab_sample_df = cds.sample_data[cds.sample_data[compare_condition].isin(condition_categories)]
            if filter_categories:
                ab_sample_df = ab_sample_df[ab_sample_df[filter_condition].isin(filter_categories)]
            gene_df = (gene_df.melt(id_vars=[cds.barcode_col, cds.gene_name_col],
                                    value_name='log2CPM', var_name=cds.sample_id_col)
                       .merge(ab_sample_df, how='inner', on=cds.sample_id_col)
                       .sort_values(compare_condition))

            col1, col2 = st.columns(2)
            groupBy = col1.radio('Group by', [cds.gene_name_col, compare_condition])
            colorBy = [c for c in [cds.gene_name_col, compare_condition] if c != groupBy][0]

            # Choose Plot Type
            box = "Box"
            violin = "Violin"
            plotType = col2.radio('Plot Type', (box, violin))
            if plotType == box:
                fig = cds.barcode_abundance_plot(gene_df, groupBy, colorBy, ALL_COLORS)
            if plotType == violin:
                fig = cds.barcode_abundance_plot(gene_df, groupBy, colorBy, ALL_COLORS, box=False)
            st.plotly_chart(fig, use_container_width=True)
