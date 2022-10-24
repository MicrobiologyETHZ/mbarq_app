import numpy as np
import pandas as pd
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode


def app():
### Test dataframe
    df = pd.DataFrame(np.random.randint(0, 100, size=(100, 4)), columns=list('ABCD'))

    @st.cache
    def load_gridOptions(df):
        """Cached function to load grid options"""
        gb = GridOptionsBuilder.from_dataframe(df)
        gb.configure_pagination(paginationAutoPageSize=True)  # Add pagination
        gb.configure_side_bar()  # Add a sidebar
        gb.configure_selection('multiple',
                               use_checkbox=True)  # Enable multiselection, if nested rows add this -->  groupSelectsChildren="Group checkbox select children"
        gridOptions = gb.build()
        return gridOptions

    st.write('Custom aggrid')
    grid_response = AgGrid(
        df,
        gridOptions=load_gridOptions(df),
        data_return_mode='AS_INPUT',
        # update_mode='NO_UPDATE',
        update_mode='MODEL_CHANGED',
        # header_checkbox_selection_filtered_only=True,
        fit_columns_on_grid_load=False,
        # theme='blue', #Add theme color to the table
        enable_enterprise_modules=True,
        height=350,
        width='100%',
        key='an_unique_key',
        reload_data=False
    )

    selected = grid_response[
        'selected_rows']  ### This creates a new dataframe with the 'selected_rows' in the grid by the user
    df2 = pd.DataFrame(selected)  # Pass the selected rows to a new dataframe df
    st.write('Selected elements')
    AgGrid(df2)