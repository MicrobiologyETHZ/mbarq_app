"""

Adapted from https://towardsdatascience.com/creating-multipage-applications-using-streamlit-efficiently-b58a58134030

This file is the framework for generating multiple Streamlit applications
through an object oriented framework.

"""

import streamlit as st


class MultiPage:
    """Framework for combining multiple streamlit applications."""

    def __init__(self) -> None:
        """Constructor class to generate a list which will store all our applications as an instance variable."""
        self.pages = []

    def add_page(self, title, func) -> None:
        """Class Method to Add pages to the project
        Args:
            title ([str]): The title of page which we are adding to the list of apps

            func: Python function to render this page in Streamlit
        """

        self.pages.append({

            "title": title,
            "function": func
        })

    def run(self):
        # Radio button to select the page to run
        page = st.sidebar.radio(
            'App Navigation',
            self.pages,
            format_func=lambda page: page['title']
        )
        # Run the app function
        page['function']()