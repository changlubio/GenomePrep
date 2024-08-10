import streamlit as st
import pandas as pd

st.set_page_config(
    page_title="open snps",
    page_icon="👋",
    layout="wide"
)

st.write("# OpenSNP uploads")

# Update data

# Load data
uploads = pd.read_csv('data/opensnp_uploads.csv')
# to datetime by specified format
uploads['time'] = pd.to_datetime(uploads['time'], format='%d.%m.%Y %H:%M')
# 10.08.2024 20:35

# 2 columns
col1, col2 = st.columns([3, 2])

with col1:
    # make data frame by year-month, count the number of uploads, sort by index, and count the total number of uploads by year-month
    # then plot the data
    time_points = uploads['time'].dt.to_period('M').value_counts().sort_index()
    time_points.index = time_points.index.to_timestamp()
    # set year-month
    time_points.index = time_points.index.strftime('%Y-%m')

    # show plot, x-axis is year-month
    st.line_chart(
        time_points,
        x_label='Year-Month',
        y_label='Number of uploads',
    )

with col2:
    # show cumulative sum of uploads
    cumsum = time_points.cumsum()
    st.line_chart(
        cumsum,
        x_label='Year-Month',
        y_label='Cumulative sum of uploads',
    )