import streamlit as st
import pandas as pd
import os

st.set_page_config(
    page_title="openSNP",
    page_icon="👋",
    layout="wide"
)

st.write("# OpenSNP uploads")

# Update data
from opensnp_update import make_updates
# check if file exists
if not os.path.exists('opensnp_uploads.csv'):
    make_updates()

# Load data
uploads = pd.read_csv('opensnp_uploads.csv')
# to datetime by specified format
uploads['time'] = pd.to_datetime(uploads['time'], format='%d.%m.%Y %H:%M')
# 10.08.2024 20:35

# print total number of entries, distinct users, and number entries in last 30 days
st.write(f"Total number of entries: {len(uploads)}. Distinct users: {len(uploads['userID'].unique())}")

# 2 columns
col1, _, col2 = st.columns([2.5,0.5, 2])

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