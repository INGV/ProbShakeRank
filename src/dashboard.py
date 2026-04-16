import streamlit as st
import os
import pandas as pd
from PIL import Image, ImageChops


# -----------------------------
# SETUP
# -----------------------------
out_dir = "./OUTPUT"
events = [ev for ev in os.listdir(out_dir) if ev != ".DS_Store"]

st.set_page_config(layout="wide")
col1, col2 = st.columns([2, 3], gap="large")


# -----------------------------
# LEFT PANEL
# -----------------------------
with col1:

    selected_event = st.selectbox("Event", events)
    out_dir_ev = os.path.join(out_dir, selected_event)

    # Metadata
    metadata_path = os.path.join(out_dir_ev, "metadata.txt")

    with open(metadata_path, "r") as f:
        lines = f.readlines()

    event_info = {
        line.split(":")[0].strip(): line.split(":", 1)[1].strip()
        for line in lines if ":" in line
    }

    st.markdown(
        f"**Mw**: {event_info.get('Mw', 'N/A')} | "
        f"**Time**: {event_info.get('Time', 'N/A')} | "
        f"**LAT**: {event_info.get('LatEvent', 'N/A')} "
        f"**LON**: {event_info.get('LonEvent', 'N/A')}"
    )

    # GMPE / IM selection
    out_folders = [
        f.replace("OUTPUT_", "")
        for f in os.listdir(out_dir_ev)
        if f.startswith("OUTPUT_")
    ]

    gmpes = sorted(set(f.split("_")[0] for f in out_folders))
    imts = sorted(set(f.split("_")[1] for f in out_folders))

    selected_gmpe = st.selectbox("GMPE", gmpes)
    selected_imt = st.selectbox("IM", imts)

    # Statistic selection
    stat_dir = os.path.join(out_dir_ev, f"OUTPUT_{selected_gmpe}_{selected_imt}", "STATISTICS")

    stat_files = [
        f for f in os.listdir(stat_dir)
        if f.startswith("Stat_") and f.endswith(".png")
    ]

    stat_names = [f.replace("Stat_", "").replace(".png", "") for f in stat_files]

    selected_stat = st.selectbox("Statistic", stat_names)

    img = os.path.join(stat_dir, f"Stat_{selected_stat}.png")

    st.image(img)


# -----------------------------
# RIGHT PANEL
# -----------------------------
with col2:

    st.subheader("GMPEs Ranking")

    selected_metric = st.selectbox(
        "Metric for ranking",
        ("LLH_Score", "Gambling_Score")
    )

    multi_rank_file = os.path.join(out_dir_ev, "RANK/MultiIMs_Ranking.txt")

    fig_dir = os.path.join(out_dir_ev, f"OUTPUT_{selected_gmpe}_{selected_imt}", "RANK_FIGURES")

    img1 = os.path.join(fig_dir, "Normalized_Residuals.png")


    # -----------------------------
    # LLH SCORE MODE
    # -----------------------------
    if selected_metric == "LLH_Score":

        score_file = os.path.join(out_dir_ev, f"RANK/LLH_Score_{selected_imt}.txt")

        df = pd.read_csv(
            score_file,
            delim_whitespace=True,
            skiprows=1,
            header=None,
            names=["GMPE", "LLH_Score"]
        )

        df = df.sort_values("LLH_Score")
        df["LLH_Score"] = df["LLH_Score"].map(lambda x: f"{x:.3f}")
        df_display = df.reset_index(drop=True)
        df_display.index = df_display.index + 1

        # Multivariate
        df_multi = pd.read_csv(
            multi_rank_file,
            delim_whitespace=True,
            skiprows=1,
            header=None,
            names=["GMPE", "Multivariate_LLH", "AIC", "BIC"]
        )

        df_multi = df_multi.sort_values("Multivariate_LLH")
        df_multi["Multivariate_LLH"] = df_multi["Multivariate_LLH"].map(lambda x: f"{x:.3f}")
        df_multi["AIC"] = df_multi["AIC"].map(lambda x: f"{x:.1f}")
        df_multi["BIC"] = df_multi["BIC"].map(lambda x: f"{x:.1f}")

        df_multi_display = df_multi.reset_index(drop=True)
        df_multi_display.index = df_multi_display.index + 1

        # Tables
        t1, t2 = st.columns(2)

        with t1:
            st.markdown(f"{selected_imt} LLH Scores")
            st.dataframe(df_display, use_container_width=True)

        with t2:
            st.markdown("Multivariate Rankings")
            st.dataframe(df_multi_display, use_container_width=True)

        # Images (side-by-side, controlled size)
        img2 = os.path.join(fig_dir, "POI_LLH.png")

        c1, c2 = st.columns(2)

        with c1:
            st.image(img1)

        with c2:
            st.image(img2)


    # -----------------------------
    # GAMBLING MODE
    # -----------------------------
    else:

        score_file = os.path.join(out_dir_ev, f"RANK/Gambling_Score_{selected_imt}.txt")

        df = pd.read_csv(
            score_file,
            delim_whitespace=True,
            skiprows=1,
            header=None,
            names=["GMPE", "Gambling_Score"]
        )

        df = df.sort_values("Gambling_Score", ascending=False)
        df["Gambling_Score"] = df["Gambling_Score"].map(lambda x: f"{x:.3f}")

        df_display = df.reset_index(drop=True)
        df_display.index = df_display.index + 1

        st.dataframe(df_display, use_container_width=True)

        # Images
        img2 = os.path.join(fig_dir, "POI_Gambling.png")

        c1, c2 = st.columns(2)

        with c1:
            st.image(img1)

        with c2:
            st.image(img2)