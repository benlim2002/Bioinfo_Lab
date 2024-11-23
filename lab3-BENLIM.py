import streamlit as st
import numpy as np


#needleman func
def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2): #set default perim first
    n, m = len(seq1) + 1, len(seq2) + 1
    dp, traceback = np.zeros((n, m), dtype=int), np.zeros((n, m), dtype=object)

    for i in range(1, n):
        dp[i][0], traceback[i][0] = dp[i - 1][0] + gap, (i - 1, 0)
    for j in range(1, m):
        dp[0][j], traceback[0][j] = dp[0][j - 1] + gap, (0, j - 1)

    for i in range(1, n):
        for j in range(1, m):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal, up, left = dp[i - 1][j - 1] + score, dp[i - 1][j] + gap, dp[i][j - 1] + gap
            dp[i][j] = max(diagonal, up, left)
            traceback[i][j] = (i - 1, j - 1) if dp[i][j] == diagonal else ((i - 1, j) if dp[i][j] == up else (i, j - 1))

    return dp, traceback


#waterman func
def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1) + 1, len(seq2) + 1
    dp, traceback = np.zeros((n, m), dtype=int), np.zeros((n, m), dtype=object)
    max_score, max_pos = 0, (0, 0)

    for i in range(1, n):
        for j in range(1, m):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal, up, left = dp[i - 1][j - 1] + score, dp[i - 1][j] + gap, dp[i][j - 1] + gap
            dp[i][j] = max(0, diagonal, up, left)
            traceback[i][j] = (i - 1, j - 1) if dp[i][j] == diagonal else ((i - 1, j) if dp[i][j] == up else (i, j - 1))
            if dp[i][j] > max_score: max_score, max_pos = dp[i][j], (i, j)

    return dp, traceback, max_pos


#alignment
def traceback_alignment(seq1, seq2, traceback, global_alignment=True, max_pos=None):
    aligned_seq1, aligned_seq2 = [], []
    i, j = (len(seq1), len(seq2)) if global_alignment else max_pos

    while i > 0 or j > 0:
        prev_i, prev_j = traceback[i][j]
        if (prev_i, prev_j) == (i - 1, j - 1): 
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
        elif prev_i == i - 1: 
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
        elif prev_j == j - 1: 
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
        i, j = prev_i, prev_j
        if not global_alignment and dp[i][j] == 0: break

    return "".join(reversed(aligned_seq1)), "".join(reversed(aligned_seq2))


#path
def traceback_path(traceback, global_alignment=True, max_pos=None):
    path = []
    i, j = (len(traceback) - 1, len(traceback[0]) - 1) if global_alignment else max_pos

    while i > 0 or j > 0:
        path.append((i, j))
        i, j = traceback[i][j]
        if not global_alignment and dp[i][j] == 0: break
    path.append((i, j))
    return path


#used html to create table since pandas table does not accept duplicates on headers, hence it would not look nice
def create_html_table(dp, seq1, seq2, path):
    html = "<table border='1' style='border-collapse: collapse;'>"
    html += "<tr><th style='background-color: lightblue;'>-</th><th style='background-color: lightblue;'>-</th>" 
    html += ''.join([f"<th style='background-color: green;'>{char}</th>" for char in seq2]) + "</tr>"
    
    for i in range(len(seq1) + 1):
        html += "<tr>"
        
        if i == 0:
            html += f"<th style='background-color: lightblue;'>-</th>"
        else:
            html += f"<th style='background-color: red;'>{seq1[i - 1]}</th>"
        
        # Matrix value
        for j in range(len(seq2) + 1):
            color = 'yellow' if (i, j) in path else 'white'
            html += f"<td style='background-color: {color}; color: black'>{dp[i][j]}</td>"
        html += "</tr>"
    
    html += "</table>"
    return html


#streamlit start
st.title("Lab 3 - Ben Lim Choong Chuen")
st.subheader("Simple Sequence Alignment Tool")

seq1 = st.text_input("Enter Sequence 1:")
seq2 = st.text_input("Enter Sequence 2:")

sel1, sel2 = st.columns([1,1])


alignment_type = st.radio("Choose Alignment Type:", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"])

with sel1:
    with st.expander("Set Alignment Parameters", expanded=False):
        match_score = st.number_input("Match Score:", value=1, step=1, help="Score for a match between two characters.")
        mismatch_score = st.number_input("Mismatch Penalty:", value=-1, step=1, help="Penalty for a mismatch between two characters.")
        gap_penalty = st.number_input("Gap Penalty:", value=-2, step=1, help="Penalty for a gap in the alignment.")

#error check for inputs
if st.button("Align"):
    if not seq1 or not seq2:
        st.error("Both sequences must be provided for alignment. Please enter sequences in both boxes.")
    else:
        #user can change alignment param
        if alignment_type == "Global (Needleman-Wunsch)":
            dp, traceback = needleman_wunsch(seq1, seq2, match=match_score, mismatch=mismatch_score, gap=gap_penalty)
            path = traceback_path(traceback, global_alignment=True)
            aligned_seq1, aligned_seq2 = traceback_alignment(seq1, seq2, traceback, global_alignment=True)
        else:
            dp, traceback, max_pos = smith_waterman(seq1, seq2, match=match_score, mismatch=mismatch_score, gap=gap_penalty)
            path = traceback_path(traceback, global_alignment=False, max_pos=max_pos)
            aligned_seq1, aligned_seq2 = traceback_alignment(seq1, seq2, traceback, global_alignment=False, max_pos=max_pos)

        #get score
        score = dp[len(seq1)][len(seq2)] if alignment_type == "Global (Needleman-Wunsch)" else dp[max_pos]
        st.markdown(" ")  #new line
        st.subheader("Alignment Results")
        
        score_text = f"<span style='color:orange;'>{score}</span>"
        aligned_seq1_text = f"<span style='color:red;'>{aligned_seq1}</span>"
        aligned_seq2_text = f"<span style='color:green;'>{aligned_seq2}</span>"

        st.markdown(f"##### Score: {score_text}", unsafe_allow_html=True)
        st.markdown(f"##### Aligned Sequence 1 : {aligned_seq1_text}", unsafe_allow_html=True)
        st.markdown(f"##### Aligned Sequence 2 : {aligned_seq2_text}", unsafe_allow_html=True)

        tback_text = f"<span style='color:yellow;'>Traceback</span>"
        st.markdown(f"### <br>Scoring Matrix with {tback_text}", unsafe_allow_html=True)
        html_table = create_html_table(dp, seq1, seq2, path)
        st.markdown(html_table, unsafe_allow_html=True)