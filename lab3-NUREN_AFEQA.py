import streamlit as st
import numpy as np
import pandas as pd

# Title and description
st.title("Sequence Alignment Application")

# Input section
st.subheader("Input Sequences")
st.write("Please enter two sequences and select the type of alignment to proceed.")
# Input sequences
seq1 = st.text_input("Enter Sequence 1:")
seq2 = st.text_input("Enter Sequence 2:")

# Alignment type selection
alignment_type = st.selectbox("Choose alignment type:", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"])

# Process button
if st.button("Align Sequences"):
    if not seq1 or not seq2:
        st.error("Please enter both sequences.")
    else:
        # Hardcoded scoring parameters
        match = 1
        mismatch = -1
        gap = -2

        # Helper function to create unique labels
        def make_unique_labels(sequence):
            """Add an index to each character to make labels unique."""
            unique_labels = []
            counts = {}
            for char in sequence:
                if char not in counts:
                    counts[char] = 0
                counts[char] += 1
                unique_labels.append(f"{char}{counts[char]}")
            return unique_labels

        def initialize_matrix(n, m, global_align=True):
            """Initialize alignment matrix."""
            matrix = np.zeros((n, m), dtype=int)
            if global_align:
                for i in range(1, n):
                    matrix[i][0] = matrix[i - 1][0] + gap
                for j in range(1, m):
                    matrix[0][j] = matrix[0][j - 1] + gap
            return matrix

        def calculate_score(matrix, global_align=True):
            """Fill the alignment matrix and store the backtrace."""
            backtrace = np.zeros_like(matrix, dtype=int)
            max_score, max_pos = 0, (0, 0)
            
            for i in range(1, len(seq1) + 1):
                for j in range(1, len(seq2) + 1):
                    match_score = matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
                    delete_score = matrix[i - 1][j] + gap
                    insert_score = matrix[i][j - 1] + gap
                    
                    if global_align:
                        matrix[i][j] = max(match_score, delete_score, insert_score)
                    else:
                        matrix[i][j] = max(0, match_score, delete_score, insert_score)
                    
                    if not global_align and matrix[i][j] >= max_score:
                        max_score = matrix[i][j]
                        max_pos = (i, j)
                        
                    if matrix[i][j] == match_score:
                        backtrace[i][j] = 1  # Diagonal
                    elif matrix[i][j] == delete_score:
                        backtrace[i][j] = 2  # Up
                    else:
                        backtrace[i][j] = 3  # Left
            
            return matrix, backtrace, max_pos if not global_align else None

        def traceback(matrix, backtrace, max_pos=None, global_align=True):
            """Perform the traceback."""
            aligned_seq1, aligned_seq2 = "", ""
            i, j = max_pos if not global_align else (len(seq1), len(seq2))
            
            while (i > 0 or j > 0) if global_align else matrix[i][j] != 0:
                if backtrace[i][j] == 1:
                    aligned_seq1 = seq1[i - 1] + aligned_seq1
                    aligned_seq2 = seq2[j - 1] + aligned_seq2
                    i, j = i - 1, j - 1
                elif backtrace[i][j] == 2:
                    aligned_seq1 = seq1[i - 1] + aligned_seq1
                    aligned_seq2 = "-" + aligned_seq2
                    i -= 1
                elif backtrace[i][j] == 3:
                    aligned_seq1 = "-" + aligned_seq1
                    aligned_seq2 = seq2[j - 1] + aligned_seq2
                    j -= 1
            
            return aligned_seq1, aligned_seq2

        # Alignment logic
        n, m = len(seq1) + 1, len(seq2) + 1
        global_align = alignment_type == "Global (Needleman-Wunsch)"
        
        # Initialize and calculate alignment matrix
        matrix = initialize_matrix(n, m, global_align)
        matrix, backtrace, max_pos = calculate_score(matrix, global_align)
        aligned_seq1, aligned_seq2 = traceback(matrix, backtrace, max_pos, global_align)
        
        # Add unique labels to rows and columns
        unique_seq1_labels = make_unique_labels(seq1)
        unique_seq2_labels = make_unique_labels(seq2)

        labeled_matrix = pd.DataFrame(
            matrix,
            index=["-"] + unique_seq1_labels,  # Add '-' for alignment with Sequence 1
            columns=["-"] + unique_seq2_labels  # Add '-' for alignment with Sequence 2
        )

        # Display outputs
        st.subheader("Alignment Results")
        st.write("Alignment Score:", matrix[max_pos] if not global_align else matrix[-1][-1])
        st.write("Aligned Sequences:")
        st.text("Sequence 1: " + aligned_seq1)
        st.text("Sequence 2: " + aligned_seq2)
        
        st.write("Alignment Matrix:")
        st.dataframe(labeled_matrix)
        
        # Bonus: Backtrace path visualization
        st.write("Backtrace Path:")
        path_matrix = pd.DataFrame(
            np.zeros_like(matrix, dtype=str),
            index=["-"] + unique_seq1_labels,
            columns=["-"] + unique_seq2_labels
        )
        i, j = max_pos if not global_align else (len(seq1), len(seq2))
        
        while (i > 0 or j > 0) if global_align else matrix[i][j] != 0:
            path_matrix.iloc[i, j] = "X"
            if backtrace[i][j] == 1:
                i, j = i - 1, j - 1
            elif backtrace[i][j] == 2:
                i -= 1
            elif backtrace[i][j] == 3:
                j -= 1
        
        st.write(path_matrix)