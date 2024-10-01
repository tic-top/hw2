## HW2
### 1. **Verification**
Sum and sum of square

**Example**:
> We verified the correctness of our parallel program by running both serial and parallel versions on small matrix sizes (e.g., 4x4, 8x8). The results of the decomposed matrices were compared, and they matched exactly. We also tested our implementation on edge cases, such as matrices filled with zeros and identity matrices.

### 2. **Matrix Decomposition Description**
Explain how you decomposed the matrix for parallel processing. If you used domain decomposition, be sure to specify the type (e.g., block decomposition, row-wise, column-wise).

**Example**:
> The matrix was decomposed using a 2D block decomposition approach. The matrix was divided into submatrices, each assigned to a separate processor. Each process was responsible for operating on its submatrix, and data exchanges occurred at the boundaries where necessary (e.g., sharing row and column data with neighboring processes). For a matrix of size \(m \times n\) and \(P\) processors, the matrix was divided into \(\frac{m}{\sqrt{P}} \times \frac{n}{\sqrt{P}}\) submatrices.

### 3. **Communication Description**
 Communication between processes was performed using point-to-point communication (MPI_Send and MPI_Recv). Each process exchanged its boundary rows and columns with its neighboring processes. The top row was sent to the process above, and the bottom row was sent to the process below. Similarly, left and right columns were exchanged. Diagonal communication was also performed for processes at the corners of the submatrices. The communication pattern ensured that each process had the required neighboring data for computation.

### 4. **Timing Results**
Provide the timing results of your program. Compare the time taken for different matrix sizes and varying numbers of processors.

### 5. **Analysis of Speedup and Scaling**
