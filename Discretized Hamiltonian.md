# Discretized Hamiltonian

$k\cdot p$ Hamiltonian from <https://aip.scitation.org/doi/abs/10.1063/1.4790804>

+ Basis we used: $\left\{|P1_z^+,\uparrow>,|P2_z^-,\uparrow>,|P1_z^+,\downarrow>,|P2_z^-,\downarrow>\right\}$

+ Origin Hamiltonian of 3D TI:

$$H(k)=\epsilon_0(k)+
\begin{pmatrix}
M(k) & A_1k_z & 0 & A_2k_-\\
A_1k_z & -M(k) & A_2k_- & 0 \\
0 & A_2k_+ & M(k) & -A_1k_z\\
A_2k_+ & 0 & -A_1k_z & -M(k)
\end{pmatrix}$$
,where $\epsilon(k)=C+D_1k_z^2+D_2k_\bot^2$, $M(k)=M-B_1k_z^2-B_2k_\bot^2$

+ Discretized $H$ in the $z$ direction, with discretization length $\Delta$. The diagonal term is:

$$
\begin{aligned}
h_{ii}(k) & =\epsilon(k_x,k_y) +U(z_i)+\frac {2G_{zz}}{\Delta^2}+\\ &
\begin{pmatrix}
M(k_x,k_y) & 0 & 0 & A_2k_-\\
0 & -M(k_x,k_y) & A_2k_- & 0 \\
0 & A_2k_+ & M(k_x,k_y) & 0\\
A_2k_+ & 0 & 0 & -M(k_x,k_y)
\end{pmatrix}
\end{aligned}
$$
and the coupling term between nearest layer:
$$
t_{i,i+1}=-\frac {G_{zz}}{\Delta^2}-\frac {G_{z}}{2\Delta}$$$$
t_{i,i-1}=-\frac {G_{zz}}{\Delta^2}+\frac {G_{z}}{2\Delta}
$$
,where
$$G_{zz}=
\begin{pmatrix}
D_1-B_1 & 0 & 0 & 0\\
0 & D_1+B_1 & 0 & 0 \\
0 & 0 & D_1-B_1 & 0\\
0 & 0 & 0 & D_1+B_1
\end{pmatrix}
$$
and
$$G_{z}=
\begin{pmatrix}
0 & +iA_1 & 0 & 0\\
+iA_1 & 0 & 0 & 0 \\
0 & 0 & 0 & -iA_1\\
0 & 0 & -iA_1 & 0
\end{pmatrix}
$$

+ The total Schrodinger Equation is $H_{ij}(k_x,k_y)\psi_j=E\psi_i$, and the large Hamiltonian is

$$H(k_x,k_y)=
\begin{pmatrix}
h_{1,1} & t_{1,2} & 0 & 0\\
t_{2,1} & h_{2,2} & t_{2,3} & 0 \\
0 & t_{3,2} & \ddots & \ddots\\
0 & 0 & \ddots & h_{N,N}
\end{pmatrix}
$$
