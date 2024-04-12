#%%
import sympy as sp
from sympy.combinatorics import Permutation
#%% 
# Set up the inputs.

eta = [0,1,0,0,1,0]
tau: Permutation = Permutation([1,2,0,3,4,5])

CGL_LENGTH = len(eta)

# Compute $\tau_\bullet$.
level_sets = [[] for _ in range(CGL_LENGTH)]
for i, v in enumerate(eta):
    level_sets[v].append(i)
tau_bullet = [0 for _ in range(CGL_LENGTH)]
tau_inv_level_sets = [[] for _ in range(CGL_LENGTH)]
for i, l in enumerate(level_sets):
    for j in l:
        tau_inv_level_sets[i].append(tau.__invert__()(j))
for l in  tau_inv_level_sets:
    l.sort()
for l in level_sets:
    l.sort()
for i, l in enumerate(level_sets):
    for j,v in enumerate(l):
        tau_bullet[int(tau.apply(tau_inv_level_sets[i][j]))] = int(v)

# Number of elements with a predecessor.
CGL_RANK = CGL_LENGTH - sum(map(lambda l: int(len(l) > 0), level_sets))

def to_one_indexed(l: list[int]) -> list[int]:
    return [i + 1 for i in l]
print(f"$\\tau = {to_one_indexed(tau.list())}$")
print(f"$\\tau_{{\\bullet}} = {to_one_indexed(tau_bullet)}$")
# Permutations are composed right to left in sympy!
print(f"$\\tau_{{\\bullet}}\\tau = {to_one_indexed((tau * Permutation(tau_bullet)).list())}$") # type: ignore
#%%


# The transpose of T and TB is taken because sympy works with permutations applied from right to left!
T: sp.Matrix = sp.PermutationMatrix(tau).as_explicit().transpose()
TB: sp.Matrix = sp.PermutationMatrix(Permutation(tau_bullet)).as_explicit().transpose()

# Add the conditions to make it behave like nonroot of unity
q = sp.Symbol('q', positive=True, noninteger=True)

bbl: sp.Matrix = sp.Matrix([
    [1,q,q**(-1),q,q**(-1),q**(-2)],
    [q**(-1), 1, q, q**(-1), q**(-2), q**(-1)],
    [q, q**(-1), 1, q**(-2), q**(-1), q],
    [q**(-1), q, q**2, 1, q, q**(-1)],
    [q, q**2, q, q**(-1), 1, q],
    [q**2, q, q**(-1), q, q**(-1), 1]
    ])
# Permute according to tau.
# Note that sympy applies permutations in the reverse order!
bbl = T.transpose() * bbl * T

# The "exponent matrix" of the $\lambda$ matrix
L: sp.Matrix = bbl.applyfunc(lambda x: sp.log(x, q).simplify())
assert L.is_anti_symmetric()
#%%

# Want to make matrix with as columns $\overline{e}_i$


E: sp.Matrix = sp.Matrix(CGL_LENGTH, CGL_LENGTH, lambda i,j: 0 if i > j else int(eta[tau(i)] == eta[tau(j)])) # type: ignore

#%%

# The matrix $A$ is the "exponent matrix" of the $\alpha$ matrix
A: sp.Matrix =   L * E 
# sp.pretty_print(A)

# The matrix $\widehat{Q}$ is the "exponent matrix" of the $\widehat{\mathbf{q}}$ matrix
QHAT: sp.Matrix =  E.transpose() * A 
assert QHAT.is_anti_symmetric()

Q:sp.Matrix = (TB * T) * QHAT * (TB * T)**(-1)
assert Q.is_anti_symmetric()

#%% 

D: sp.Matrix = sp.Matrix(CGL_LENGTH, CGL_RANK, lambda i, j: 2 if i == j else 0)
B: sp.Matrix = -Q.inv() * D
print(sp.latex(B))
# %%
