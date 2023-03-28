# author: Gennadiy Averkov

# In this code, a mixed-integer point is
# a point from set Z^{d-1} \times R.

def mixed_integer_points(Q):
    """
        For a given polytope Q of dimension Q
        returns the list of polyhedra
        whose union is the set of all 
        mixed-integer points of Q. 
    """
    assert ( Q.is_compact() )
    n = Q.ambient_dim()
    vertical = (n-1)*[0] + [1] 
    V = Q.vertices()
    Q_proj = Polyhedron(vertices = [v[:-1] for v in V])
    Z = Q_proj.integral_points()
    Z = [list(z)+[0] for z in Z]
    result = [] 
    for z in Z: 
        L = Polyhedron(vertices = [z], lines = [vertical])
        I = Q.intersection(L)
        result.append(I)
    return result

def the_defining_matrix(e):
    """
        Returns the inequality system
        describing polyhedron P_e.
    """
    a = AA(2).sqrt()
    return matrix(AA, [(0,1,0,0,-1), (0,0,0,1,-1), (1,-e,-1,(e-1)/(1+a), (e-1)*a/(1+a)), (0,0,0,1,a), (1,-1,1+e,-1,1)]  )

def to_radical_expression(entry):
    """
        Returns radical representation
        of a given number entry.
    """
    return entry.radical_expression()

######################
# MAIN PART OF PROOF
######################

# define the inequality system of polyhedron P_{\frac{1}{8}}
A = the_defining_matrix(1/8)

print("The system of inequalities is defined by the rows of this matrix")
print(A.apply_map(to_radical_expression))

# define polyhedron P_{\frac{1}{8}}
Q = Polyhedron(ieqs = A.rows())

# get the mixed-integer points in P_{\frac{1}{8}}
result = mixed_integer_points(Q)

print("Checking if the set of mixed-integer points in the respective polyhedron is finite:")
print(all([p.dim()==0 for p in result] ))

print("The mixed-integer solutions are the columns of this matrix:")
result = matrix([ p.vertices()[0] for p in result]).transpose()
result = result.apply_map(to_radical_expression)
print(result)
