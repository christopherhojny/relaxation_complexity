# author: Matthias Schymura

## Simplex T3 and the list of its integer points

T3 = Polyhedron(ieqs = [[1,-1,0,0],[1,0,-1,0],[1,0,0,-1],[0,1,1,1]])
L = list(T3.integral_points())

## removing the 4 points of Delta_3 = {0,e_1,e_2,e_3} from L

pt0 = Polyhedron([(0,0,0)])
pt1 = Polyhedron([(1,0,0)])
pt2 = Polyhedron([(0,1,0)])
pt3 = Polyhedron([(0,0,1)])
L.remove(pt0.integral_points()[0])
L.remove(pt1.integral_points()[0])
L.remove(pt2.integral_points()[0])
L.remove(pt3.integral_points()[0])

## augmenting L by all-one row to prepare 4x4 minor check

M = Matrix(L)
M = M.augment(Matrix([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]).transpose())
M = M.transpose()

## Delta_3 in homogeneous coordinates

D3 = Matrix([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
D3 = D3.augment(Matrix([1,1,1,1]).transpose())
D3 = D3.transpose()


def augmentT3(M, i):
    """
    M   - homogenized list of integer points in T3 without Delta_3
    i   - number of elements from L to augment Delta_3 with
    returns all augmentations of Delta_3 by i elements from L and with all 4x4 minors in {-1,0,1}
    """
    
    X_list = []
    C = Combinations(range(M.ncols()),i)
    
    for c in C :
      A = D3.augment(M.matrix_from_columns(c))
      mns = A.minors(4)
      if max(mns) == 1 and min(mns) >= -1 :
        X_list.append(A)
    
    return X_list

## compute all possibilities to augment L with 1 , 2 or 3 elements from T3 \ Delta_3

X1 = augmentT3(M, 1)
X2 = augmentT3(M, 2)
X3 = augmentT3(M, 3)
    
#####################################################################################################
# affine normal form (c/ Chris Borger) for checking
# affine unimodular equivalence
# source: https://github.com/christopherborger/mixed_volume_classification/blob/master/polytopes.sage
#####################################################################################################
    
def affine_normal_form(P):
    V=[vector(v) for v in P.vertices()]
    s=sum(V)
    N=len(V)
    Q=N*P-s
    Q=Polyhedron(Q.lattice_polytope().normal_form())
    Q=Q-vector(min(Q.vertices()))
    return Q/N

## compute affine normal forms of polytopes corresponding to sets in X1

AFN1 = []
for X in X1 :
    m = X.matrix_from_rows([0,1,2])
    P = Polyhedron(m.transpose())
    AFN1.append(affine_normal_form(P))

S1 = set(AFN1)
print("\tRESULT: there are %d augmentations of Delta_3 by 1 element of T3 \\ Delta_3." % len(X1))
print("\tRESULT: %d among them are pairwise unimodularly inequivalent.\n" % len(S1))
print("\trepresentatives are")
print(X1[0].matrix_from_rows([0,1,2]), "\n")

## compute affine normal forms of polytopes corresponding to sets in X2

AFN2 = []
for X in X2 :
    m = X.matrix_from_rows([0,1,2])
    P = Polyhedron(m.transpose())
    AFN2.append(affine_normal_form(P))

S2 = set(AFN2)
print("\tRESULT: there are %d augmentations of Delta_3 by 2 elements of T3 \\ Delta_3." % len(X2))
print("\tRESULT: %d among them are pairwise unimodularly inequivalent.\n" % len(S2))
print("\trepresentatives are")
print(X2[0].matrix_from_rows([0,1,2]), "\n")
print(X2[2].matrix_from_rows([0,1,2]), "\n")

## no three elements from L can be used to extend Delta_3

if len(X3) == 0:
    print("\tRESULT: no three integer points from T3 \\ Delta_3 can be added to Delta_3 without introducing a simplex of normalized volume >= 2.")
else:
    print("\tERROR: there is a triple of integer points from T3 \\ Delta_3 that can be added to Delta_3 without introducing a simplex of normalized volume >= 2.")

