from pyscipopt import *
from sage.arith.functions import LCM_list

# author: Gennadiy Averkov
# author: Christopher Hojny

def pairwise_intersection(FA,FB):
    """
        FA - first list of polyhedra
        FB - second list of polyhedra
        returns A cap B for A in FA and B in FB
    """

    intersections = []
    for A in FA:
        for B in FB:
            inter = A.intersection(B)

            if not inter.is_empty():
                intersections.append(inter)
    return intersections

def remove_redundancies(list_of_polyhedra):
    """
        list of polyhedra
        returns list of inclusionwise maximal polyhedra
    """

    max_polyhedra = []

    while len(list_of_polyhedra) > 0:
        P = list_of_polyhedra.pop()
        vertices_P = P.vertices_list()
        rays_P = P.rays_list()
        lines_P = P.lines_list()

        found_superset = False

        for Q in list_of_polyhedra + max_polyhedra:
            # check whether Q contains P
            contained = True

            # check vertices
            for v in vertices_P:
                if not Q.contains(v):
                    contained = False
                    break

            # there exists a vertex of P not contained in Q -> P \nsubseteq Q
            if not contained:
                continue

            # if P is bounded, it is contained in Q
            if len(rays_P) == 0 and len(lines_P) == 0:
                found_superset = True
                break

            # if P is unbounded, check whether rec(P) \subseteq rec(Q)
            if len(Q.rays_list()) == 0 and len(Q.lines_list()) == 0:
                continue

            # build recession cone of Q, have to turn lines l into rays l and -l
            C = Cone(Q.rays_list() + Q.lines_list() + [list(map(lambda t: -t, x)) for x in Q.lines_list()])

            for r in rays_P + lines_P + [list(map(lambda t: -t, x)) for x in lines_P]:
                if not C.contains(r):
                    contained = False
                    break

            if contained:
                found_superset = True
                break

        if not found_superset:
            max_polyhedra.append(P)

    return max_polyhedra
            
            
def radial_cone(P,z):
    """
        P - polyhedron
        z - point
        returns z+cone(z-P) 
    """
    z=vector(z)
    return Polyhedron(vertices=[z],rays=[z-vector(v) for v in P.vertices()])

def lat_complement_of_polyhedron(P): 
    """
        P - lattice polyhedron
        return list of halfspaces that cover all lattice points outside the polyhedron P
    """
    complements=[]
    for ineq in P.inequalities_list():
        b=ineq[0]
        a=ineq[1:]
        # turn inequality into an "integer" inequality if necessary
        denoms = [Rational(ai).denominator() for ai in a] + [Rational(b).denominator()]
        denomlcm = lcm(denoms)
        
        # the inequality has the form b+a*x>=0 
        # and we turn it to b+a*x<0
        # or equivalently (if x is integer) to 1+b+a*x<=0
        # or equivanelty to -1-b-a*x>=0
        compl_ineq = [(-1-b) * denomlcm]+[-t * denomlcm for t in a]
        complements.append( Polyhedron(ieqs= [ compl_ineq ]) )
    return complements

def is_observer(P,z):
    """
        P - lattice polytope
        z - lattice point
        checks if z is an observer for the lattice polytope P
    """
    assert P.is_lattice_polytope()
    Q=Polyhedron(P.vertices_list() + [z])
    return P.integral_points_count()+1==Q.integral_points_count()


def parity(x):
    """
        x - lattice point
        returns parity vector corresponding to x
    """
    return tuple(map( lambda t: t%2 , x))

def is_parity_complete(X):
    """
        X - lattice polytope or list of integer points
        returns whether X is parity-complete
    """
    if hasattr(X,'integral_points'):
        X=X.integral_points()
    parities=[parity(x) for x in X]
    d=len(parities[0])
    return len(set(parities))==2**d

def observers_for_parity_complete_case(P):
    """"
        P - parity-complete lattice polytope
        returns list of observers of P
    """
    X=tuple(map(vector,P.integral_points()))
    candidates=set([tuple(2*x-y) for x in X for y in X if x!=y])
    return [z for z in candidates if z not in P and is_observer(P,z)]                

def observers_in_box(P, boxlen):
    """
        P      - lattice polytope
        boxlen - length of box
        returns all observers of P in [-boxlen, boxlen]^d
    """
    X=tuple(map(vector,P.integral_points()))
    boxv = []
    for i in range(P.ambient_dim()):
        ine = [0 for j in range(P.ambient_dim())]
        ine[i] = 1
        ine = [boxlen] + ine
        boxv.append(ine)
        ine = [0 for j in range(P.ambient_dim())]
        ine[i] = -1
        ine = [boxlen] + ine
        boxv.append(ine)
    box = Polyhedron(ieqs=boxv)
    candidates = box.integral_points()

    return [z for z in candidates if z not in P and is_observer(P,z)]                

def find_closest_integer_point_in_region(Q, P):
    """
    Q - polyhedron
    P - full-dimensional lattice polyhedron disjoint from Q
    returns integer point q in Q with smallest lattice height w.r.t. P or "empty" if Q has no integer point
    """
    assert( P.intersection(Q).is_empty() )
    assert( P.ambient_dim() == Q.ambient_dim() )

    # collect information about problem
    dim = P.ambient_dim()
    conss_P = P.inequalities_list()
    conss_Q = Q.inequalities_list()
    equations_Q = Q.equations_list()

    # create MIP to find closest point
    #
    # assume P := Ax <= b
    #        Q := Cx <= d
    #
    # min h
    #     Aq <= b + h
    #     Cq <= d
    #      q integral
    model = Model("closest point")

    # hide output
    model.hideOutput()

    # create variables
    qvars = [model.addVar(vtype = "I", name = "q%d" % i, obj = 0.0, lb = -infinity) for i in range(dim)]
    hvar = model.addVar(vtype = "C", name = "h", obj = 1.0, lb = 0.0)

    # create constraints

    # qvars is contained in Q, need to add constraints for both inequalities and equations of Q
    cnt = 0
    for cons in conss_Q:
        model.addCons( cons[0] + quicksum( cons[i + 1] * qvars[i] for i in range(dim) ) >= 0, name = "Qcons%d" % cnt )
        cnt += 1
    cnt = 0
    for eq in equations_Q:
        model.addCons( eq[0] + quicksum( eq[i + 1] * qvars[i] for i in range(dim) ) == 0, name = "Qeqs%d" % cnt )
        cnt += 1

    # qvars is contained in P expanded by lattice height h
    cnt = 0
    for cons in conss_P:
        model.addCons( hvar + cons[0] + quicksum( cons[i + 1] * qvars[i] for i in range(dim) ) >= 0, name = "Pcons%d" % cnt )
        cnt += 1

    # solve the problem
    model.optimize()

    # get optimal solution
    if model.getStatus() == "infeasible":
        return "empty", 0

    q = [model.getVal(qvars[i]) for i in range(len(qvars))]
    return q, model.getVal(hvar)


def get_region_with_minimal_lattice_height(Regions, P):
    """ 
        Regions - list of polyhedral regions disjoint from P
        P       - full-dimensional lattice polyhedron that has finitely many observers 
        returns a lattice point z with minimal lattice height h contained on one of the search regions
    """
    minlh = infinity
    minregion = None
    minz = None

    for R in Regions:
        # solve a MIP to find closest integer point in region R
        z, h = find_closest_integer_point_in_region(R, P)

        if z == "empty":
            continue

        # if we have found a point with lattice height 1, this point is as close as possible
        if round(h) == 1:
            return R, z, 1
        elif round(h) < minlh:
            minlh = round(h)
            minregion = R
            minz = z

    return minregion, minz, minlh

def observers_through_direct_procedure(P):
    """ 
        P - full-dimensional lattice polyhedron that has finitely many observers 
        returns list of observers for (the set of integer points of) P
    """
    ObserversFound=[]
    SearchRegions=lat_complement_of_polyhedron(P)
    # Invariant of this iterative: every observer is either in ObserversFound or in one of the SearchRegions
    cnt = 0
    while len(SearchRegions)>0:
        # find a lattice point z with minimal lattice height h contained on one of the search regions
        Region, z, h = get_region_with_minimal_lattice_height(SearchRegions, P)

        # all regions are lattice-free
        if h == infinity:
            break
        SearchRegions.remove(Region)

        # make z integral and turn z into a rational vector
        for i in range(len(z)):
            z[i] = Rational(round(z[i] + 0.25))
        assert(is_observer(P,z))

        # we have found a new observer
        ObserversFound.append(z)

        # remove its special cone from the search regions intersecting the cone
        cone = radial_cone(P,z)
        compl=lat_complement_of_polyhedron(cone)
        regions_intersecting_cone = [S for S in SearchRegions if not S.intersection(cone).is_empty()]
        for S in regions_intersecting_cone:
            SearchRegions.remove(S)
        
        regions_intersecting_cone.append(Region)
        reduced_regions = pairwise_intersection(regions_intersecting_cone, compl)
        SearchRegions.extend(reduced_regions)

        SearchRegions = remove_redundancies(SearchRegions)

    return ObserversFound


