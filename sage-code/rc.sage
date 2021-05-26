from pyscipopt import *
from itertools import *
import time

# author: Christopher Hojny

def extract_solution(model, avars, bvars, ninequs, dim, usecons):
    """
    model   - optimization model to find rc
    avars   - variables encoding left-hand side coefficients
    bvars   - variables encoding rhs
    dim     - dimension of ambient space of lattice-convex set
    usecons - list encoding which inequalities have been used in the formulation
    returns inequality system defining relaxation in sage format
    """

    solution = []
    for i in range(ninequs):
        if usecons[i] == 0:
            continue
        inequality = [model.getVal(bvars[i])] + [-model.getVal(avars[i][j]) for j in range(dim)]
        solution.append(inequality)

    return solution


def verify_solution(P, relaxinequs):
    """
    P           - integral polyhedron
    relaxinequs - system to be checked to define a relaxation
    returns whether relaxinequs define a relaxation of X
    """

    relax = Polyhedron(ieqs=relaxinequs)

    if P.integral_points_count() == relax.integral_points_count():
        return True
    print("relaxation has %d integer points, but P just %d" % (relax.integral_points_count(), P.integral_points_count()))
    return False


def scale_inequalities(ineqs):
    """
    ineqs - list of inequalities
    returns a list of equivalent inequalities (but with better scaling)
    """

    newineqs = []
    for ineq in ineqs:
        minentry = infinity
        maxentry = 0
        newineq = [i for i in ineq]
        for j in range(len(ineq)):
            newineq[j] = round(newineq[j], 6)
            if minentry > abs(newineq[j]):
                minentry = abs(newineq[j])
            if maxentry < abs(newineq[j]):
                maxentry = abs(newineq[j])

        # scale inequality
        if maxentry / minentry <= 100 and minentry < 1:
            for j in range(len(newineq)):
                newineq[j] = round(newineq[j] / minentry,6)

        for j in range(len(newineq)):
            newineq[j] = Rational(newineq[j])
        newineqs.append(newineq)

    return newineqs



def has_integer_formulation(X, Y, ninequs, conflictcuts=[], feasibilitymodel=True, returnsolution=False, handlesymmetries=True):
    """
        X                - set of integer points
        Y                - set of integer points disjoint from X
        ninequs          - number of inequalities in integer formulation
        conflictcuts     - indices of pairs of entries of Y which cannot be separated simultaneously by one inequality
        feasibilitymodel - if True, we use a feasibility based model, otherwise an optimization based model
        handlesymmetries - whether symmetries shall be handled
        returns True if there exists a set of ninequs linear inequalities separating X and Y
    """

    # find dimension of problem and maximal sup-norm of point in X
    dim = len(X[0])
    maxval = -infinity
    for i in range(len(X)):
        for j in range(dim):
            if maxval < abs(X[i][j]):
                maxval = abs(X[i][j])

    model = Model("formulationExists")

    model.hideOutput()

    model.setObjective(0, "maximize")

    # create variables
    if feasibilitymodel:
        muvar = model.addVar(vtype="C", lb=0.0001, ub=1.0, obj=0.0, name="mu")
    else:
        muvar = model.addVar(vtype="C", lb=0.0, ub=1.0, obj=1.0, name="mu")
    avars = [[model.addVar(vtype="C", lb=-1.0, ub=1.0, obj=0.0, name="a#%d#%d" % (i,j)) for j in range(dim)] for i in range(ninequs)]
    bvars = [model.addVar(vtype="C", lb=-dim*maxval, ub=dim*maxval, obj=0.0, name="b#%d" % i) for i in range(ninequs)]
    svars = [[model.addVar(vtype="B", obj=0.0, name="s#%d#%d" % (i,j)) for j in range(ninequs)] for i in range(len(Y))]

    # create constraints

    # each inequality (a,b) is satisfied by all points in X
    for i in range(ninequs):
        for x in range(len(X)):
            model.addCons(quicksum(avars[i][j] * X[x][j] for j in range(dim)) <= bvars[i], name="validinequ#%d#%d" % (i,x))

    # for each y in Y, there exists at least one constraint that is violated
    for y in range(len(Y)):
        model.addCons(quicksum(svars[y][i] for i in range(ninequs)) >= 1, name="packs#%d" % y)

    # link constraint-violation-indicator svars with constraints
    #
    # To avoid big-M constraints, we use an indicator constraint: If binvar=1, then the given constraint has to be satisfied
    for y in range(len(Y)):
        for i in range(ninequs):
            model.addConsIndicator(quicksum(avars[i][j] * Y[y][j] for j in range(dim)) >= bvars[i] + muvar, binvar=svars[y][i], name="violated#%d#%d" % (y,i))

    # handle symmetries
    if handlesymmetries :

        # first constraint has sorted coefficients
        for j in range(dim-1):
            model.addCons( avars[0][j] >= avars[0][j + 1], name="sym1#%d" % j)

        # coefficients are sorted w.r.t. their first coefficient
        for i in range(ninequs - 1):
            model.addCons( avars[i][1] >= avars[i + 1][1], name="sym2#%d" % i)

    P = Polyhedron(X)

    # add hiding set cuts
    for i in range(len(conflictcuts)):
        p1 = conflictcuts[i][0]
        p2 = conflictcuts[i][1]

        for k in range(ninequs):
            model.addCons( svars[p1][k] + svars[p2][k] <= 1, name="hidingset#%d#%d#%d" % (p1,p2,k), initial=False)

    model.data = ninequs, X, Y, P, svars

    model.hideOutput()
    model.optimize()
    model.printStatistics()

    returncode = False

    if feasibilitymodel:
        if model.getStatus() == "optimal":
            returncode = True
    else:
        if model.getStatus() == "infeasible" or model.getStatus() == "inforunbd":
            returncode =  False
        elif model.getObjVal() > 0:
            returncode = True

    if returnsolution:
        if returncode:
            return True, extract_solution(model, avars, bvars, ninequs, dim, [1 for i in range(ninequs)])
        else:
            return False, []
    else:
        return returncode


def find_minimum_formulation(X, Y, nmininequs, nmaxinequs, conflictcuts=[], returnsolution=False, handlesymmetries=True):
    """
        X                - set of integer points
        Y                - set of integer points disjoint from X
        nmininequs       - minimum number of inequalities in integer formulation
        nmininequs       - maximum number of inequalities in integer formulation
        handlesymmetries - whether symmetries shall be handled
        returns rc
    """

    # find dimension of problem and maximal sup-norm of point in X
    dim = len(X[0])
    maxval = -infinity
    for i in range(len(X)):
        for j in range(dim):
            if maxval < abs(X[i][j]):
                maxval = abs(X[i][j])


    model = Model("formulationExists")

    model.hideOutput()

    model.setObjective(0, "minimize")

    # create variables
    avars = [[model.addVar(vtype="C", lb=-1.0, ub=1.0, obj=0.0, name="a#%d#%d" % (i,j)) for j in range(dim)] for i in range(nmaxinequs)]
    bvars = [model.addVar(vtype="C", lb=-dim*maxval, ub=dim*maxval, obj=0.0, name="b#%d" % i) for i in range(nmaxinequs)]
    svars = [[model.addVar(vtype="B", obj=0.0, name="s#%d#%d" % (i,j)) for j in range(nmaxinequs)] for i in range(len(Y))]
    zvars = [model.addVar(vtype="B", obj=1.0, name="z#%d" % i) for i in range(nmaxinequs)]

    # create constraints

    # each inequality (a,b) is satisfied by all points in X
    for i in range(nmaxinequs):
        for x in range(len(X)):
            model.addCons(quicksum(avars[i][j] * X[x][j] for j in range(dim)) <= bvars[i], name="validinequ#%d#%d" % (i,x))

    # for each y in Y, there exists at least one constraint that is violated
    for y in range(len(Y)):
        model.addCons(quicksum(svars[y][i] for i in range(nmaxinequs)) >= 1, name="packs#%d" % y)

    # link constraint-violation-indicator svars with constraints
    #
    # To avoid big-M constraints, we use an indicator constraint: If binvar=1, then the given constraint has to be satisfied
    for y in range(len(Y)):
        for i in range(nmaxinequs):
            model.addConsIndicator(quicksum(avars[i][j] * Y[y][j] for j in range(dim)) >= bvars[i] + 0.001, binvar=svars[y][i], name="violated#%d#%d" % (y,i))

    # link svars and zvars
    for c in range(nmaxinequs):
        for y in range(len(Y)):
            model.addCons(svars[y][c] <= zvars[c], name="s%dz%d" % (y,c))

    # add constraint for minimum number of inequalities
    model.addCons( quicksum(zvars[c] for c in range(nmaxinequs)) >= nmininequs, "minnconss")

    # handle symmetries
    if handlesymmetries :
        for j in range(nmininequs-1):
            model.addCons( zvars[j] >= zvars[j + 1], name="sym0#%d" % j)

        # coefficients are sorted w.r.t. their first coefficient
        for i in range(nmaxinequs - 1):
            model.addCons( avars[i+1][0] + 2*zvars[i+1] <= avars[i][0] + 2*zvars[i], name="sym2#%d" % i)

        # inequality becomes trivial if not used
        for i in range(nmaxinequs):
            for j in range(dim):
                model.addCons( avars[i][j] <= uvars[i], name="trivialA#%d#%d" % (i,j) )
                model.addCons( avars[i][j] >= -uvars[i], name="trivialB#%d#%d" % (i,j) )
            model.addCons( dim * maxval <= bvars[i] + 2 * maxval * uvars[i], name="trivialC#%d" % i )

    P = Polyhedron(X)

    # add hiding set cuts
    for i in range(len(conflictcuts)):
        p1 = conflictcuts[i][0]
        p2 = conflictcuts[i][1]

        for k in range(nmaxinequs):
            model.addCons( svars[p1][k] + svars[p2][k] <= 1, name="hidingset#%d#%d#%d" % (p1,p2,k), initial=False)

    model.data = nmaxinequs, X, Y, P, svars

    model.optimize()
    model.printStatistics()

    if returnsolution:
        usecons = [0 for i in range(nmaxinequs)]
        for c in range(nmaxinequs):
            if model.getVal(zvars[c]) > 0.5:
                usecons[c] = 1
        return round(model.getObjVal()), extract_solution(model, avars, bvars, nmaxinequs, dim, usecons)
    else:
        return round(model.getObjVal())


def find_max_hidingset(P, observers):
    """
        P         - full-dimensional lattice polyhedron that has finitely many observers
        observers - list of observers of P
        returns the maximum size of a hiding set
    """

    conflicts = []
    for i in range(len(observers)):
        for j in range(i+1, len(observers)):
            line = Polyhedron([observers[i], observers[j]])
            if P.intersection(line).is_empty():
                conflicts.append([i,j])

    model = Model("hidingSets")
    model.setObjective(0, "maximize")

    xvars = [model.addVar(vtype="B", obj=1.0, name="x#%d" % i) for i in range(len(observers))]

    for conflict in conflicts:
        model.addCons(xvars[conflict[0]] + xvars[conflict[1]] <= 1, name="conflict%d#%d" % (conflict[0],conflict[1]), initial=False)

    # add conflict cliques based on facet defining inequalities
    facets = P.inequalities_list()
    for f in range(len(facets)):
        cutpoints = []

        # find points that are separated by the facet defining inequality
        for i in range(len(observers)):
            val = facets[f][0];
            for j in range(P.ambient_dim()):
                val += observers[i][j] * facets[f][j + 1]

            if val < 0:
                cutpoints.append(i)

        # all points in cutpoints are separated by facet defining inequality: at most one of them can
        # be contained in hiding set
        model.addCons( quicksum(xvars[point] for point in cutpoints) <= 1, name="facetcut#%d" % f )

    model.hideOutput()
    model.optimize()

    return round(model.getObjVal())



def computeRC(P, observers, silent=True, returnsolution=False):
    """
        P              - full-dimensional lattice polyhedron that has finitely many observers
        observers      - list of observers of P
        silent         - whether different steps shall be displayed
        returnsolution - whether also a relaxation shall be returned
        returns the relaxation complexity of P (and a relaxation)
    """

    X = P.vertices_list()

    # if dimension <= 4, we know all integer formulations are bounded
    rc_lb = 1
    if P.ambient_dim() <= 4:
        rc_lb = P.ambient_dim() + 1

    if not silent:
        print("Start binary search with lower bound %d and upper bound %d." % (rc_lb, rc_ub))
    cnt = 0

    # compute conflicts based on hiding sets
    conflicts = []
    for i in range(len(observers)):
        for j in range(i+1, len(observers)):

            # if observers[i] and observers[j] form a hiding set, no inequality can cut off both simultaneously
            points = [observers[i], observers[j]]
            lineseg = Polyhedron(points)
            if not P.intersection(lineseg).is_empty():
                conflicts.append([i,j])

    # perform a binary search to check whether an integer formulation with test_rc many inequalities exist
    while rc_lb < rc_ub:

        test_rc = int((rc_ub + rc_lb)/2)

        if returnsolution:
            formulationexists, solution = has_integer_formulation(X, observers, test_rc, conflictcuts=conflicts, returnsolution=True)
        else:
            formulationexists = has_integer_formulation(X, observers, test_rc, conflictcuts=conflicts)

        if formulationexists:
            rc_ub = test_rc
            if returnsolution:
                foundsolution = solution
            if not silent:
                print("Iteration %d: update ub. New range [%d, %d]." % (cnt, rc_lb, rc_ub))
        else:
            rc_lb = test_rc + 1
            if not silent:
                print("Iteration %d: update lb. New range [%d, %d]." % (cnt, rc_lb, rc_ub))
        cnt += 1

    if returnsolution:
        return rc_lb, foundsolution
    return rc_lb


def computeRCboxapprox(P, silent=True, returnsolution=False):
    """
        P - full-dimensional lattice polyhedron
        returns the relaxation complexity of P computed by iteratively increasing a bounding box for P
    """

    X = P.vertices_list()

    # if dimension <= 4, we know all integer formulations are bounded
    rc_lb = 1
    if P.ambient_dim() <= 4:
        rc_lb = P.ambient_dim() + 1
    rc_ub = P.n_facets()

    # can we stop early?
    if rc_lb == rc_ub:
        print("rc is %d because of matching upper and lower bounds" % rc_ub)
        if returnsolution:
            return rc_ub, []
        return rc_ub

    l1_radius = 2

    starttime = time.time()
    while True:

        boxvertices = []
        for i in range(P.ambient_dim()):
            vertex = [0 for j in range(P.ambient_dim())]
            vertex[i] = l1_radius
            boxvertices.append(vertex)

            vertex = [0 for j in range(P.ambient_dim())]
            vertex[i] = -l1_radius
            boxvertices.append(vertex)
        box = Polyhedron(boxvertices)

        # compute the infeasible points within l1_radius
        Y = []
        for x in X:
            for y in box.integral_points():
                elem = [x[i] + y[i] for i in range(P.ambient_dim())]
                if not P.contains(elem)  and not elem in Y:
                    Y.append(elem)

        # possibly generate conflicts
        conflicts = []
        for i in range(min(100,len(Y))):
            for j in range(i+1, min(100,len(Y))):
                if not (Polyhedron([Y[i], Y[j]]).intersection(P)).is_empty():
                    conflicts.append([i,j])

        # compute rc(X,Y)
        rc_guess, relax_guess = find_minimum_formulation(X, Y, rc_lb, rc_ub, conflictcuts=conflicts,
                                                         returnsolution=True, handlesymmetries=True)

        # we can stop if the we have found a relaxation
        relax = scale_inequalities(relax_guess)
        relax_polyhedron = Polyhedron(ieqs=relax)

        # update upper bound
        if relax_polyhedron.is_compact():
            if len(relax) < rc_ub and verify_solution(P, relax):
                rc_ub = len(relax)
        else:
            print("candidate relaxation is unbounded")

        # update lower bound
        if len(relax) > rc_lb:
            rc_lb = len(relax)

        if len(relax) >= rc_ub:
            print("box-rc equals upper bound")
            endtime = time.time()
            print("running time: %f" % (endtime - starttime))
            if returnsolution:
                return rc_guess, relax
            return rc_guess

        # the box is too small, no relaxation found so far
        l1_radius += 1
        print("box is too small: increase radius to %d" % l1_radius)

