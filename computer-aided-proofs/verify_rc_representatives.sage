# author: Christopher Hojny

def verify_certificate_simplex3():
    '''
    verifies that the relaxation complexity of
    X = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    is at least 4
    '''

    simplex3 = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])

    print("Step 1: Find a certificate for the claimed lower bound of 4")
    print("Test whether the following set is a certificate:")

    certificate = [(0, 1, 1), (1, 0, 1), (0, 2, 0), (1, 1, -1), (1, 1, -2), (0, 0, -1),
                   (-1, 0, 1), (-1, 0, 0), (2, 0, 0), (0, 0, 2), (2, 0, -1), (-1, 0, 2),
                   (0, -1, 1), (0, -1, 0), (1, 1, 1), (-1, 1, 0), (1, 1, 0), (0, 2, -1),
                   (0, -1, 2), (-2, 1, 1), (-1, 1, 1), (2, -1, 0), (-1, 2, 0), (1, 0, -1),
                   (1, -1, 0), (1, -1, 1), (1, -2, 1), (0, 1, -1)]
    print(certificate)

    # verify that none of the certificate points is contained in the simplex
    print("Verify that claimed certificate does not contain a point from the 3-simplex.")
    for point in certificate:
        if simplex3.contains(point):
            print("\tERROR: certificate point {} is contained in the 3-simplex.".format(point))
            return False
    print("\tSUCCESS: all certificate points are not contained in the 3-simplex.")

    edges = []
    for i in range(len(certificate)):
        for j in range(i+1, len(certificate)):
            Q = Polyhedron([certificate[i], certificate[j]])
            if not Q.intersection(simplex3).is_empty():
                edges.append([i,j])

    G = Graph(edges)

    # compute the chromatic number of the hiding graph G
    print("Compute the hiding graph's chromatic number.")
    print("If it is 4, we have found a finite certificate for rc(3-simplex) >= 4.")
    chromatic_number = G.chromatic_number()

    if chromatic_number == 4:
        print("\tSUCCESS: The hiding graph's chromatic number is 4.")
    else:
        print("\tERROR: The hiding graph's chromatic number is %d." % chromatic_number)
        return False

    # verify matching upper bound by providing a relaxation
    print("\nStep 2: The simplex itself is a certificate for the claimed upper bound of 4")

    return True


def verify_certificate_polytope1():
    '''
    verifies that the relaxation complexity of
    X = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1)]
    is at least 4
    '''

    X = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1)]
    P = Polyhedron(X)

    print("Step 1: Find a certificate for the claimed lower bound of 4")
    print("Test whether the following set is a certificate:")

    certificate = [(2, 0, -1), (-1, 1, 0), (-1, 1, 1), (0, 1, -1), (0, 1, 1), (1, -1, 0),
                   (1, 0, -1), (1, 0, 1), (1, 1, -1), (1, 1, 0), (2, 0, 0), (2, 1, -1),
                   (-2, 0, 1), (-2, 0, 2), (-2, 1, 1), (-1, -1, 2), (-1, 0, 0), (-1, 0, 2)]
    print(certificate)

    # verify that none of the certificate points is contained in the polytope
    print("Verify that claimed certificate does not contain a point from polytope 1.")
    for point in certificate:
        if P.contains(point):
            print("\tERROR: certificate point {} is contained in polytope 1.".format(point))
            return False
    print("\tSUCCESS: all certificate points are not contained in polytope 1.")

    edges = []
    for i in range(len(certificate)):
        for j in range(i+1, len(certificate)):
            Q = Polyhedron([certificate[i], certificate[j]])
            if not Q.intersection(P).is_empty():
                edges.append([i,j])

    G = Graph(edges)

    # compute the chromatic number of the hiding graph G
    print("Compute the hiding graph's chromatic number.")
    print("If it is 4, we have found a finite certificate for rc(polytope 1) >= 4.")
    chromatic_number = G.chromatic_number()

    if chromatic_number == 4:
        print("\tSUCCESS: The hiding graph's chromatic number is 4.")
    else:
        print("\tERROR: The hiding graph's chromatic number is %d." % chromatic_number)
        return False

    # verify matching upper bound by providing a relaxation
    print("\nStep 2: Find a certificate for the claimed upper bound of 4")
    relax = Polyhedron(ieqs=[[0,0.001,1,0.002],[1,-0.5005,-1,-1],[1,3,-1,2],[0.499,-0.333333,-0.166667,1]])

    # check whether the lattice points in the relaxation match with X
    if relax.integral_points_count() != len(X):
        print("\tERROR: polytope 1 has %d lattice points, but the claimed relaxation has %d" % (len(X), relax.integral_points_count()))
        return False

    for point in X:
        if not relax.contains(point):
            print("\tERROR: the claimed relaxation does not contain point {}".format(point))
            return False

    if relax.n_facets() != 4:
        print("\tERROR: the claimed relaxation is not simplicial")
    print("\tSUCCESS: polytope 1 has a relaxation with 4 facets, given by {}.".format(relax.inequalities_list()))


    return True

def verify_certificate_polytope2():
    '''
    verifies that the relaxation complexity of
    X = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(-1,1,0)]
    is at least 4
    '''

    X = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(-1,1,0)]
    P = Polyhedron(X)

    print("Step 1: Find a certificate for the claimed lower bound of 4")
    print("Test whether the following set is a certificate:")

    certificate = [(-2,0, 1), (0, -1, 1), (0, 1, -1), (0, 1, 1), (1, 1, 0), (2, -1, 1),
                   (2, 0, 0), (2, 1, -1), (3, 0, 1), (-2, 0, 0), (-2, 1, -1), (-2, 1, 0), (-2, 1, 1)]
    print(certificate)

    # verify that none of the certificate points is contained in the polytope
    print("Verify that claimed certificate does not contain a point from polytope 2.")
    for point in certificate:
        if P.contains(point):
            print("\tERROR: certificate point {} is contained in polytope 2.".format(point))
            return False
    print("\tSUCCESS: all certificate points are not contained in polytope 2.")

    edges = []
    for i in range(len(certificate)):
        for j in range(i+1, len(certificate)):
            Q = Polyhedron([certificate[i], certificate[j]])
            if not Q.intersection(P).is_empty():
                edges.append([i,j])

    G = Graph(edges)

    # compute the chromatic number of the hiding graph G
    print("Compute the hiding graph's chromatic number.")
    print("If it is 4, we have found a finite certificate for rc(polytope 2) >= 4.")
    chromatic_number = G.chromatic_number()

    if chromatic_number == 4:
        print("\tSUCCESS: The hiding graph's chromatic number is 4.")
    else:
        print("\tERROR: The hiding graph's chromatic number is %d." % chromatic_number)
        return False

    # verify matching upper bound by providing a relaxation
    print("\nStep 2: Find a certificate for the claimed upper bound of 4")
    relax = Polyhedron(ieqs=[[0,1,1,1],[0,0,0,1],[1,-0.24,-1,-1],[1,-1,1,-1]])

    # check whether the lattice points in the relaxation match with X
    if relax.integral_points_count() != len(X):
        print("\tERROR: polytope 2 has %d lattice points, but the claimed relaxation has %d" % (len(X), relax.integral_points_count()))
        return False

    for point in X:
        if not relax.contains(point):
            print("\tERROR: the claimed relaxation does not contain point {}".format(point))
            return False

    if relax.n_facets() != 4:
        print("\tERROR: the claimed relaxation is not simplicial")
    print("\tSUCCESS: polytope 2 has a relaxation with 4 facets, given by {}.".format(relax.inequalities_list()))

    return True

def verify_certificate_polytope3():
    '''
    verifies that the relaxation complexity of
    X = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(0,-1,1)]
    is at least 4
    '''

    X = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(0,-1,1)]
    P = Polyhedron(X)

    print("Step 1: Find a certificate for the claimed lower bound of 4")
    print("Test whether the following set is a certificate:")

    certificate = [(-1, -1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]
    print(certificate)

    # verify that none of the certificate points is contained in the polytope
    print("Verify that claimed certificate does not contain a point from polytope 3.")
    for point in certificate:
        if P.contains(point):
            print("\tERROR: certificate point {} is contained in polytope 3.".format(point))
            return False
    print("\tSUCCESS: all certificate points are not contained in polytope 3.")

    edges = []
    for i in range(len(certificate)):
        for j in range(i+1, len(certificate)):
            Q = Polyhedron([certificate[i], certificate[j]])
            if not Q.intersection(P).is_empty():
                edges.append([i,j])

    G = Graph(edges)

    # compute the chromatic number of the hiding graph G
    print("Compute the hiding graph's chromatic number.")
    print("If it is 4, we have found a finite certificate for rc(polytope 3) >= 4.")
    chromatic_number = G.chromatic_number()

    if chromatic_number == 4:
        print("\tSUCCESS: The hiding graph's chromatic number is 4.")
    else:
        print("\tERROR: The hiding graph's chromatic number is %d." % chromatic_number)
        return False

    # verify matching upper bound by providing a relaxation
    print("\nStep 2: Find a certificate for the claimed upper bound of 4")
    relax = Polyhedron(ieqs=[[0,2,1,2],[1,-0.499,-0.5,-1],[1.499,-0.5,-1,1],[3,-2,4,1]])

    # check whether the lattice points in the relaxation match with X
    if relax.integral_points_count() != len(X):
        print("\tERROR: polytope 3 has %d lattice points, but the claimed relaxation has %d" % (len(X), relax.integral_points_count()))
        return False

    for point in X:
        if not relax.contains(point):
            print("\tERROR: the claimed relaxation does not contain point {}".format(point))
            return False

    if relax.n_facets() != 4:
        print("\tERROR: the claimed relaxation is not simplicial")
    print("\tSUCCESS: polytope 3 has a relaxation with 4 facets, given by {}.".format(relax.inequalities_list()))

    return True


def verify_certificates():
    '''
    verifies that the relaxation complexities of
    X0 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)],
    X1 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1)],
    X2 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(-1,1,0)], and
    X3 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(0,-1,1)]
    is at least 4
    '''

    print("=====================================================================================================================")
    print("Verify certificate for standard simplex with vertices {}.".format([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))
    print("=====================================================================================================================")
    verify_certificate_simplex3()
    print("\n\n")

    print("=====================================================================================================================")
    print("Verify certificate for polytope 1 with vertices {}.".format([(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1)]))
    print("=====================================================================================================================")
    verify_certificate_polytope1()
    print("\n\n")

    print("=====================================================================================================================")
    print("Verify certificate for polytope 2 with vertices {}.".format([(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(-1,1,0)]))
    print("=====================================================================================================================")
    verify_certificate_polytope2()
    print("\n\n")

    print("=====================================================================================================================")
    print("Verify certificate for polytope 3 with vertices {}.".format([(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(0,-1,1)]))
    print("=====================================================================================================================")
    verify_certificate_polytope3()
