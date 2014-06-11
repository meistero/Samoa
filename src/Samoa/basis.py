#!/usr/bin/python

from sympy import *
from sympy.plotting import plot3d

class basis:
    def __init__(self, order, domain):
        self.order = order
        self.domain = domain

        x, y = symbols('x,y', real=True)

        self.functions = [self.basis(order, index, domain, x, y) for index in xrange(0, (self.order + 1) * (self.order + 2) / 2)]

    def __repr__(self):
        if len(self.functions) == 0:
            return "basis()"

        s = "basis("
        for f in self.functions[:-1]:
            s += "%s, " % f
        s += "%s)" % self.functions[-1]

        return s

    def __getitem__(self, index):
        return self.functions[index]

class basis_function:
    def __init__(self, expr, domain):
        self.expr = expr
        self.domain = domain

    def __repr__(self):
        return "function(%s, %s)" % (self.expr, self.domain)

    def __getitem__(self, index):
        assert(isinstance(self.expr, tuple))
        return basis_function(self.expr[index], self.domain)

    def __add__(f, g):
        return f, g

    def __mul__(f, g):
        if isinstance(g, dirac_basis_function):
            return NotImplemented

        domain = intersect(f.domain, g.domain)

        return basis_function(f.expr * g.expr, domain)

    def grad(f):
        x, y = symbols('x,y', real=True)
        return basis_function((diff(f.expr, x), diff(f.expr, y)), f.domain)

    def inner(f, g):
        expr = sum([a * b for (a,b) in zip(f.expr, g.expr)])
        domain = intersect(f.domain, g.domain)

        return basis_function(expr, domain)

    def volume_integrate(f):
        x, alpha = symbols('x, alpha', real=True)

        F = basis_function((Integral(f.expr.subs({x:alpha}), (alpha, 0, x)), 0), f.domain)
        return F.boundary_integrate()

    def boundary_integrate(f):
        x, y = symbols('x, y', real=True)
        t = symbols('t', real=True, nonnegative=True)

        bnd_int = S(0)

        if isinstance(f.domain, Triangle) or isinstance(f.domain, Polygon):
            p = f.domain.vertices[-1]

            if (isinstance(f.expr, tuple)):
                for p_next in f.domain.vertices:
                    v = p_next - p
                    bnd_int = bnd_int + Integral((v.y * f.expr[0] - v.x * f.expr[1]).subs({x:p.x + t * v.x,y:p.y + t * v.y}), (t, 0, 1))
                    p = p_next
            else:
                for p_next in f.domain.vertices:
                    v = p_next - p
                    bnd_int = bnd_int + Integral(sqrt(v.x ** 2 + v.y ** 2) * f.expr.subs({x:p.x + t * v.x,y:p.y + t * v.y}), (t, 0, 1))
                    p = p_next
        elif isinstance(f.domain, tuple):
            p1 = f.domain[0]
            p2 = f.domain[1]

            if (isinstance(f.expr, tuple)):
                v = p2 - p1
                bnd_int = bnd_int + Integral((v.y * f.expr[0] - v.x * f.expr[1]).subs({x:p1.x + t * v.x,y:p1.y + t * v.y}), (t, 0, 1))
            else:
                v = p2 - p1
                bnd_int = bnd_int + Integral(sqrt(v.x ** 2 + v.y ** 2) * f.expr.subs({x:p1.x + t * v.x,y:p1.y + t * v.y}), (t, 0, 1))
        else:
            bnd_int = bnd_int + sum([boundary_integrate(basis_function(f.expr, subdomain)) for subdomain in f.domain])

        return bnd_int.doit()

class dirac_basis_function(basis_function):
    def __init__(self, expr, point, domain):
        self.expr = expr
        self.point = point
        self.domain = domain

    def __mul__(f, g):
        point = f.point
        domain = intersect(f.domain, g.domain)

        return dirac_basis_function(f.expr * g.expr, point, domain)

    def __rmul__(f, g):
        return f * g

    def grad(f):
        raise NotImplementedError("Gradient of Dirac function is not available")

    def inner(f, g):
        expr = sum([a * b for (a,b) in zip(f.expr, g.expr)])

        point = intersect(f.point, g.point)
        domain = intersect(f.domain, g.domain)

        return dirac_basis_function(expr, point, domain)

    def epsilon_area(f, epsilon):
        d = intersect(f.domain, Polygon(
            Point(f.point.x - epsilon / 2, f.point.y - epsilon / 2),
            Point(f.point.x + epsilon / 2, f.point.y - epsilon / 2),
            Point(f.point.x + epsilon / 2, f.point.y + epsilon / 2),
            Point(f.point.x - epsilon / 2, f.point.y + epsilon / 2)))

        b = basis_function(S(1), d)
        area = volume_integrate(b)

        return area / (epsilon ** 2)

    def volume_integrate(f):
        x,y = symbols('x,y', real=True)

        if isinstance(f.point, Point):
            #compute the area ratio covered by the domain around the dirac point
            epsilon = Rational(1, 1e99)
            ratio = f.epsilon_area(epsilon)

            #multiply by the value, this gives an additive integral
            return ratio * f.expr.subs({x:f.point.x, y:f.point.y})

        else:
            return 0

    def boundary_integrate(f):
        if isinstance(f.domain, Point) and (f.domain.x == 0 or f.domain.y == 0 or 1 - f.domain.x - f.domain.y == 0):
            return f.expr.subs({x:f.domain.x, y:f.domain.y})
        else:
            return 0

def diverg(F):
    return grad(F[0])[0] + grad(F[1])[1]

def dx(f):
    return grad(f)[0]

def dy(f):
    return grad(f)[1]

def grad(f):
    if isinstance(f, basis_function):
        return f.grad()
    else:
        return sum([grad(f_sub) for f_sub in f])

def inner(f, g):
    if isinstance(g, dirac_basis_function):
        return g.inner(f)
    elif isinstance(f, basis_function):
        return f.inner(g)
    else:
        return sum([f_sub.inner(g_sub) for (f_sub, g_sub) in zip(f, g)])

def volume_integrate(f):
    if isinstance(f, basis_function):
        return f.volume_integrate()
    else:
        return sum([volume_integrate(f_sub) for f_sub in f])

def boundary_integrate(f):
    if isinstance(f, basis_function):
        return f.boundary_integrate()
    else:
        #this is problematic. a boundary integral over a function
        #is the sum of subintegrals only if the function is continuous
        raise NotImplementedError("Boundary integral can not be split")

def intersect(A, B):
    isection = decompose(intersection(A, B))

    if isinstance(A, Polygon) or isinstance(A, Segment):
        v_B = decompose(B)

        for b in v_B:
            if A.encloses(b):
                isection.append(b)

    if isinstance(B, Polygon) or isinstance(B, Segment):
        v_A = decompose(A)

        for a in v_A:
            if B.encloses(a):
                isection.append(a)

    if isinstance(A, Point) and isinstance(B, Point):
        if A == B:
            isection.append(A)

    if len(isection) == 0:
        return []
    else:
        return convex_hull(*isection)

def decompose(A):
    if isinstance(A, Polygon):
        v_A = A.vertices
    elif isinstance(A, Segment):
        v_A = [A.p1, A.p2]
    elif isinstance(A, Point):
        v_A = [A]
    else:
        v_A = [x for a in A for x in decompose(a)]

    return v_A

def index_to_coords(order, index):
    assert(order >= 0)
    assert(index < (order + 1) * (order + 2) / 2)

    if index == 0:
        i = order
        j = 0
    elif index == 1:
        i = 0
        j = 0
    elif index == 2:
        i = 0
        j = order
    elif index < 3 + (order - 1):
        i = 0
        j = index - 2
    elif index < 3 + 2 * (order - 1):
        i = 1 + index - (3 + (order - 1))
        j = order - i
    elif index < 3 + 3 * (order - 1):
        i = index - 2 * order
        j = 0
    else:
        j = int(sqrt(2 * (index - 3 * order) + 0.25) - 0.5)
        i = index - 3 * order - (j * (j + 1)) / 2
        j = order - 3 - j
        j += 1
        i += 1

    return i, j

def dirac_test_function(order, index, triangle, x, y):
    if order == 0:
        point = triangle.vertices[1] + (triangle.vertices[0] - triangle.vertices[1]) / 3 + (triangle.vertices[2] - triangle.vertices[1]) / 3
    else:
        i,j = index_to_coords(order, index)
        point = triangle.vertices[1] + (triangle.vertices[0] - triangle.vertices[1]) * Rational(i,order) + (triangle.vertices[2] - triangle.vertices[1]) * Rational(j,order)

    return dirac_basis_function(S(1), point, triangle)

def fv_dual_basis_function(order, index, triangle, x, y):
    assert(order == 1)
    assert(index < 3)

    domains = [
        Triangle(triangle.vertices[0], (triangle.vertices[0] + triangle.vertices[2]) / 2, (triangle.vertices[0] + triangle.vertices[1]) / 2),
        Polygon(triangle.vertices[1], (triangle.vertices[0] + triangle.vertices[1]) / 2, (triangle.vertices[0] + triangle.vertices[2]) / 2, (triangle.vertices[1] + triangle.vertices[2]) / 2),
        Triangle(triangle.vertices[2], (triangle.vertices[1] + triangle.vertices[2]) / 2, (triangle.vertices[0] + triangle.vertices[2]) / 2)
        ]

    return basis_function(S(1), domains[index])

def lagrange_basis_function_cb(order, i, j, triangle, x, y):
    epsilon, xi = symbols('epsilon, xi')

    base_i = [S(1) for m in xrange(order + 1)]
    base_j = [S(1) for m in xrange(order + 1)]
    base_k = [S(1) for m in xrange(order + 1)]

    for m in xrange(0, order):
        base_i[m + 1] = base_i[m] * (epsilon - Rational(m, order))
        base_j[m + 1] = base_j[m] * (xi - Rational(m, order))
        base_k[m + 1] = base_k[m] * (Rational(order - m, order) - epsilon - xi)

    k = order - i - j
    expr = base_i[i] * base_j[j] * base_k[k]
    expr /= expr.subs({epsilon:Rational(i,order), xi:Rational(j,order)})

    p = triangle.vertices[1] + (triangle.vertices[0] - triangle.vertices[1]) * epsilon + (triangle.vertices[2] - triangle.vertices[1]) * xi
    repl_rule = solve([p.x - x, p.y - y], epsilon, xi)
    expr = expr.subs(repl_rule)

    return basis_function(expr, triangle)

def lagrange_basis_function(order, index, triangle, x, y):
    i,j = index_to_coords(order, index)

    return lagrange_basis_function_cb(order, i, j, triangle, x, y)

def hierarchical_basis_function(order, index, triangle, x, y):
    i,j = index_to_coords(order, index)
    divisor = gcd(order, gcd(i, j))

    return lagrange_basis_function_cb(Rational(order, divisor), Rational(i, divisor), Rational(j, divisor), triangle, x, y)

class lagrange_basis(basis):
    def __init__(self, order, domain):
        self.basis = lagrange_basis_function
        basis.__init__(self, order, domain)

class hierarchical_basis(basis):
    def __init__(self, order, domain):
        self.basis = hierarchical_basis_function
        basis.__init__(self, order, domain)

class fv_dual_basis(basis):
    def __init__(self, order, domain):
        self.basis = fv_dual_basis_function
        basis.__init__(self, order, domain)

class dirac_basis(basis):
    def __init__(self, order, domain):
        self.basis = dirac_test_function
        basis.__init__(self, order, domain)

def mass_matrix(P, Q):
    return ImmutableMatrix([[volume_integrate(p * q) for p in P] for q in Q])

def stiffness_matrix(P, Q, T=Identity(2)):
    A_xx = ImmutableMatrix([[volume_integrate(dx(p) * dx(q)) for p in P] for q in Q])
    A_xy = ImmutableMatrix([[volume_integrate(dx(p) * dy(q)) for p in P] for q in Q])
    A_yx = ImmutableMatrix([[volume_integrate(dy(p) * dx(q)) for p in P] for q in Q])
    A_yy = ImmutableMatrix([[volume_integrate(dy(p) * dy(q)) for p in P] for q in Q])

    return T[0, 0] * A_xx + T[0, 1] * A_xy + T[1, 0] * A_yx + T[1, 1] * A_yy

def stiffness_matrix_masked(P, Q, M):
    A = [ImmutableMatrix([[volume_integrate(m * inner(grad(p), grad(q))) for p in P] for q in Q]) for m in M]

    return A

def deriv_matrices(P, Q, T=Identity(2)):
    Ax = ImmutableMatrix([[volume_integrate(dx(p) * q) for p in P] for q in Q])
    Ay = ImmutableMatrix([[volume_integrate(dy(p) * q) for p in P] for q in Q])

    return T[0, 0] * Ax + T[0, 1] * Ay, T[1, 0] * Ax + T[1, 1] * Ay

def boundary_matrices(P, Q):
    one_half = Rational(1, 2)

    domain = [  [(Point(1, 0), Point(one_half, one_half)), (Point(one_half, 0), Point(1, 0))],
                [(Point(0, one_half), Point(0, 0)), (Point(0, 0), Point(one_half, 0))],
                [(Point(one_half, one_half), Point(0, 1)), (Point(0, 1), Point(0, one_half))]
             ]

    domain = [  [(Point(one_half, one_half), Point(one_half, 0))],
                [(Point(one_half, 0), Point(one_half, one_half)), (Point(one_half, one_half), Point(0, one_half))],
                [(Point(0, one_half), Point(one_half, one_half))]
             ]

    Ax = ImmutableMatrix([[boundary_integrate(basis_function((p.expr * q.expr, 0), d)) for p in P] for (q, d) in zip(Q, domain)])
    Ay = ImmutableMatrix([[boundary_integrate(basis_function((0, p.expr * q.expr), d)) for p in P] for (q, d) in zip(Q, domain)])

    return Ax, Ay

def main():
    T1 = Triangle((1, 0), (0, 0), (0, 1))
    T2 = Triangle((1, 0), (Rational(1,2), Rational(1,2)), (0, 0))

    K = MatrixSymbol('K', 2, 2)
    p = lagrange_basis(1, T1)
    q = lagrange_basis(1, T1)
    s = fv_dual_basis(1, T1)

    M = mass_matrix(p, q)
    pprint(Eq(Symbol('M'), M))
    print fcode(M)

    Dx, Dy = deriv_matrices(p, q, K)
    pprint(Eq(Symbol('Dx'), Dx))
    pprint(Eq(Symbol('Dy'), Dy))
    print fcode(Dx)
    print fcode(Dy)

    A = stiffness_matrix(p, q, K)
    pprint(Eq(Symbol('A'), A))
    print fcode(A)

    S = stiffness_matrix_masked(p, q, s)
    pprint(S)
    print fcode(S)
main()

