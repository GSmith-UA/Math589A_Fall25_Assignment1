import cubic_solver
import math
import cmath

def solve_quartic(a, b, c, d, e):
    """Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    Returns a list of 1..4 roots (real numbers or complex numbers).
    If the leading coefficients are zero the function will
    handle lower-degree polynomials automatically.
    """
    roots = []

    # Going to try and handle some edge cases first...
    if e == 0:
        # We can factor an x and turn into a cubic...
        roots.append(0)
        cubicRoots = cubic_solver.solve_cubic(a,b,c,d)
        roots.append(cubicRoots[0])
        roots.append(cubicRoots[1])
        roots.append(cubicRoots[2])
        return roots
    
    if a == 0:
        # This is not a quartic...
        # I will solve it for you anyway...
        roots = cubic_solver.solve_cubic(b,c,d,e)
        return roots
    
    # Ok... lets get a leading coeff of 1 and depress the quartic...
    c4 = 1
    c3 = b/a
    c2 = c/a
    c1 = d/a
    c0 = e/a

    # Going to put it in the form t^4 +qt^2 +rt +s...
    q = c2 - (3/8)*(c3**2)
    r = c1 - (1/2)*c3*c2 + (1/8)*(c3**3)
    s = c0 - (1/4)*c3*c1 + (1/16)*(c3**2)*c2 - (3/256)*(c3**4)

    # If r is zero then this is actually a quadratic problem...
    if r ==0:
        quadRoots = cubic_solver.solve_quad(1,q,s)
        r1 = quadRoots[0]
        r2 = quadRoots[1]
        multVec = [1,1,1,1]
        if r1<0:
            multVec[0] = 0+1j
            multVec[1] = multVec[0]
        if r2<0:
            multVec[2] = 0+1j
            multVec[3] = multVec[2]
        root1 = (-1)*(r1*r1)**(0.25) - (c3/4)
        root2 = (r1*r1)**(0.25) - (c3/4)
        root3 = (-1)*(r2*r2)**(0.25) - (c3/4)
        root4 = (r2*r2)**(0.25) - (c3/4)
        roots.append(root1*multVec[0])
        roots.append(root2*multVec[1])
        roots.append(root3*multVec[2])
        roots.append(root4*multVec[3])
        return roots
    
    # Now we have to construct the resolvent cubic...
    resolveCoeff = [8,8*q,(2*q*q - 8*s),-1*r*r]
    resolveCubicRoots = cubic_solver.solve_cubic(resolveCoeff[0],resolveCoeff[1],resolveCoeff[2],resolveCoeff[3])
    # One of the roots have to be real....
    for root in resolveCubicRoots:
        if math.abs(root.imag)<(10**(-6)):
            chosenRoot = root
            break
    
    linearTerm = (4*chosenRoot*chosenRoot)**(0.25)
    quad1Coeff = [1,-1*linearTerm,(q/2) + chosenRoot + (r/(2*linearTerm))]
    quad2Coeff = [1,linearTerm,(q/2) + chosenRoot - (r/(2*linearTerm))]

    qroot1 = cubic_solver.solve_quad(quad1Coeff[0],quad1Coeff[1],quad1Coeff[2])
    qroot2 = cubic_solver.solve_quad(quad2Coeff[0],quad2Coeff[1],quad2Coeff[2])

    roots.append(qroot1[0] - (c3/4))
    roots.append(qroot1[1] - (c3/4))
    roots.append(qroot2[0] - (c3/4))
    roots.append(qroot2[1] - (c3/4))
    return roots



def main():
    tests = [
        (1, 0, 0, 0, -1),  # roots of x^4 - 1 = 0 (2 real and 2 complex roots)
        (1, 0, 1, 0, -1),  # roots of x^4 + x^2 - 1 = 0 (2 real and 2 complex roots)
        (0, 1, -3, 2, 0),  # a=0 => cubic: x^3 - 3x^2 + 2x = 0  (roots 0,1,2)
    ]
    for a, b, c, d, e in tests:
        roots = solve_quartic(a, b, c, d, e)
        print(f"solve_quartic({a}, {b}, {c}, {d}, {e}) -> {roots}")

if __name__ == "__main__":
    main()
