import math, cmath

def solve_cubic(a, b, c, d):
    """Solve a*x^3 + b*x^2 + c*x + d = 0
    Returns list of 1..3 roots (complex if needed).
    """
    roots = []
    # Going to try and handle some degenerate cases...
    if d == 0:
        # We can factor an x out and reduce to a quadratic...
        roots.append(0)
        quadRoots = solve_quad(a,b,c)
        roots.append(quadRoots[0])
        roots.append(quadRoots[1])
        return roots
    
    if a == 0:
        # This is not even a real cubic...
        # This is a quadratic problem...
        # I will solve it for you though...
        print('User has entered a degenerate cubic...')
        roots = solve_quad(b,c,d)
        return roots

    if (b == 0) & (c == 0):
        # This is really a roots of unity problem...
        c3 = 1
        c0 = d/a
        if c0<0:
            crho = (abs(c0))**(1/3)
            phase0 = 0
            phase1 = (2/3)*math.pi
            phase2 = (4/3)*math.pi
            j = 0+1j
            roots.append(crho*cmath.exp(j*phase0))
            roots.append(crho*cmath.exp(j*phase1))
            roots.append(crho*cmath.exp(j*phase2))
            return roots
        else:
            phi = (-1/3)*cmath.pi
            crho = (abs(c0))**(1/3)
            phase0 = phi
            phase1 = phi + (2/3)*cmath.pi
            phase2 = phi + (4/3)*cmath.pi
            j = 0+1j
            roots.append(crho*cmath.exp(j*phase0))
            roots.append(crho*cmath.exp(j*phase1))
            roots.append(crho*cmath.exp(j*phase2))
            return roots
    
    # Ok.... so we have a somewhat interesting cubic if we reach this point...
    # I am going to divide through by leading coeff...
    c3 = 1
    c2 = b/a
    c1 = c/a
    c0 = d/a
    # Now we depress it into t^3 + pt + q
    p = (-1*(c2)**2)/3 + c1
    q = (2*(c3**3)/27) -(c1*c2/3) + c0

    # Here is a helper term...
    if p<0:
        a = (-2/3)*p*(p*p/9)**(0.25)
        threeTheta = cmath.acos(-q/a)
        theta1 = threeTheta/3
        t1 = cmath.cos(theta1)*2*(p*p/9)**(0.25)
        theta2 = theta1 + (2/3)*cmath.pi
        t2 = cmath.cos(theta2)*2*(p*p/9)**(0.25)
        theta3 = theta1 - (2/3)*cmath.pi
        t3 = cmath.cos(theta3)*2*(p*p/9)**(0.25)
    else:
        a = (2/3)*p*(p*p/9)**(0.25)
        threeTheta = cmath.asinh(-q/a)
        theta1 = threeTheta/3
        t1 = cmath.sinh(theta1)*2*(p*p/9)**(0.25)
        theta2 = theta1 + (2/3)*cmath.pi*(0+1j)
        t2 = cmath.sinh(theta2)*2*(p*p/9)**(0.25)
        theta3 = theta1 - (2/3)*cmath.pi*(0+1j)
        t3 = cmath.sinh(theta3)*2*(p*p/9)**(0.25)
    
    roots.append(t1-c2/3)
    roots.append(t2-c2/3)
    roots.append(t3-c2/3)
    return roots
    



def solve_linear(a,b):
    roots = [-b/a]
    return roots

def solve_quad(a,b,c):
    roots = []
    if a == 0:
        # This is just a linear problem
        roots = solve_linear(b,c)
        return roots
    if c == 0:
        # We can factor out a zero...
        roots.append(0)
        roots.append(solve_linear(a,b))
        return roots
    # Now we should have a non-trivial quadratic...
    if a != 1:
        # We want to recondition...
        c0 = c/a
        c1 = b/a
        c2 = 1
    else:
        c0 = c
        c1 = b
        c2 = a

    discrim = c1**2 - 4*c0

    if discrim == 0:
        roots.append(-1*c1/2)
        roots.append(-1*c1/2)
        return roots
    
    if (discrim < 8)&( discrim >0):
        # We should be able to hit it with arcos
        alpha1 = math.acos((c1**2 -4*c0-4)/4)
        alpha2 = 2*math.pi - alpha1

        alpha1/=2
        alpha2/=2

        y1 = math.cos(alpha1)
        y2 = math.cos(alpha2)

        t1 = 4**(0.25)*y1
        t2 = 4**(0.25)*y2

        roots.append(t1 - (c1/2))
        roots.append(t2 - (c1/2))
        return roots
    else:
        alpha1 = math.acosh((c1**2 -4*c0-4)/4)
        alpha1/=2
        y1 = math.cosh(alpha1)
        t1 = 4**(0.25)*y1
        roots.append(t1-(c1/2))
        roots.append((c0/(roots[0])))
        return roots

        



def main():
    tests = [
        (1, 0, 0, -1),     # roots of x^3 - 1 = 0 (1 and two complex cube roots)
        (1, -6, 11, -6),   # roots [1.0, 2.0, 3.0]
    ]
    for a, b, c, d in tests:
        roots = solve_cubic(a, b, c, d)
        print(f"solve_cubic({a}, {b}, {c}, {d}) -> {roots}")

if __name__ == "__main__":
    main()
