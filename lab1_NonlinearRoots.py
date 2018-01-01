EPSILON = 0.000001
DX = 0.0000000001


def f(x):
    return 4*x**5 - 3*x**4 + x ** 3 + 2*x**2 - 4*x + 3


def derivative(x):
    return (f(x + DX) - f(x - DX)) / (2 * DX)


class Bisection:
    iteration = 0

    def process(self, a, b):
        self.iteration += 1
        value = (b + a) / 2
        print("iter " + repr(self.iteration).rjust(2) + ":   a = %f; b = %.10f; f(a) = %.10f; f(b) = %.10f; f(value) = %.10f" % (a, b, f(a), f(b), f(value)))

        if abs(b - a) < EPSILON and abs(f(value)) < EPSILON:
            return value

        elif f(value) * f(a) < 0:
            return self.process(a, value)

        else:
            return self.process(value, b)


class Chord:
    iteration = 0

    def process(self, a, b, value0):
        self.iteration += 1
        value = (a * f(b) - b * f(a)) / (f(b) - f(a))
        print("iter " + repr(self.iteration).rjust(2) + ":   a = %.10f; b = %.10f; f(a) = %.10f; f(b) = %.10f; f(value) = %.10f" % (a, b, f(a), f(b), f(value)))

        if abs(value0 - value) < EPSILON and abs(f(value)) < EPSILON:
            return value

        elif f(value) * f(a) < 0:
            return self.process(a, value, value)

        else:
            return self.process(value, b, value)


class Newton:
    iteration = 0

    def process(self, value0):
        self.iteration += 1
        value = value0 - f(value0) / derivative(value0)
        print("iter " + repr(self.iteration).rjust(2) + ":   value = %.10f; f(value) = %.10f;" % (value, f(value)))

        if abs(value - value0) < EPSILON and abs(f(value)) < EPSILON:
            return value

        else:
            return self.process(value)

bisection = Bisection()
chord = Chord()
newton = Newton()

try:
    print("\n\tPROCESSING BISECTION METHOD\n")
    print("Result: expectable root value = %.8s \n" % repr(bisection.process(-2, -1)).rjust(5))
except OverflowError:
    print("BISECTION DIVERGES")

try:
    print("\n\tPROCESSING CHORD METHOD\n")
    print("Result: expectable root value = %.8s \n" % repr(chord.process(-2, -1, -1.5)).rjust(5))
except OverflowError:
    print("CHORD DIVERGES")

try:
    print("\n\tPROCESSING NEWTON METHOD\n")
    print("Result: expectable root value = %.8s \n" % repr(newton.process(-1.5)).rjust(5))
except OverflowError:
    print("NEWTON DIVERGES")