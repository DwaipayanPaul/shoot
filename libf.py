# Library Module
import math
import csv
import random

############## ASSIGNMENT 7 ###################
# Euler's method
def eulerm(x0,y0,st,h,f,N):
    x,y=[],[]     # to save the data for x and y(x)
    # initial values
    y.append(y0)
    x.append(x0)
    # N=1000
    for i in range(1,N):
        z=y[i-1]+(h*f(y[i-1],x[i-1]))  # y(x_n + h) = y(x_n) + h f(y(x_n), x_n)
        y.append(z)    # storing y(x_n + h)
        x.append(x0+(i*h))  # and x
    write_csv(st, y, x)     # writing the arrays in a csv file

def runge_kutta(f1,f2,x0,y0,z0,h,N,st):    # Runge-Kutta 4
    ax,ay=[],[]   # arrays to store x and y(x) respectively
    # stores the initial value
    ax.append(x0)
    ay.append(y0)
    z=z0

    for i in range(N):   # values of k1,k2,k3,k4 for f1 and f2
        #k1
        k1=f1(z,ay[i],ax[i])
        k1t=f2(z,ay[i],ax[i])
        #k2
        k2=f1(z+(k1t*h)/2, ay[i]+(k1*h)/2, ax[i]+h/2)
        k2t=f2(z+(k1t*h)/2, ay[i]+(k1*h)/2, ax[i]+h/2)
        #k3
        k3=f1(z+(k2t*h)/2, ay[i]+(k2*h)/2, ax[i]+h/2)
        k3t=f2(z+(k2t*h)/2, ay[i]+(k2*h)/2, ax[i]+h/2)
        #k4
        k4=f1(z+(k3t*h), ay[i]+(k3*h), ax[i]+h)
        k4t=f2(z+(k3t*h), ay[i]+(k3*h), ax[i]+h)
        # z and y
        z=z+((k1t+2*k2t+2*k3t+k4t)*h)/6
        yp=ay[i]+((k1+2*k2+2*k3+k4)*h)/6
        # storing x and y in ax and ay
        ay.append(yp)
        xp=x0+(h*(i+1))
        ax.append(xp)

    write_csv(st, ay, ax)  # writing x and y(x) in a csv file

    return(yp)      # returning y(x=xn)

def shooting_method(f1,f2,x0,y0,yn,h,N,st,g):

    y = runge_kutta(f1, f2, x0, y0, g, h, N, "test.csv")    # yn for the first guess
    print(f"Value of y(x=xn) for the above guess {g}=",y)

    # initialising
    lower,upper=0.0,0.0
    # checking if y overshoots or undershoots
    # if y overshoots
    if y > yn and abs(y - yn) > 10E-4:
        upper = g
        # we got upper bracket of y
        # to find lower bound
        while y > yn:
            g = float(input(f"Guess a value of y\'({x0}) lower than the previous guess\n"))
            y = runge_kutta(f1, f2, x0, y0, g, h, N, "test.csv")
            print(f"Value of y(x=xn) for the above guess {g}=", y)

        if abs(y - yn) < 10E-4:                # if yn for the guess is equal to or very near to actual yn
            runge_kutta(f1, f2, x0, y0, g, h, N, st)     # writing the final data file st
            print("Value of y(x=xn) found, integration successful")
            return 0
        else:                                # if yn of guess is less than actual yn
            lower = g                        # then we have found the lower bracket
            lagrange_interpolation(upper,lower,f1,f2,x0,y0,yn,h,N,st)

        # if y undershoots
    elif y < yn and abs(y - yn) > 10E-4:
        lower = g     # got the lower bracket
        # now to find upper bound
        while y < yn:
            g = float(input(f"Guess a value of y\'({x0}) greater than the previous guess\n"))
            y = runge_kutta(f1, f2, x0, y0, g, h, N, "test.csv")
            print(f"Value of y(x=xn) for the above guess {g}=", y)

        if abs(y - yn) < 10E-4:           # if yn for the guess is equal to or very near to actual yn
            runge_kutta(f1, f2, x0, y0, g, h, N, st)      # writing the final data file st
            print("Value of y(x=xn) found, integration successful")
        else:
            upper = g
            lagrange_interpolation(upper, lower, f1, f2, x0, y0, yn, h, N, st)
        #
    elif abs(y - yn) < 10E-4:           # if yn for the guess is equal to or very near to actual yn
        runge_kutta(f1, f2, x0, y0, g, h, N, st)  # if guess gives perfect value of yn at xn, then solution is obtained
        print("Value of y(x=xn) found, integration successful")



def lagrange_interpolation(upper,lower,f1, f2, x0, y0, yn, h, N, st):
    yl = runge_kutta(f1, f2, x0, y0, lower, h, N, st)    # yn for lower bracket

    yh = runge_kutta(f1, f2, x0, y0, upper, h, N, st)    # yn for upper bracket
    # for next y'(x0)
    g = lower + ((upper - lower) / (yh - yl)) * (yn - yl)

    y = runge_kutta(f1, f2, x0, y0, g, h, N, st)   # yn for the new y'(x0)
    print("Value of y(x=xn) found, integration successful")


###################### Assignment 6  #########################################################
# Methods for numerical integration
def midpoint(a,b,n,f):     # Midpoint/Rectangle method
    # width of N equal parts
    h=(b-a)/n
    # x stores values of x1, x2, x3........
    x=0.0
    sum =0.0
    for i in range(n): # values of x
        x=((a+i*h)+(a+((i+1)*h)))/2
        sum += h * f(x)    # summing the values of h*f(x1),h*f(x2),.....

    return sum

def trapezoidal(a,b,n,f):    # Trapezoidal method
    h = (b - a) / n
    # x stores values of x0, x1, x2, x3........ xN
    x = 0
    sum = 0
    w = 1  # intialising weight function w=1 or 2
    for i in range(n+1):
        x = a + (i * h)   # values of x
        if i==0 or i==n:   # w(x0)=w(xN)=1
            w=1
        else: w=2          # w(x1)=w(x2)=.....w(x{N-1})=2
        sum += h * w * f(x)/2  # summing values of h * w * f(xi)/2

    return(sum)

def simpsons(a,b,n,f):    # simpsons method
    h = (b - a) / n

    x = [0 for i in range(n + 1)]   # x stores x0, x1, x2.....
    for i in range(0,n+1,2):
        x[i] = a + (i * h)      # putting the values of x0, x2, x4,....x2n

    for i in range(1,n,2):      # putting the values of x1,x3,x5....
        x[i]=(x[i-1]+x[i+1])/2  # x1 is avg of x0 and x2, x3 is avg of x2 and x4 and so on..
    sum = 0
    w = 1  # intialising weight function w=1,2 or 4

    for i in range(n + 1):
        if i == 0 or i == n:  # w(x0)=w(xN)=1
            w = 1
        elif i%2==0:         # w(xi)=2 for even i
            w = 2
        else: w=4             # w(xi)=4 for odd i
        sum += h * w * f(x[i]) / 3    # summing values of h * w * f(xi)/3

    return sum

def monte_carlo(a,b,n,f):     # monte carlo method
    X=[]   # to store the random variables Xi

    for i in range(n):
        r=random.random()      # random number generator from [0,1]
        r1=a+((b-a)*r)         # r is converted to be in range [a,b]
        X.append(r1)           # storing in X

    # calculation of integral
    sum = 0.0
    for i in range(n):
        sum+=f(X[i])        # f(x1) + f(x2) +.... f(xN)

    p=((b-a)*sum)/n          # value of integral

    return p


def mperr(a, b, e, fd2):        # Error for midpoint method, gives the value of N
    N = float(((((b - a) ** 3) * fd2) / (24 * e)) ** (1 / 2))
    return math.ceil(N)           # returns smallest integer bigger than N


def tperr(a, b, e, fd2):        # Error for trapezoidal method, gives the value of N
    N = float(((((b - a) ** 3) * fd2) / (12 * e)) ** (1 / 2))
    return math.ceil(N)        # returns smallest integer bigger than N


def sperr(a, b, e, fd4):        # Error for simpsons method, gives the value of N
    N = float(((((b - a) ** 5) * fd4) / (180 * e)) ** (1 / 4))
    if N % 2.0 == 0.0:         # only returns even smallest integer bigger than N
        return math.ceil(N)
    else:
        return math.ceil(N + 1)

# To write an array in a csv file
def write_csv(str,er,er1):
    with open(str, 'w', newline='') as file:  # str: the name of file
        writer = csv.writer(file)
        writer.writerow(["X", "y(x)"])  # The first row of the file
        for i in range(len(er)):
            writer.writerow([er1[i], er[i]])     # No of iteration in first column and absolute error in second column

#################################################################################################################
#
def funcdt(func,x):    # deritive of function
    h = 0.0001
    y = (func(x + h) - func(x-h)) / (2*h)
    return y

def bracket(func,a,b):  # bracketing of the root


    for i in range(12):   # iteration limit: 12

        if func(a)*func(b)<0: # roots are on either sides of root(bracketing done)
            print("Bracketing complete: a=",a,"and b=",b)
            return a,b

        elif func(a)*func(b)>0:  # roots are on same side w.r.t.the root

            if abs(func(a))<abs(func(b)):  # need to shift a
                a-=1.5*(b-a)

            elif abs(func(a))>abs(func(b)): # need to shift b
                b=b+1.5*(b-a)

def bisecting(func,a,b):    # Bisection method
    i=0             # counter for iterations
    error=[]         # stores the error values for each iterations
    print("Bisecting method:")
    while(b-a>1E-6 and i<=100):    # iteration limit: 100 and abs error limit: 1E-6

        error.append(round((b - a), 7))  # error value gets appended in the array

        c=(a+b)/2             # midpoint of a and b
        if func(a)*func(c)<0:   #   if a and c are on either sides of root
            b=c               #  b is shifted to c
        elif func(a)*func(c)>0: # if a and c are on same sides of root
            a=c                # a is shifted to c
        else:
            print("Solution found:",c)   # if c is the root:f(c)=0
            return 0
        i+=1
    # print solution
    print("The solution lies between a=",a,"b=",b)
    return error  # return the error array

def falsi(func,a,b):      # Regular falsi method
    i=0            # counter for iterations
    error=[]       # stores the error values for each iterations
    print()
    print("Regular Falsi method:")

    x1,x2=a,b   # counters for calculating error: C_i+1-C_i
    while (abs(x2 - x1) > 1E-6 and i <= 200):   # iteration limit: 200 and abs error limit: 1E-6

        error.append(round(abs(x2 - x1), 7))     # error value gets appended in the array
        # False root c
        c = b - ((b - a) * func(b) / (func(b) - func(a)))

        if func(a) * func(c) < 0:  #   if a and c are on either sides of root
            b = c               #  b is shifted to c
        elif func(a) * func(c) > 0:   # if a and c are on same sides of root
            a = c              # a is shifted to c
        else:
            print("Solution found:", c)   # if c is the root:f(c)=0
            return 0

        if i%2==0:     # if the iteration no. is even
            x2=c       # C_2n=c
        else:
            x1=c       # else C_2n+1=c
        i += 1
    # print output
    print("The solution lies in range ",x1,"and ",x2)
    return error    # return the error array

def newton(func,x1):      # Newton-Raphson method
    i=0         # counter for iterations
    error=[]     # stores the error values for each iterations
    x2=0      # x1 and x2 are counters for calculating error: X_i+1-X_i
    print()
    print("Newton-Raphson method:")

    while(abs(x2-x1)>1E-6 and i<=200):    # iteration limit: 200 and abs error limit: 1E-6

        error.append(round(abs(x2 - x1), 7))    # error value gets appended in the array

        if i%2==0:       # if the iteration no. is even
            x2=x1-(func(x1)/funcdt(func,x1))   # X_2n=X_2n+1 -[f(X_2n)/f'(X_2n+1)]
        else:                                 # else
            x1=x2-(func(x2)/funcdt(func,x2))    # X_2n+1=X_2n -[f(X_2n+1)/f'(X_2n)]
        i+=1
    # print the solution
    print("The solution lies in range ", x2, "and", x1)
    return error    # return the error array




# For Q2
def funct(x,a):  # function to calculate f(x) of a polynomial for a value x
                 # a is the array of co-efficients starting with the constant
    n=len(a)
    sum=0.0
    #  a[i] corresponds to the coefficient of x^i
    for i in range(n-1,-1,-1):
        sum+=a[i]*(x**i)   # stores the value of f(x)
    return sum

def functd1(x,a):  # 1st order derivative of a function (in case of polynomial)
    h=0.001
    y=(funct(x+h,a)-funct(x-h,a))/(2*h)    # f'(x)=[f(x+h)-f(x)]/2h
    return y

def functd2(x,a):  # 2nd order derivative of a function (in case of polynomial)
    h = 0.001
    y = (funct(x + h, a) + funct(x - h, a)-2*funct(x,a)) / (h*h)  # f"(x)=[f(x+h)-f(x-h)]/h^2
    return y

def deflate(sol,a):  # deflation
                     # a is the array of co-efficients starting with the constant
    n=len(a)
    q=[0 for i in range(n-1)]  # intialization of q(x)=p(x)/(x-x_o)
    q[n-2]=a[n-1]      # coefficient of x^n in p is x^n-1 in q
    # synthetic division
    for i in range(n-3,-1,-1):  # from x^n-2 in q
        q[i]=a[i+1]+(sol*q[i+1])

    return q   # final q

def solut(a,i): # to find solutions: i is the guess
                # a is the array of co-efficients starting with the constant
    n=len(a)  # n-1 is no. of roots

    if n!=2:  # when f is not of form: f(x)=x-c

        j1,j2=i,0  # counters for error: \alpha_i+1-\alpha_i
        j = i  # takes the guess for the solution
        a1=0   # a1 is for calculation of a=n/[G(+-)math.sqrt((n-1)*(nH-G^2))]
        k=1   # counter for iterations
        if funct(i,a)!=0:  # when i is not the root of f(x)

            while abs(j2-j1)>1E-6 and k<200:  # iteration limit: 200 and abs error limit: 1E-6
                # calculation G and H
                g=functd1(j,a)/funct(j,a)
                h=g**2-(functd2(j,a)/funct(j,a))
                # denominators : d1 and d2
                det1=g+math.sqrt((n-1)*(n*h-g**2))
                det2=g-math.sqrt((n-1)*(n*h-g**2))

                if abs(det1)>abs(det2):  # if absolute value of det1 is max
                    a1=n/det1          # a=n/[G(+)math.sqrt((n-1)*(nH-G^2))]
                else:
                    a1=n/det2          # a=n/[G(-)math.sqrt((n-1)*(nH-G^2))]

                if k%2==0:          # for even no. iteration
                    j1=j2-a1         # \alpha_2n+1=\alpha_2n - a
                    j=j1            # for next iteration: \alpha_2n+1
                            # else
                else:
                    j2=j1-a1        # \alpha_2n=\alpha_2n+1 - a
                    j=j2            # for next iteration: \alpha_2n
                k+=1

        # The iteration ended in even no.
        if k%2==0:
            print(j1)   # \alpha_2n+1 is the nearest solution
            # deflation and saving the new polynomial q as a(j1 is solution)
            a=deflate(j1,a)
        else:          # else
            print(j2)   # \alpha_2n is the nearest solution
            # deflation and saving the new polynomial q as a(j2 is solution)
            a = deflate(j2, a)
        # return the new polynomial array a
        return a

    else:  # when f is of form: f(x)=x-c
        if a[1]*a[0]<0 or a[1]<0:  # if eq is of form: x-c=0 or -x+c=0
            print(a[0]) # print         x=c (solution)
        else:                    # if eq is of form: -x+c=0
            print(-a[0]) # print        x=-c (solution)

        return 0


# THE LIBRARY MODULE ( with functions involving matrices)

def read_write(st):        # reading and writing matrix
    a=[]
    # Reading matrices from the files
    f1 = open(st, 'r')
    for line in f1.readlines():
        a.append([float(x) for x in line.split(',')])  # adding rows
    return a

def part_pivot(a,b):      # partial pivoting
    n=len(a)
    # initialise
    (c,d)=(0,0)
    for k in range(n-1):
        if a[k][k]==0:     # checking if the diagonal element is zero
            for r in range(k+1,n):
                if abs(a[r][k])>abs(a[k][k]):   # swapping
                    for i in range(n):
                        # swapping in matrix b
                        c = b[r][i]
                        b[r][i] = b[k][i]
                        b[k][i] = c
                        # swapping in matrix a
                        d=a[k][i]
                        a[k][i]=a[r][i]
                        a[r][i]=d


def matrix_mult(m,n):    # multiply two matrices
    l=len(m)
    r=[[ 0.0 for i in range(l) ] for j in range(l)]
    for i in range(l):
        for j in range(l):
            for k in range(l):
                r[i][j] = r[i][j] + (m[i][k] * n[k][j])
    return r

def print_mat(a):      # print a matrix
    n=len(a)
    for i in range(n):
        for j in range(n):
            print(a[i][j]," ",end="")
        print()


