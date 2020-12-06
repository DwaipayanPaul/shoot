import libf
import math

f1=lambda z,y,x:z
f2=lambda z,y,x:z+1
e = 2.71828
x0,y0=0,1
xn,yn=1,2*(e-1)
for i in range(4):
    h=float(input("Enter the value of step size h: \n"))
    print("For h=",h)
    print("-----------")
    n = int((xn - x0) / h)
    g = float(input(f"Guess a value of y\'({x0}) \n"))
    libf.shooting_method(f1,f2,x0,y0,yn,h,n,f"q3data(h={h}).csv",g)

# Data written in csv file
# OUTPUT:-
'''
Enter the value of step size h: 
0.02
For h= 0.02
-----------
Guess a value of y'(0) 
0.4
Value of y(x=xn) for the above guess 0.4= 2.4055945548523856
Guess a value of y'(0) greater than the previous guess
0.8
Value of y(x=xn) for the above guess 0.8= 3.09290728481021
Guess a value of y'(0) greater than the previous guess
1.3
Value of y(x=xn) for the above guess 1.3= 3.952048197257492
Value of y(x=xn) found, integration successful 
Enter the value of step size h: 
0.08
For h= 0.08
-----------
Guess a value of y'(0) 
1.1
Value of y(x=xn) for the above guess 1.1= 3.424560912805738
Guess a value of y'(0) greater than the previous guess
1.3
Value of y(x=xn) for the above guess 1.3= 3.7469000473586647
Value of y(x=xn) found, integration successful 
Enter the value of step size h: 
0.1
For h= 0.1
-----------
Guess a value of y'(0) 
1.3
Value of y(x=xn) for the above guess 1.3= 3.9520434115108816
Guess a value of y'(0) lower than the previous guess
1.1
Value of y(x=xn) for the above guess 1.1= 3.6083874626838486
Guess a value of y'(0) lower than the previous guess
0.9
Value of y(x=xn) for the above guess 0.9= 3.2647315138568147
Value of y(x=xn) found, integration successful 
Enter the value of step size h: 
0.3
For h= 0.3
-----------
Guess a value of y'(0) 
1
Value of y(x=xn) for the above guess 1.0= 3.018973276382043
Guess a value of y'(0) greater than the previous guess
1.3
Value of y(x=xn) for the above guess 1.3= 3.4568192678393492
Value of y(x=xn) found, integration successful
'''