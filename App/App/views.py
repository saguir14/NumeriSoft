from django.shortcuts import render
from django.http import HttpResponse
#from .forms import CreateNewList
import numpy as np
from numpy import array, zeros, diag, diagflat, dot
from numpy import linalg as LA
import sys
import scipy
import scipy.linalg
from prettytable import PrettyTable
from py_expression_eval import Parser
import ast
import sympy as sp
import sympy as sym

def index(request):
   return render(request, "index.html", {})

def IncrementalSearch(request):
   return render(request, "IncrementalSearchRes.html", {}) 

def IncrementalSearchRes(request):
   #x0=float(input("Ingrese valor inicial: "))
   x0=request.POST.get('x0', float)
   #delta=float(input("Ingrese valor delta de x: "))
   delta = request.POST.get('delt', float)
   #niter=float(input("Ingrese número de iteraciones: "))
   niter = request.POST.get('iter', int)
   #func=str(input("Ingrese la función: "))
   func = request.POST.get('func', str)
   x0 = float(x0)
   delta = float(delta)
   niter = int(niter)
   #t2 = PrettyTable(['iterations','x0','f(x)', 'message'])
   #t2.border = False
   t = PrettyTable(['iterations','x0','f(x)', 'message'])
   t.border = False
   t.header = False
   parser = Parser()
   f = parser.parse(func).evaluate({'x':x0})
   fmelo = "{:f}".format(f)
   if f==0:
      print("EL valor inicial ingresado es una raiz")
   else:
      x1=x0+delta
      cont=0
      fx1=parser.parse(func).evaluate({'x':x1})
      print("| Iteración | x | f(x) |")
      print(" | ",cont," | ",x0," | ",fmelo, " | ")
      #context = {'iter':cont, 'x0':x0, 'fx':fmelo}
      t.add_row([" | "+str(cont)+" | ",str(x0)+" | ",str(fmelo)+" | ", ' '])
      #tmelo = str(t).replace('|',' ').replace('+',' ')
      while cont < niter:
         x0=x1
         f=fx1
         x1=x0+delta
         fx1=parser.parse(func).evaluate({'x':x1})
         fmelo = "{:f}".format(f)
         cont=cont+1
         #context = {'iter':cont, 'x0':x0, 'fx':fmelo}
         if fx1==0:
            print("x1 is a root")
         else:
            if f*fx1 < 0:
               print(" | ",cont," | ",x0," | ",fmelo," | "," There is a root between",x0," ",x1)
               t.add_row([" | "+str(cont)+" | ", str(x0)+" | ", str(fmelo)+" | ", 'There is a root between'+' '+str(x0)+' and '+str(x1) ])
               #tmelo = str(t).replace('|',' ').replace('+',' ')
               #context = {'iter':cont, 'x0':x0, 'fx':fmelo}
            else:
               print(" | ",cont," | ",x0," | ",fmelo," | ")
               t.add_row([" | "+str(cont)+" | ", str(x0)+" | ", str(fmelo)+" | ", ' '])
               #tmelo = str(t).replace('|',' ').replace('+',' ')
               #context = {'iter':cont, 'x0':x0, 'fx':fmelo}
   print(t)
   #k = t.from_html_one
   context = {'vec':t}
   #context2 = {'enca':t2}
   return render(request,"IncrementalSearchRes.html",context)

def Biseccion(request):
   return render(request, "biseccionres.html", {}) 

def BiseccionRes(request):
   #xi=float(input("Ingrese límite inferior: "))
   xi = request.POST.get('xi', None)
   #xs=float(input("Ingrese límite superior: "))
   xs = request.POST.get('xs', None)
   #niter=float(input("Ingrese número de iteraciones: "))
   tol = request.POST.get('tol', None)
   niter = request.POST.get('iter', None)
   #tol=float(input("Ingrese tolerancia: "))
   #func=str(input("Ingrese la función: "))
   func = request.POST.get('func', None)
   uno=float(1)

   xi = float(xi)
   xs = float(xs)
   niter = int(niter)
   tol = float(tol)
   t = PrettyTable(['iterations','xi','xm','xs','f(x)','Error'])
   t.border = False
   t.header = False
   parser = Parser()
   fxi = parser.parse(func).evaluate({'x':xi})
   fxs = parser.parse(func).evaluate({'x':xs})
   xm = (xs+xi)/2
   fxm = parser.parse(func).evaluate({'x':xm})
   ximelo = "{:f}".format(xi)
   xsmelo = "{:f}".format(xs)
   xmmelo = "{:f}".format(xm)
   fxmmelo = "{:e}".format(fxm)
   print(" | "," Iterations "," | ", " xi ", " | ", " xm ", " | ", " xs ", " | ", " fx ", " | ", " error (E) ", " | ")
   print(" | ", "1"," | ", ximelo, " | ", xmmelo, " | ", xsmelo, " | ", fxmmelo, " |      ", "      | ")
   t.add_row([" | "+str( 1 )+" | ",str(ximelo)+" | ",str(xmmelo)+" | ",str(xsmelo)+" | ",str(fxmmelo)+" | ","       | "])
   if fxi==0:
      print("The value {} is a root".format(xi))
   else:
      if fxs==0:
         print("The value {} is a root".format(xs))
      else:
         if (fxi*fxs) < 0:
            xm=(xi+xs)/2
            fxm = parser.parse(func).evaluate({'x':xm})
            cont=1
            error=(tol)+uno
            while error > tol and fxm != 0 and cont < niter:
               if (fxi*fxm) < 0:
                  xs=xm
                  fxs=fxm
               else:
                  xi=xm
                  fxi=fxm
               xaux=xm
               xm=(xi+xs)/2
               fxm = parser.parse(func).evaluate({'x':xm})
               error=abs(xm-xaux)
               cont=cont+1
               ximelo = "{:f}".format(xi)
               xsmelo = "{:f}".format(xs)
               xmmelo = "{:f}".format(xm)
               fxmmelo = "{:e}".format(fxm)
               errmelo = "{:e}".format(error)
               print(" | ", cont," | ", ximelo, " | ", xmmelo, " | ", xsmelo, " | ", fxmmelo, " | ", errmelo, " | ")
               t.add_row([" | "+str(cont)+" | ",str(ximelo)+" | ",str(xmmelo)+" | ",str(xsmelo)+" | ",str(fxmmelo)+" | ",str(errmelo)+" | ",])

            if fxm==0:
               print("The value {} is a root".format(xmmelo))
            else:
               if error < tol:
                  print("{} is an aproximation of a root with {} tolerance".format(xmmelo,errmelo))
               else:
                  print("The method fails in {} iterations".format(niter))
         else:
            print("The interval is inappropriate")
   context2 = {'vec':t}
   return render(request, "biseccionres.html", context2)

def fixedpoint(request):
   return render(request, "fixedpointres.html", {})

def fixedpointres(request):
   import math
   from py_expression_eval import Parser

   #tol=float(input("Ingrese tolerancia: "))
   #xa=float(input("Ingrese valor inicial: "))
   #niter=float(input("Ingrese número de iteraciones: "))
   #func=str(input("Ingrese la función f: "))
   #fung=str(input("Ingrese función g: "))
   xa=request.POST.get('x0', None)
   niter=request.POST.get('iter', None)
   tol=request.POST.get('tol', None)
   func=request.POST.get('func', None)
   fung=request.POST.get('fung', None)

   xa=float(xa)
   niter = int(niter)
   tol = float(tol)
   t = PrettyTable(['Iterations','xi','g(xi)','f(xi)','Error'])
   t.border = False
   t.header = False
   parser = Parser()
   fx=parser.parse(func).evaluate({'x':xa})
   cont=0
   error=tol+1
   print("| Iter | xi | g(xi) | f(xi) | E |")
   while fx != 0 and error > tol and cont < niter:
      xn=parser.parse(fung).evaluate({'x':xa})
      fx=parser.parse(func).evaluate({'x':xn})
      error=abs(xn-xa)
      xa=xn
      cont=cont+1
      print(" | ",cont," | ",xa," | ",xn, " | ",fx," | ",error," | ")
      t.add_row([" | "+str(cont)+" | ",str(xa)+" | ",str(xn)+" | ", str(fx)+" | ",str(error)+" | "])
   if fx==0:
      print("{} es una raiz".format(xa))
   else:
      if error < tol:
         print("{} es una aproximación con una tolerancia de {}".format(xa,error))
      else:
         print("El método fracasó en {} iteraciones".format(niter))
   context6={'vec':t}
   return render(request, "fixedpointres.html", context6)

def falseposition(request):
   return render(request, "falsepositionres.html", {})

def falsepositionres(request):
   import math
   from py_expression_eval import Parser

   #xi=float(input("Ingrese límite inferior: "))
   #xs=float(input("Ingrese límite superior: "))
   #niter=float(input("Ingrese número de iteraciones: "))
   #tol=float(input("Ingrese tolerancia: "))
   #func=str(input("Ingrese la función: "))
   xi=request.POST.get('xi', None)
   xs=request.POST.get('xs', None)
   niter=request.POST.get('iter', None)
   tol=request.POST.get('tol', None)
   func=request.POST.get('func', None)
   
   xi=float(xi)
   xs=float(xs)
   niter = int(niter)
   tol = float(tol)
   t = PrettyTable(['Iterations','xi','xm','xs','f(xm)','Error'])
   t.border = False
   t.header = False
   parser = Parser()
   #fmelo = "{:f}".format(f)

   yi=parser.parse(func).evaluate({'x':xi})
   ys=parser.parse(func).evaluate({'x':xs})
   if yi==0:
      print("{} es una raíz".format(xs))
   else:
      if (yi*ys) < 0:
         xm=xi-((yi*(xs-xi))/(ys-yi))
         cont=1
         ym=parser.parse(func).evaluate({'x':xm})
         error=tol+1
         print("| Iteración | xi | xm | xs | f(xm) | E |")
         print(" | ",cont," | ",xi," | ",xm, " | ",xs," | ",ym," | ",error," | ")
         t.add_row([" | "+str(cont)+" | ",str(xi)+" | ",str(xm)+" | ", str(xs)+" | ",str(ym)+" | ",str(error)+" | "])
         while error > tol and ym != 0 and cont < niter:
               if yi*ym < 0:
                  xs=xm
                  ys=ym
               else:
                  xi=xm
                  yi=ym
               xaux=xm
               xm=xi-((yi*(xs-xi))/(ys-yi))
               ym=parser.parse(func).evaluate({'x':xm})
               error=abs(xm-xaux)
               cont=cont+1
               print(" | ",cont," | ",xi," | ",xm, " | ",xs," | ",ym," | ",error," | ")
               t.add_row([" | "+str(cont)+" | ",str(xi)+" | ",str(xm)+" | ", str(xs)+" | ",str(ym)+" | ",str(error)+" | "])
         if ym==0:
               print("{} es una raiz".format(xm))
         else:
            print("The interval is inappropriate")
   context5 = {'vec':t}
   return render(request, "falsepositionres.html", context5)

def newton(request):
   return render(request, "newtonres.html", {})

def newtonres(request):
   import math
   from py_expression_eval import Parser

   x0=request.POST.get('x0', None)
   niter=request.POST.get('iter', None)
   tol=request.POST.get('tol', None)
   func=request.POST.get('func', None)
   fprim=request.POST.get('fprim', None)

   x0 = float(x0)
   niter = int(niter)
   tol = float(tol)
   t = PrettyTable(['iterations','xi','f(x)','Error'])
   t.border = False
   t.header = False
   parser = Parser()
   fx = parser.parse(func).evaluate({'x':x0})
   derif = parser.parse(fprim).evaluate({'x':x0})
   cont=0
   error=tol+1
   x0melo = "{:f}".format(x0)
   fximelo = "{:e}".format(fx)
   print(" | "," Iter "," |   ", " xi ", "   |    ", " f(xi) ", "    |  ", " error (E) ", " | ")
   print(" |    ", cont,"   | ", x0melo, " | ", fximelo, " |     ","          | ")
   t.add_row([" | "+str(cont)+" | ",str(x0melo)+" | ",str(fximelo)+" | ","       | "])
   while error > tol and fx != 0 and derif != 0 and cont < niter:
      x1=x0-(fx/derif)
      fx=parser.parse(func).evaluate({'x':x1})
      derif=parser.parse(fprim).evaluate({'x':x1})
      error=abs(x1-x0)
      x0=x1
      cont=cont+1
      x0melo = "{:f}".format(x0)
      fximelo = "{:e}".format(fx)
      errmelo = "{:e}".format(error)
      print(" |    ", cont,"   | ", x0melo, " | ", fximelo, " | ", errmelo, " | ")
      t.add_row([" | "+str(cont)+" | ",str(x0melo)+" | ",str(fximelo)+" | ",str(errmelo)+" | "])
   if fx==0:
      print("{} es una raiz".format(x0))
   else:
      if error < tol:
         print("{} es una aproximación a una raíz con una tolerancia de {}".format(x1,error))
      else:
         if derif==0:
            print("{} es una posible raíz múltiple".format(x1))
         else:
            print("El método fracasó en {} de iteraciones".format(niter))
   context3 = {'vec':t}
   return render(request, "newtonres.html", context3)

def secant(request):
   return render(request, "secantres.html", {})

def secantres(request):
   import math
   from py_expression_eval import Parser

   #x0=float(input("Ingrese límite inferior: "))
   x0=request.POST.get('x0', None)
   #x1=float(input("Ingrese límite superior: "))
   x1=request.POST.get('x1', None)
   #niter=float(input("Ingrese número de iteraciones: "))
   niter=request.POST.get('iter', None)
   #tol=float(input("Ingrese tolerancia: "))
   tol=request.POST.get('tol', None)
   #func=str(input("Ingrese la función: "))
   func=request.POST.get('func', None)

   x0 = float(x0)
   x1 = float(x1)
   niter = int(niter)
   tol = float(tol)
   t = PrettyTable(['iterations','xi','f(x)','Error'])
   t.border = False
   t.header = False
   parser = Parser()
   fx0=parser.parse(func).evaluate({'x':x0})
   if fx0==0:
      print("{} es raíz".format(x0))
   else:
      fx1=parser.parse(func).evaluate({'x':x1})
      cont=1
      error=tol+1
      x0melo = "{:f}".format(x0)
      fximelo = "{:e}".format(fx0)
      x1melo = "{:f}".format(x1)
      fx2melo = "{:e}".format(fx1)
      print(" | "," Iter "," |   ", " xi ", "   |    ", " f(xi) ", "    |  ", " error (E) ", " | ")
      print(" |    ", "0","   | ", x0melo, " | ", fximelo, " |     ","          | ")
      t.add_row([" | "+str(0)+" | ",str(x0melo)+" | ",str(fximelo)+" | ","       | "])
      print(" |    ", "1","   | ", x1melo, " | ", fx2melo, "  |     ","          | ")
      t.add_row([" | "+str(1)+" | ",str(x0melo)+" | ",str(fximelo)+" | ","       | "])
      while error > tol and fx1 != 0 and (fx1-fx0) != 0 and cont < niter:
         fx0=parser.parse(func).evaluate({'x':x0})
         fx1=parser.parse(func).evaluate({'x':x1})
         x2=x1-fx1*(x1-x0)/(fx1-fx0)
         fx2=parser.parse(func).evaluate({'x':x2})
         error=abs(x2-x1)
         cont=cont+1
         x0melo = "{:f}".format(x2)
         fximelo = "{:e}".format(fx2)
         errmelo = "{:e}".format(error)
         print(" |    ", cont,"   | ", x0melo, " | ", fximelo, "  | ", errmelo, " | ")
         t.add_row([" | "+str(cont)+" | ",str(x0melo)+" | ",str(fximelo)+" | ",str(errmelo)+" | "])
         x0=x1
         x1=x2
      if fx1==0:
         print("{} es una raiz".format(x1))
      else:
         if error < tol:
            print("{} es aproximación a una raiz con una tolerancia de {}".format(x1,error))
         else:
            if (fx1-fx0)==0:
               print("Hay una posible raíz múltiple")
            else:
               print("El método fracasó en {} iteraciones".format(niter))
   context4 = {'vec':t}
   return render(request, "secantres.html", context4)

def multipleroots(request):
   return render(request, "multiplerootsres.html", {})

def multiplerootsres(request):
   import math
   from py_expression_eval import Parser

   x0=request.POST.get('x0', None)
   niter=request.POST.get('iter', None)
   tol=request.POST.get('tol', None)
   func=request.POST.get('func', None)
   df=request.POST.get('df', None)
   d2f=request.POST.get('d2f', None)
   t = PrettyTable(['Iterations','xi','f(xi)','Error'])
   t.border = False
   t.header = False

   x0 = float(x0)
   niter = int(niter)
   tol = float(tol)

   parser = Parser()

   print("| Iter | xi | f(xi) | E |")
   for i in range(1,niter):
      fx=parser.parse(func).evaluate({'x':x0})
      f1x=parser.parse(df).evaluate({'x':x0})
      f2x=parser.parse(d2f).evaluate({'x':x0})
      nuevo=x0-((fx*f1x)/((f1x**2)-(fx*f2x)))
      error=abs(nuevo-f1x)
      if error < tol:
         print(x0)
         break 
      x0=nuevo
      print(" | ",i," | ",x0," | ",fx, " | ",error," | ")
      t.add_row([" | "+str(i)+" | ",str(x0)+" | ",str(fx)+" | ",str(error)+" | "])
   context7 = {'vec':t}
   return render(request, "multiplerootsres.html", context7)

#Views Metodos Lineales xD
def CroutRender(request):
   return render(request, "Crout.html", {})


#Codigos Izi
def Crout(request):
   Matrix = request.POST.get('matriz', None)
   b = request.POST.get('vector', None)
   b = np.matrix(b)
   b = np.array(b)

   Matrix = np.matrix(Matrix)
   Matrix = np.array(Matrix)
   n, m = Matrix.shape
   U = np.ones((n,m))
   L = np.zeros((n,m))

   for k in range(0, n):
       suma=0
       for p in range(0,k):
           suma = suma + L[k][p]*U[p][k]
       L[k][k] = Matrix[k][k]-suma
       for i in range(k+1, n):
           suma2 = 0
           for p in range(0, k):
               suma2 = suma2+L[i][p]*U[p][k]
           L[i][k]= (Matrix[i][k]-suma2)/U[k][k]
       for j in range(k+1, n):
           suma3 = 0
           for p in range(0, k):
               suma3 = suma3+L[k][p]*U[p][j]
           U[k][j]=(Matrix[k][j]-suma3)/L[k][k]
   UF = np.triu(U)

   z = ProSubstitution(L, b)
   z1 = np.vstack(z)
   x = RegSubstitution(UF, z1)

   UF = np.array(UF)
   L = np.array(L)

   print(L)
   print(UF)
   context = {'L': L, 'U' : UF, 'Z': z, 'X': x}
   return render(request, "CroutR.html", context)

def CholeskyRender(request):
   return render(request, "Cholesky.html", {})

def Cholesky(request):
   Matrix = request.POST.get('matriz', None)
   b = request.POST.get('vector', None)
   b = np.matrix(b)
   b = np.array(b)

   Matrix = np.matrix(Matrix)
   Matrix = np.array(Matrix)
   if np.all(np.linalg.eigvals(Matrix)) == True:
    print('|----------Cholesky ERROR-------------|')
    print('| A is not a definite positive matrix |')
    print('|-------------------------------------|')
    mensaje1 = '|----------Cholesky ERROR-------------|'
    mensaje2 = '| A is not a definite positive matrix |'
    mensaje3 = '|-------------------------------------|'
    context = {'M1': mensaje1, 'M2': mensaje2, 'M3': mensaje3}  
    return render(request, "CholeskyR.html", context)
   n, m = Matrix.shape
   U = np.ones((n,m))
   L = np.ones((n,m))
   for k in range(0, n):
      suma=0
      for p in range(0,k):
         suma = suma + L[k][p]*U[p][k]
         L[k][k] = (Matrix[k][k]-suma)**(1/2)
         U[k][k] = L[k][k]
         for i in range(k+1, n):
            suma2 = 0
            for p in range(0, k):
               suma2 = suma2+L[i][p]*U[p][k]
            L[i][k]= (Matrix[i][k]-suma2)/U[k][k]
         for j in range(k+1, n):
           suma3 = 0
           for p in range(k):
              suma3 = suma3+L[k][p]*U[p][j]
           U[k][j]=(Matrix[k][j]-suma3)/L[k][k]

   z = ProSubstitution(L, b)
   z1 = np.vstack(z)
   x = RegSubstitution(U, z1)

   U = np.array(U)
   L = np.array(L)   

   context = {'L': L, 'U' : U, 'Z': z, 'X': x}
   return render(request, "CholeskyR.html", context)

def DoolittleRender(request):
   return render(request, "Doolittle.html", {})

def Doolittle(request):
   Matrix = request.POST.get('matriz', None)
   b = request.POST.get('vector', None)
   b = np.matrix(b)
   b = np.array(b)

   Matrix = np.matrix(Matrix)
   Matrix = np.array(Matrix)
   n, m = Matrix.shape
   U = np.zeros((n,m))
   L = np.ones((n,m))
   for k in range(0, n):
      suma=0
      for p in range(0,k):
         suma = suma + L[k][p]*U[p][k]
      U[k][k] = Matrix[k][k]-suma
      for i in range(k+1, n):
         suma2 = 0
         for p in range(0, k):
            suma2 = suma2+L[i][p]*U[p][k]
         L[i][k]= (Matrix[i][k]-suma2)/U[k][k]
      for j in range(k+1, n):
        suma3 = 0
        for p in range(0, k):
           suma3 = suma3+L[k][p]*U[p][j]
        U[k][j]=(Matrix[k][j]-suma3)/L[k][k]
   LF = np.tril(L)

   z = ProSubstitution(LF, b)
   z1 = np.vstack(z)
   x = RegSubstitution(U, z1)

   U = np.array(U)
   LF = np.array(L) 
   context = {'L': LF, 'U' : U, 'Z': z, 'X': x}
   return render(request, "DoolittleR.html", context)

def LUParcialRender(request):
   return render(request, "LUParcial.html", {})

def LUParcial(request):
   Matrix = request.POST.get('matriz', None)
   b = request.POST.get('vector', None)
   b = np.matrix(b)
   b = np.array(b)
   Matrix = np.matrix(Matrix)
   Matrix = np.array(Matrix)


   P, L, U = scipy.linalg.lu(Matrix)
   z = ProSubstitution(L, b)
   z1 = np.vstack(z)
   x = RegSubstitution(U, z1)
   P = np.array(P)
   L = np.array(L)
   U = np.array(U)

   context = {'P': P , 'U' : U, 'L': L, 'Z': z, 'X': x}
   return render(request, "LUParcialR.html", context)

#Codigos inter Cris
def newtoninterp(request):
   return render(request, "newtoninterpres.html",{})

def newtoninterpres(request):
   
   x = request.POST.get('x', None)
   y = request.POST.get('y', None)

   x = np.array(x)
   y = np.array(y).T

   x = x.tolist()
   y = y.tolist()

   x1 = ast.literal_eval(x)
   y1 = ast.literal_eval(y)

   L = np.ma.size(x1)
   n = L-1
   nprev = n+1
   s = (nprev,L)
   dd = np.ma.zeros(s)
   dd[:,0] = y1
   print (dd)

   for k in range(1, (n+1)):
      for j in range(k,(n+1)):
         d1=[(dd[j,(k-1)])-(dd[(j-1),(k-1)])] 
         d2=[(x1[j])-(x1[j-k])]
         dd[j,k]= np.divide(d1,d2)
   xa = x1
   ya = y1
   y2 = len(y1)
   y3 = len(y1)
   z = (y2,y3)
   d = np.ma.zeros(z)
   d[:,0] = y1
   tam = (len(x1))
   tam2 = tam
   for k in range(1, tam):
      for j in range(0, (tam2-k)):
         e1=[(d[(j+1),(k-1)])-(d[j,(k-1)])]
         e2=[(x1[j+k])-(x1[j])]
         d[j,k]=np.divide(e1,e2)

   l = len(x1)
   polact = []
   aux = 0
   sg1 = ""
   sg2 = ""
   ini = "("
   fin = ")"
   for w in range(0,l):
      ds=np.array2string(np.absolute(d[0,w]))
      dsok=''.join(map(str,ds))
      if w>0:
         if x1[w-1]<0:
            sg1="+"
      else:
         sg1="-"
      if d[0,w]<0:
         sg2="-"
      else:
         sg2="+"
      if w==0:
         acum=np.array2string(d[0,0])
         acumok=''.join(map(str,acum))
      elif w==1:
         aux=np.array2string(np.absolute(x1[w-1]))
         auxok=''.join(map(str,aux))
         polact=[ini,"x",sg1,auxok,fin]
         polactok=''.join(map(str,polact))
         actual=[dsok,"*",polactok]
         actualok=''.join(map(str,actual))
         acum=[acum,sg2,actualok]
         acumok=''.join(map(str,acum))
      else:
         aux=np.array2string(np.absolute(x1[w-1]))
         polact=[polact,"*",ini,"x",sg1,aux,fin]
         polactok=''.join(map(str,polact))         
         actual=[dsok,"*",polactok]
         actualok=''.join(map(str,actual))
         acum=[acum,sg2,actualok]
         acumok=''.join(map(str,acum))

   print("The Matrix's divided differences:")
   print(dd)
   print('\n')
   print("Newton's interpolation polynomial: ")
   acumcasi = str(acum).replace(' ','').replace(',','').replace('[','').replace(']','').replace("'",'').replace('"','').replace('.)',')')
   print (acumcasi)
   print('\n')
   context7 = {'vec':acumcasi}
   return render(request, "newtoninterpres.html", context7)

def vandermonde(request):
   return render(request, "vandermonderes.html", {})

def vandermonderes(request):
   x = request.POST.get('x', None)
   y = request.POST.get('y', None)

   x = np.array(x)
   print(x)
   #y = np.array(y).T
   y = np.array(y)

   x = x.tolist()
   y = y.tolist()

   x1 = ast.literal_eval(x)
   y1 = ast.literal_eval(y)

   n = len(x1)-1
   A = np.vander(x1)
   print("The Vandermonde matrix is:")
   print(A)
   Ainv = np.linalg.inv(A)
   resul = np.matmul(Ainv,y1)
   invarray = resul[::-1]
   print('\n')

   m=0
   melo2 =[]
   melo=[]
   for i in range(0, n+1): 
      xi = invarray[m]
      if xi < 0:
         todo = (xi,"*","x^",i)
      else:
         todo = ("+",xi,"*","x^",i)
      m = m+1
      melo2.append(xi)
      melo.append(todo)
   print("The polynomial coefficients are:")
   invmelo2 = melo2[::-1]
   print(invmelo2)
   print('\n')
   print("The interpolant polynomial is:")
   acumcasi = str(melo).replace(' ','').replace(',','').replace('[','').replace(']','').replace("'",'').replace('"','').replace('.)',')').replace('+-','-').replace('(-','-(').replace('(+','+(')
   print(acumcasi)
   context8 = {'vec':acumcasi}
   return render(request, "vandermonderes.html", context8)

def lagrange(request):
   return render(request, "lagrangeres.html", {})

def lagrangeres(request):
   x = request.POST.get('x', None)
   y = request.POST.get('y', None)

   x = np.array(x)
   y = np.array(y)

   x = x.tolist()
   y = y.tolist()

   x1 = ast.literal_eval(x)
   y1 = ast.literal_eval(y)
   n = len(x1)
   x = sym.Symbol('x')
   polinomio = 0
   divisorL = np.zeros(n, dtype = float)
   for i in range(0,n,1):
      # Termino de Lagrange
      numerador = 1
      denominador = 1
      for j  in range(0,n,1):
         if (j!=i):
            numerador = numerador*(x-x1[j])
            denominador = denominador*(x1[i]-x1[j])
      terminoLi = numerador/denominador
      polinomio = polinomio + terminoLi*y1[i]
      divisorL[i] = denominador
   # simplifica el polinomio
   polisimple = polinomio.expand()
   # para evaluación numérica
   px = sym.lambdify(x,polisimple)

   print('    valores de fi: ',y1)
   print('divisores en L(i): ',divisorL)

   print()
   print('Polinomio de Lagrange: ')
   print(polisimple)
   acum = str(polisimple).replace('**','^')
   context9 = {'vec':acum}
   return render(request, "lagrangeres.html", context9)

#Codigos izin't

def PartialPivotingRender(request):
   return render(request, "PartialPivoting.html", {})

def PartialPivoting(request):
   A = request.POST.get('matriz', None)
   A = np.matrix(A)
   A = np.array(A)
   n, m = A.shape
   # Create zero matrices for L and U 
   U = np.zeros((n,m))
   L = np.zeros((n,m))  
   # Create the pivot matrix P and the multipled matrix PA     
   P = np.identity(n)
   PA = np.identity(n)
   for j in range(n):
     # All diagonal entries of L are set to unity                                                                                                                                                                                                   
      L[j][j] = 1.0
      for i in range(j+1):
         s1 = sum(U[k][j] * L[i][k] for k in range(i))
         U[i][j] = PA[i][j] - s1
      for i in range(j, n):
         s2 = sum(U[k][j] * L[i][k] for k in range(j))
         L[i][j] = (PA[i][j] - s2) / U[j][j]

   PA = np.array(P)
   L = np.array(L)
   U = np.array(U)
   context = {'P': PA , 'U' : U, 'L': L}
   return render(request, "PartialPivotingR.html", context)

def JacobiRender(request):
   return render(request, "Jacobi.html", {})

def Jacobi(request):
   A = request.POST.get('matriz', None)
   b = request.POST.get('vectorB', None)
   guess = request.POST.get('vectorGuess', None)
   iter = request.POST.get('iteraciones', None)
   tol = request.POST.get('tolerancia', None)

   tol = float(tol)
   iter = int(iter)
   
   guess = np.matrix(guess)
   guess = np.array(guess)

   b = np.matrix(b)
   b = np.array(b)

   A = np.matrix(A)
   A = np.array(A)

   t = PrettyTable(['iterations','Error','xi'])
   t.border = False
   t.header = False
   norma = LA.norm(A, ord=2, axis=None, keepdims=False)
   inv = LA.linalg.inv(A)
   norminv = LA.norm(inv, ord=2, axis=None, keepdims=False)
   cond = norma*norminv

   deter=LA.linalg.det(A)
   if deter==0:
      print('\n'"The determinant of the matrix is 0, the problem hasn't solution")
      sys.exit()
   
   n = len(b)
   d = np.diag(np.diag(A))
   l = d - np.tril(A)
   u = d - np.triu(A)
   print("The transition matrix is:")
   dinv = LA.linalg.inv(d)
   lu=l+u
   T = np.matmul(dinv,lu)

   vecprops = LA.eigvals(T)
   absvec = np.absolute(vecprops)
   re = np.amax(absvec)
   print("The spectral radius is:")
   print(re,'\n')
   if re>1:
      print("The spectral raduis is greater than 1, the method cannot converge")
      sys.exit()
   
   print("The constant vector is:")
   c = np.matmul(dinv, b)
   print(c, '\n')

   def jacob (A,b, iter,xa=guess):
    err=tol+1
    space = ("      "*2)
    barr = ("|")
    print("0", barr, space, barr, guess)
    t.add_row([" | "+str(0)+" | ","     |     ", str(guess)+" | "])
    for i in range(1, iter):
        prev = np.matmul(T,xa)
        xi = np.add(prev,c)
        err = LA.norm(xi-xa)
        errf = "{:e}".format(err)
        print(i, barr, errf, barr, xi)
        t.add_row([" | "+str(i)+" | ",str(errf)+" | ", str(xi)+" | "])
        xa = xi
        if err<tol:
            break
    return xi
   sol = jacob(A,b,iter,xa=guess)
   print("El vector solución es:")
   print(sol)
   context = {'vec':t}
   return render(request, "JacobiR.html", context)

def SorRender(request):
   return render(request, "Sor.html", {})

def Sor(request):
   A = request.POST.get('matriz', None)
   b = request.POST.get('vectorB', None)
   guess = request.POST.get('vectorGuess', None)
   iter = request.POST.get('iteraciones', None)
   tol = request.POST.get('tolerancia', None)
   w = request.POST.get('relajacion', None)

   tol = float(tol)
   iter = int(iter)
   w = float(w)

   if w<=1:
    print('The ralaxing factor must be greater than 1\n')
    sys.exit()
   
   guess = np.matrix(guess)
   guess = np.array(guess)

   b = np.matrix(b)
   b = np.array(b)

   A = np.matrix(A)
   A = np.array(A)

   t = PrettyTable(['iterations','Error','xi'])
   t.border = False
   t.header = False
   norma = LA.norm(A, ord=2, axis=None, keepdims=False)
   inv = LA.linalg.inv(A)
   norminv = LA.norm(inv, ord=2, axis=None, keepdims=False)
   cond = norma*norminv

   deter=LA.linalg.det(A)
   if deter==0:
      print('\n'"The determinant of the matrix is 0, the problem hasn't solution")
      sys.exit()
   
   n = len(b)
   d = np.diag(np.diag(A))
   l = d - np.tril(A)
   u = d - np.triu(A)
   wl = w*l
   print("The transition matrix is:")
   dwl = np.subtract(d, wl)
   dinv = LA.linalg.inv(dwl)
   dprev = (1-w)
   dpos = dprev * d
   wu = w*u
   tprev = np.add(dpos, wu)
   T = np.matmul(dinv,tprev)
   print(T,'\n')

   vecprops = LA.eigvals(T)
   absvec = np.absolute(vecprops)
   re = np.amax(absvec)
   print("The spectral radius is:")
   print(re,'\n')
   if re>1:
      print("The spectral raduis is greater than 1, the method cannot converge")
      sys.exit()
   
   print("The constant vector is:")
   cprev = w * dinv
   c = np.matmul(cprev, b)
   print(c, '\n')

   def Sorr (A,b, iter,xa=guess):
    err=tol+1
    space = ("      "*2)
    barr = ("|")
    print("0", barr, space, barr, guess)
    t.add_row([" | "+str(0)+" | ","     |     ", str(guess)+" | "])
    for i in range(1, iter):
        prev = np.matmul(T,xa)
        xi = np.add(prev,c)
        err = LA.norm(xi-xa)
        errf = "{:e}".format(err)
        print(i, barr, errf, barr, xi)
        t.add_row([" | "+str(i)+" | ",str(errf)+" | ", str(xi)+" | "])
        xa = xi
        if err<tol:
            break
    return xi
   sor = Sorr(A,b,iter,xa=guess)
   print("El vector solución es:")
   print(sor)
   context = {'vec':t}
   return render(request, "SorR.html", context)

def GaussSeidelRender(request):
   return render(request, "GaussSeidel.html", {})

def GaussSeidel(request):
   A = request.POST.get('matriz', None)
   b = request.POST.get('vectorB', None)
   guess = request.POST.get('vectorGuess', None)
   iter = request.POST.get('iteraciones', None)
   tol = request.POST.get('tolerancia', None)

   tol = float(tol)
   iter = int(iter)
   
   guess = np.matrix(guess)
   guess = np.array(guess)

   b = np.matrix(b)
   b = np.array(b)

   A = np.matrix(A)
   A = np.array(A)

   t = PrettyTable(['iterations','Error','xi'])
   t.border = False
   t.header = False
   norma = LA.norm(A, ord=2, axis=None, keepdims=False)
   inv = LA.linalg.inv(A)
   norminv = LA.norm(inv, ord=2, axis=None, keepdims=False)
   cond = norma*norminv

   deter=LA.linalg.det(A)
   if deter==0:
      print('\n'"The determinant of the matrix is 0, the problem hasn't solution")
      sys.exit()
   
   n = len(b)
   d = np.diag(np.diag(A))
   l = d - np.tril(A)
   u = d - np.triu(A)
   print("The transition matrix is:")
   dinv = LA.linalg.inv(d-l)
   T = np.matmul(dinv,u)
   print(T,'\n')

   vecprops = LA.eigvals(T)
   absvec = np.absolute(vecprops)
   re = np.amax(absvec)
   print("The spectral radius is:")
   print(re,'\n')
   if re>1:
      print("The spectral raduis is greater than 1, the method cannot converge")
      sys.exit()
   
   print("The constant vector is:")
   c = np.matmul(dinv, b)
   print(c, '\n')

   def Gauss (A,b, iter,xa=guess):
    err=tol+1
    space = ("      "*2)
    barr = ("|")
    print("0", barr, space, barr, guess)
    t.add_row([" | "+str(0)+" | ","     |     ", str(guess)+" | "])
    for i in range(1, iter):
        prev = np.matmul(T,xa)
        xi = np.add(prev,c)
        err = LA.norm(xi-xa)
        errf = "{:e}".format(err)
        print(i, barr, errf, barr, xi)
        t.add_row([" | "+str(i)+" | ",str(errf)+" | ", str(xi)+" | "])
        xa = xi
        if err<tol:
            break
    return xi
   sol = Gauss(A,b,iter,xa=guess)
   print("El vector solución es:")
   print(sol)
   context = {'vec':t}
   return render(request, "GaussSeidelR.html", context)

#Splinesssssssssssssssssssssssssssssssssssssssssssss

def CubicSplineRender(request):
   return render(request, "CubicSpline.html", {})

def CubicSpline(request):
   x = request.POST.get('vectorX', None)
   f = request.POST.get('vectorY', None)
   

   x = np.array(x)
   f = np.array(f)

   x = x.tolist()
   f = f.tolist()

   xi = ast.literal_eval(x)
   fi = ast.literal_eval(f)
   resolution = 10



   def Cubic(xi, yi):
    n = len(xi)
    # h values
    h = np.zeros([n-1])
    for j in range(0, n-1, 1):
        h[j] = xi[j+1]-xi[j]
    # Equations system
    A = np.zeros([n-2, n-2])
    B = np.zeros([n-2])
    S = np.zeros([n])
    A[0, 0] = 2*(h[0]+h[1])
    A[0, 1] = h[1]
    B[0] = 6*((yi[2]-yi[1])/h[1] - (yi[1]-yi[0])/h[0])
    for i in range(1, n-3, 1):
        A[i, i-1] = h[i]
        A[i, i] = 2*(h[i]+h[i+1])
        A[i, i+1] = h[i+1]
        B[i] = 6*((yi[i+2]-yi[i+1])/h[i+1] - (yi[i+1]-yi[i])/h[i])
    A[n-3, n-4] = h[n-3]
    A[n-3, n-3] = 2*(h[n-3]+h[n-2])
    B[n-3] = 6*((yi[n-1]-yi[n-2])/h[n-2] - (yi[n-2]-yi[n-3])/h[n-3])
    # Solve equations system
    r = np.linalg.solve(A, B)
    print("r: ", r)
    # S
    for j in range(1, n-1, 1):
        S[j] = r[j-1]
    S[0] = 0
    S[n-1] = 0

    # coefficients
    a = np.zeros([n-1])
    b = np.zeros([n-1])
    c = np.zeros([n-1])
    d = np.zeros([n-1])
    for j in range(0, n-1, 1):
        a[j] = (S[j+1]-S[j])/(6*h[j])
        b[j] = S[j]/2
        c[j] = (yi[j+1]-yi[j])/h[j] - (2*h[j]*S[j]+h[j]*S[j+1])/6
        d[j] = yi[j]
    # Polynomial spline
    x = sp.Symbol('x')
    polynomial = []
    for j in range(0, n-1, 1):
        pSection = a[j]*(x-xi[j])**3 + b[j]*(x-xi[j])**2 + c[j]*(x-xi[j]) + d[j]
        pSection = pSection.expand()
        polynomial.append(pSection)
    return(polynomial)

   n = len(xi)
   polynomial = Cubic(xi, fi)
   xSpline = np.array([])
   ySpline = np.array([])
   section = 1
   while not(section >= n):
    pSection = polynomial[section-1]
    pxSection = sp.lambdify('x', pSection)

    a = xi[section-1]
    b = xi[section]
    xSection = np.linspace(a, b, resolution)
    ySection = pxSection(xSection)

    xSpline = np.concatenate((xSpline, xSection))
    ySpline = np.concatenate((ySpline, ySection))
    section = section + 1

   t = PrettyTable(['polym','interval'])
   t.header = False
   t.border = False
   # SALIDA
   print('Polinomyal sections: ')
   for section in range(1, n, 1):
    print(' x = ['+str(xi[section-1])
          + ','+str(xi[section])+']')
    print(str(polynomial[section-1]))
    t.add_row([str(polynomial[section-1])," |  x = ["+str(xi[section-1])+","+str(xi[section])+"]"])
   
   
   context = {'vec':t}
   print(t)


   return render(request, "CubicSplineR.html", context)


def QuadraticSplineRender(request):
   return render(request, "QuadraticSpline.html", {})

def QuadraticSpline(request):
   x = request.POST.get('vectorX', None)
   y = request.POST.get('vectorY', None)

   x = np.array(x)
   y = np.array(y)

   x = x.tolist()
   y = y.tolist()

   x1 = ast.literal_eval(x)
   y1 = ast.literal_eval(y)

   points = {"X":x1,"Y":y1}

   result = compute(points)
   print(result) 
   print(type(result))
   context = {'Resultado': result} 
   return render(request, "QuadraticSplineR.html", context)

def LinearSplineRender(request):
   return render(request, "LinearSpline.html", {})

def LinearSpline(request):
   x = request.POST.get('vectorX', None)
   y = request.POST.get('vectorY', None)

   x = np.array(x)
   y = np.array(y)

   x = x.tolist()
   y = y.tolist()

   vector1 = ast.literal_eval(x)
   vector2 = ast.literal_eval(y)
   n= len(vector1)
   resul = np.zeros(n-1)
   mult = np.zeros(n-1)

   for j in range(1,n):
      resul[j-1]=(vector2[j]-vector2[j-1])/(vector1[j]-vector1[j-1])
      mult[j-1]=(resul[j-1]*(-vector1[j]))+vector2[j]
   print(resul)
   print(mult)
   final = zip(resul,mult)

   context = {'final': final} 
   return render(request, "LinearSplineR.html", context)

#Sustituciones

def RegSubstitution(U, z):
    n, m = U.shape
    Uz = np.concatenate((U, z), axis=1)
    x = []
    for j in range(0, n-1): 
        x.append(1)
    x.append(Uz[n-1][n]/Uz[n-1][n-1])
    for i in range(n-1,-1,-1):
        acum=0
        for p in range(i+1, n):
            acum += Uz[i][p]*x[p]
        x[i] = ((Uz[i][n]-acum))/Uz[i][i]
    return x

def ProSubstitution(L, b):
    n, m = L.shape
    Lb = np.concatenate((L, b), axis=1)
    x = []
    x.append( Lb[0][n] / Lb[0][0])
    for j in range(0, n-1): 
        x.append(1)
    for i in range(1, n):
        acum=0
        for p in range(0, i):
            acum += Lb[i][p]*x[p]
        x[i] = ((Lb[i][n]-acum))/Lb[i][i]
    return x

#Metodos XYZ jeje XD

def mult_matrix(M, N):
   tuple_N = zip(*N)
   return [[sum(el_m * el_n for el_m, el_n in zip(row_m, col_n)) for col_n in tuple_N] for row_m in M]

def pivot_matrix(M):
   m = len(M)
   id_mat = [[float(i ==j) for i in range(m)] for j in range(m)]                                                                                                                                                                         
   for j in range(m):
      row = max(range(j, m), key=lambda i: abs(M[i][j]))
      if j != row:
         # Swap the rows                                                                                                                                                                                                                            
         id_mat[j], id_mat[row] = id_mat[row], id_mat[j]
   return id_mat

def expandParams(X, Y):
   n = len(X)
   new_X = np.zeros(n)
   new_Y = np.zeros(n)
   for i in range(n):
      new_X[i] = X[i]
      new_Y[i] = Y[i]

   points = np.array((new_X, new_Y)).T
   return points

def compute(params):
   X = params["X"]
   Y = params["Y"]
   points = expandParams(X, Y)
   n = len(points) - 1
   points = np.array(points)
   matrix = np.zeros((n*3, n*3))
   indVector = np.zeros(n*3)

   j = 0
   k = 0
   for i in range(0, n*2, 2):
      matrix[i, j+0] = points[k, 0] ** 2
      matrix[i, j+1] = points[k, 0]
      matrix[i, j+2] = 1

      matrix[i+1, j+0] = points[k+1, 0] ** 2
      matrix[i+1, j+1] = points[k+1, 0]
      matrix[i+1, j+2] = 1

      j += 3
      k += 1

   j = 1
   k = 0
   for i in range(n*2, n*3-1):
      matrix[i][k + 0] = 2 * points[j, 0]
      matrix[i][k + 1] = 1

      matrix[i][k + 2+1] = - 2 * points[j, 0]
      matrix[i][k + 3+1] = - 1
      j += 1
      k += 3

   matrix[n*3-1, 0] = 1

   indVector[0] = points[0, 1]
   j = 1
   for i in range(1, n):
      indVector[j] = points[i, 1]
      indVector[j+1] = points[i, 1]
      j += 2


   solution = np.linalg.solve(matrix, indVector)
   function =  generateEquation(solution, points)
   return function


def generateEquation(coefficients, points):
   segmentFunction = []
   coefficients = np.round(coefficients, 2)
   n = len(points) - 1

   for i in range(0, n*3, 3):
      function = "{a}x^2 + {b}x + {c}".format(
         a=coefficients[i],
         b=coefficients[i+1],
         c=coefficients[i+2],
      )
      segmentFunction.append([function, "{x0} <= x <= {x1}".format(
         x0=points[i//3, 0],
         x1=points[i//3+1, 0]
      )])


   return segmentFunction