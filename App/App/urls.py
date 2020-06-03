"""App URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from App import views
#Listo
urlpatterns = [
    path('', views.index, name='index'),

    path('IncrementalSearch/', views.IncrementalSearch, name='IncrementalSearch'),
    path('IncrementalSearchRes/', views.IncrementalSearchRes, name='IncrementalSearchRes'),

    path('biseccion/', views.Biseccion, name='Biseccion'),
    path('biseccionres/', views.BiseccionRes, name='BiseccionRes'),

    path('fixedpoint/', views.fixedpoint, name='FixedPoint'),
    path('fixedpointres/', views.fixedpointres, name='fixedpointres'),
    
    path('falseposition/', views.falseposition, name='FalsePosition'),
    path('falsepositionres/', views.falsepositionres, name='falsepositionres'),

    path('newton/', views.newton, name='newton'),
    path('newtonres/', views.newtonres, name='newtonres'),

    path('secant/', views.secant, name='secant'),
    path('secantres/', views.secantres, name='secantres'),

    path('multipleroots/', views.multipleroots, name='MultipleRoots'),
    path('multiplerootsres/', views.multiplerootsres, name='multiplerootsres'),

    path('CroutRender/', views.CroutRender, name='CroutRender'),
    path('Crout/', views.Crout, name='Crout'),

    path('CholeskyRender/', views.CholeskyRender , name='CholeskyRender'),
    path('Cholesky/', views.Cholesky , name='Cholesky'),

    path('DoolittleRender/', views.DoolittleRender , name='DoolittleRender'),
    path('Doolittle/', views.Doolittle, name='Doolittle'),

    path('LUParcialRender/', views.LUParcialRender , name='LUParcialRender'),
    path('LUParcial/', views.LUParcial, name='LUParcial'),

    path('PartialPivotingRender/', views.PartialPivotingRender , name='PartialPivotingRender'),
    path('PartialPivoting/', views.PartialPivoting, name='PartialPivoting'),

    path('JacobiRender/', views.JacobiRender , name='JacobiRender'),
    path('Jacobi/', views.Jacobi, name='Jacobi'),

    path('SorRender/', views.SorRender , name='SorRender'),
    path('Sor/', views.Sor, name='Sor'),

    path('GaussSeidelRender/', views.GaussSeidelRender , name='GaussSeidelRender'),
    path('GaussSeidel/', views.GaussSeidel, name='GaussSeidel'),

    path('QuadraticSplineRender/', views.QuadraticSplineRender , name='QuadraticSplineRender'),
    path('QuadraticSpline/', views.QuadraticSpline, name='QuadraticSpline'),

    path('CubicSplineRender/', views.CubicSplineRender , name='CubicSplineRender'),
    path('CubicSpline/', views.CubicSpline, name='CubicSpline'),

    path('LinearSplineRender/', views.LinearSplineRender , name='LinearSplineRender'),
    path('LinearSpline/', views.LinearSpline, name='LinearSpline'),

    path('newtoninterp/', views.newtoninterp, name='newtoninterp'),
    path('newtoninterpres/', views.newtoninterpres, name='newtoninterpres'),

    path('vandermonde/', views.vandermonde, name='vandermonde'),
    path('vandermonderes/', views.vandermonderes, name='vandermonderes'),

    path('lagrange/', views.lagrange, name='lagrange'),
    path('lagrangeres/', views.lagrangeres, name='lagrangeres'),

    path('admin/', admin.site.urls),
]
