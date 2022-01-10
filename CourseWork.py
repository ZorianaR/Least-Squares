import matplotlib.pyplot as plt
from numpy import *


def create_equations(X,Y,m):
    """Повертає матрицю системи рівнянь А та вектор вільних членів B,
       що відповідають поліному степеня m
    """
    n=len(X)
    assert n==len(Y)
    A=[[0 for i in range(m+1)]for j in range(m+1)]
    B=[0 for i in range(m+1)]
    sums_of_X=[0 for i in range(2*m+1)]
    sums_of_X[0]=n

    for i in range(1,2*m+1):
        for j in range(n):
            sums_of_X[i]+=X[j]**i
   
    for i in range(m+1):
        for j in range(m+1):
            A[i][j]=sums_of_X[i+j]
        
    for i in range(m+1):
        for j in range(n):
            B[i]+=Y[j]*(X[j]**i)
    return A, B

def max_in_column(A,i,j):
    """Повертає номер рядка з найбільшим по модулю елементом матриці А
       у j-му стовпці. Пошук починається з і-того рядка
    """
    m=len(A)
    b=i
    maxj=abs(A[i][j])
    for k in range(i,m):
        if abs(A[k][j])>maxj:
            maxj=abs(A[k][j])
            b=k
    return b
   

def Gauss(A,B):
    """Повертає розв'язок системи лінійних рівнянь з матрицею коефіцієнтів А
       та стовпцем вільних членів В, знайдений методом Гаусса з вибором головного 
       елемента у стовпці
    """
    m=len(A)
    assert m==len(B)
    # Прямий хід 
    for j in range(m):
        b=max_in_column(A,j,j)
        A[b],A[j]=A[j],A[b]
        B[b],B[j]=B[j],B[b]
        for i in range(1+j,m):
            k=-A[i][j]/A[j][j]
            for l in range(j,m):
                A[i][l]=k*A[j][l]+A[i][l]
            B[i]=k*B[j]+B[i]
    # Зворотній хід
    for j in range(m-1,0,-1):
        for i in range(0,j):
            k=-A[i][j]/A[j][j]
            A[i][j]=k*A[j][j]+A[i][j]
            B[i]=k*B[j]+B[i]
    # Знаходження невідомих
    koef=[]
    for j in range(m):
        koef.append(B[j]/A[j][j])
    return koef
        
def polinom(x,koef):
    """Повертая значення полінома з коефіцієнтами koef у точці х"""
    s=koef[0]
    for i in range(1,len(koef)):
        s+=koef[i]*x**i
    return s

def polinom_to_string(koef):
    """Повертає стрічкове представлення полінома"""
    s=f'{koef[0]: ^3.3f}'
    for i in range(1,len(koef)):
        s+=f'+{koef[i]: ^3.3f}*x**{i}'
    return s
    
def square_error(Y, approximation):
    """Повертає середньоквадратичну похибку  апроксимації"""
    s=0
    n=len(Y)
    for i in range(n):
        s+=(Y[i]-approximation[i])**2
    return s/n

def max_abs_error(Y, approximation):
    """Повертає найбільшу абсолютну похибку, та номер точки у якій вона досягяється"""
    x=0
    max_y=abs(Y[0]-approximation[0])
    for i in range(1,len(Y)):
        if abs(Y[i]-approximation[i])>max_y:
            max_y=abs(Y[i]-approximation[i])
            x=i
    return max_y, x

m=3#степінь полінома
choice=input("Наближуємо \n a)таблично задану функцію \t b)аналітично задану функцію \n")
if choice=='a':
    X=list(map(float,input("Введіть масив X через пробіли ").split()))
    Y=list(map(float,input("Введіть масив Y через пробіли ").split()))
elif choice=='b':
    f=input ('Введіть функцію f(x)=')
    code="""
def f(x):
    return %s
""" % f
    exec(code)
    a=float(input("Введіть початок інтервалу "))
    b=float(input("Введіть кінець інтервалу "))
    n=float(input("Введіть кількість точок розбиття "))
    dx=(b-a)/n
    X=list(arange(a,b+dx,dx))
    Y=[f(x) for x in X]


choice=input("Зашумлюємо дані?(y/n) ")
if choice=='y':
    Y1=Y+random.normal(0.0, abs(max(Y)-min(Y))*0.005, len(Y))
else:
    Y1=Y

A, B=create_equations(X,Y1,m)
a_i=Gauss(A, B)
Y_approx=[polinom(x,a_i) for x in X]
serror=square_error(Y1, Y_approx)
abs_error, i=max_abs_error(Y1, Y_approx)


plt.scatter(X,Y1,marker='o',label='Точки по яких будується наближення',color='#d83d3d')
plt.plot(X,Y, label='Введена функція')
plt.plot(X, Y_approx, label=polinom_to_string(a_i))
plt.title(f'Середньоквадратична похибка={serror: ^.5f}, найбільша абсолютна похибка={abs_error: ^.5f} у точці x={X[i]: ^.4f}')
plt.legend()
plt.grid(True)
plt.show()




