import math
import numpy as np
import copy
from sympy import *
 
def isSimmetrial(coefficients, numberOfEquation ):
  result = True;
  
  for i  in range(0,numberOfEquation):
    for j in range(i+1, numberOfEquation-1):
      if (coefficients[i][j] != coefficients[j][i]):
        result = False;
        break    
    if (not result):
      break
      
  return result
  
def meanmanual(listt):
    lsum = 0
    lenoflist = len(listt)
    for i in listt:
        lsum += i
    mean = lsum / lenoflist
    return float(mean)


def spectr_rad(matrix):
    eps = 0.001
    y0 = y1 = [1, 1, 1, 1]
    while True:
        y0 = y1
        y1 = np.dot(matrix, y0)
        l = meanmanual(y1) / meanmanual(y0)
        y1 = y1 / max(abs(y1))
        if max(abs(np.dot(matrix, y1) - l * y1)) / max(abs(y1)) < eps:
            print("\nСпектральный радиус = ", round(l, 3))
            return

def wrachenie(coefficients, numberOfEquation, solution, precision):
  result = 1;
  maxI = None 
  maxJ = None
  max = None
  fi = None
  
  matricaPoworota = np.zeros((numberOfEquation, numberOfEquation))
  
  temp = np.zeros((numberOfEquation, numberOfEquation))
  
  fault = 0.0;
  for i in range(0, numberOfEquation):
    for j in range(i+1, numberOfEquation):
      fault = fault + coefficients[i][j]*coefficients[i][j]
  
  fault = math.sqrt(2*fault)
  
  while (fault > precision):
    max = 0.0;
    for i in range(0, numberOfEquation):
      for j in range(i+1, numberOfEquation):
        if (coefficients[i][j] > 0 and coefficients[i][j] > max ):
          max = coefficients[i][j]
          maxI = i
          maxJ = j
        elif ( coefficients[i][j] < 0 and - coefficients[i][j] > max ):
          max = - coefficients[i][j]
          maxI = i
          maxJ = j
          
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        matricaPoworota[i][j] = 0
      matricaPoworota[i][i] = 1
      
    if (coefficients[maxI][maxI] == coefficients[maxJ][maxJ]):
      matricaPoworota[maxI][maxI] = matricaPoworota[maxJ][maxJ] = matricaPoworota[maxJ][maxI] = math.sqrt(2.0) / 2.0
      matricaPoworota[maxI][maxJ] = - math.sqrt(2.0) / 2.0
    else:
      fi = 0.5 * math.atan( ( 2.0 * coefficients[maxI][maxJ] ) / ( coefficients[maxI][maxI] - coefficients[maxJ][maxJ] ) )
      matricaPoworota[maxI][maxI] = matricaPoworota[maxJ][maxJ] = math.cos(fi)
      matricaPoworota[maxI][maxJ] = - math.sin(fi)
      matricaPoworota[maxJ][maxI] = math.sin(fi)
      
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        temp[i][j] = 0.0
    
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        for k in range(0, numberOfEquation):
          temp[i][j] = temp[i][j] + matricaPoworota[k][i] * coefficients[k][j]
    
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        coefficients[i][j] = 0.0
        
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        for k in range(0, numberOfEquation):
          coefficients[i][j] = coefficients[i][j] + temp[i][k] * matricaPoworota[k][j]
    
    fault = 0.0
    for i in range(0, numberOfEquation):
      for j in range(i+1, numberOfEquation):
        fault = fault + coefficients[i][j]*coefficients[i][j]
    
    fault = math.sqrt(2*fault)
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        temp[i][j] = 0.0
    
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        for k in range(0, numberOfEquation):
          temp[i][j] = temp[i][j] + solution[i][k] * matricaPoworota[k][j]
    
    for i in range(0, numberOfEquation):
      for j in range(0, numberOfEquation):
        solution[i][j] = temp[i][j]
    
    result+=1
  return result, coefficients, solution
 
#чебышев
  
N = 14 #вариант

X = [0.3+0.1*N, 0.4+0.1*N, 0.5+0.1*N, 0.6+0.1*N, 0.7+0.1*N, 0.8+0.1*N, 0.9+0.1*N, 1.0+0.1*N, 1.1+0.1*N, 1.2+0.1*N]

Y = [0.5913, 0.63 + N / 17, 0.7162, 0.8731, 0.9574, 1.8 - math.cos(N / 11), 1.3561, 1.2738, 1.1 + N / 29, 1.1672]

m = 3
n = 9
#табличные значения
P =  np.array([[1, 9, 6, 42],
              [1, 7, 2, -14],
              [1, 5, -1, -35],
              [1, 3, -3, -31],
              [1, 1, -4, -12],
              [1, -1, -4, 12],
              [1, -3, -3, 34],
              [1, -5, -1, 35],
              [1, -7, 2, 14],
              [1, -9, 6, -42]]) 

k = [i for i in range(m+1)]
h = 0.1 #шаг
ck = np.dot(np.transpose(Y), P)
# print(ck)
sk = np.sum(P ** 2, axis = 0)
# print(sk)
bk = np.round(ck/sk,3)
# print(bk)
Qt = bk * P[0] #домножаем на наш множитель, чтобы вернуться к ортогональному полиному
# print(Qt)
x, t, N  = symbols('x t N')
N = 9
t = (x-X[0]) / h
p19 = 1 - 2/N*t
p29 = 1 - 6*t/9 + 6*t*(t-1)/(N*(N-1))
p39 = 1 - 12*t/N + 30*t*(t-1)/(N*(N-1)) - 20*t*(t-1)*(t-2)/(N*(N-1)*(N-2))
Q = Qt[0] + Qt[1]*p19 + Qt[2]*p29 + Qt[3]*p39

#округляем
Q_round = expand(Q) # расскрываем скобки
for a in preorder_traversal(expand(Q)):
    if isinstance(a, Float):
        Q_round = Q_round.subs(a, round(a, 3))
print("Многочлен Чебышева ", end = '')
print(Q_round)

print("значения функции в хi + h/2:")
for i in range(10):
    print(f"x{i} = {round(X[i]+h/2,4)} y = {round(Q_round.subs(x,X[i]+h/2),4)}")

#Вращение
size = 4 #размер матрицы

coefficients = np.array([[1.6, 1.6, 1.7, 1.8], [1.6, 2.6, 1.3, 1.3], [1.7, 1.3, 3.6, 1.4], [1.8, 1.3, 1.4, 4.6]])

solution = np.zeros((size, size))
for i in range(0, size):
  solution[i][i] = 1

precision = 0.0001
if (not isSimmetrial(coefficients, size)):
  print("Матрица не симметричная")
else:
  steps, coefficients2, solution2 = wrachenie(coefficients, size, solution, precision)

  # print(coefficients2)
  # print()
  # print(solution2)

  for i in range(0, size):
    print()
    print("Собственный вектор k=", i + 1)
    for j in range(0, size):
      print("%.4f     " % solution[j][i])

  print("\nСобственные значения:")
  for i in range(0, size):
    print("%.4f     " % coefficients[i][i])
  
  # print("\nОбщее число шагов: ", steps)
  
#спектральный радиус
spectr_rad(coefficients)
print()