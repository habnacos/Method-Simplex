#######################################################
# 
# Restriccion.py
# Python implementation of the Class Restriccion
# Generated by Enterprise Architect
# Created on:      10-nov.-2018 11:31:47 p. m.
# Original author: habna
# 
#######################################################
from sympy import Rational

class Restriccion():
    def __init__(self, coef, igualdad, valor):
        self.__valor = valor
        self.igualdad = igualdad
        self.coef = [self.__valor] + coef                       # Coeficientes de la restriccion.
        self._var_H = self.__holgura()                          # Variables de holgura.

    def __holgura(self):
        """ Método privado qu retorna la variable de holgura.
        
        A partir de la condicion:
                '<='        Variable 1S_i
                '>='        Variable -1S_i + 1R_i
                '='         Variable 1R_i
        """
        if '<=' in self.igualdad:
            return [1, None]                                      # Retorna S.
        elif '>=' in self.igualdad:
            return [-1, 1]                                        # Retorna -S y R.
        elif '=' in self.igualdad:
            return [None, 1]                                      # Retorna R.
        return [None, None]

    def holgura(self, pos, cantRes):
        """Obtiene los coeficientes de las variables de olgura.

        Recibe como parámetro la cantidad de restricciones y la posición
        donde se hubica esta. Método público.
        """
        coeficientes = self.coef                                              # Coeficientes de las variables.
        for i in (range(cantRes)):
            if i is pos:
                if self._var_H is ([1, None] or [None, 1]):
                    coeficientes.append(1)                                      # Coeficiente S o R.
                elif self._var_H is [-1, 1]:
                    coeficientes.append(-1)                                     # Coeficiente -S.
                    coeficientes.append(1)                                      # Coeficiente R.
            else:
                coeficientes.append(0)                                          # Coeficiente 0.
            if len(coeficientes) == cantRes:
                break
        return coeficientes, self._var_H                                       # Coeficientes de holgura 

    def operDiv(self, divisor):
        """Divide toda la función para obtener 1 en Aij.
        """
        for i in range(len(self.coef)):
            self.coef[i] /= Rational(divisor)
        return  self.coef

    def oper(self, columna, fila):
        """Realiza las operaciones para obtener 0 en Aij.
        """
        mult = self.coef[columna]
        for i in range(len(self.coef)):
            self.coef[i] -= mult*fila.coef[i]
        return self.coef

    def teta(self, pos):
        """Obtiene el valor de teta θ

        Divide el valor de la restricción por el valor del coeficiente Aij.
        """
        if self.coef[pos] > 0:
            return Rational(self.coef[0] / self.coef[pos])
        else:
            return 0

    def getHolgura(self):
        return self._var_H