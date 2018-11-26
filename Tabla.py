from sympy import Symbol, Rational, Poly

class Matriz():
    def __init__(self, fo, restricciones):
        self.pivote = [0, 0]
        self.fo = fo                                    # Función objetivo.
        self.restricciones = restricciones              # Restricciones.
        self.num_res = len(restricciones)               # Número de restricciones.
        self.num_var = len(fo.coef)                     # Número de variables.
        self.num_hol = self.num_holgura()               # Número de variables de holgura.
        self.is_M, self.zj, self.tabla = self.matriz()  # Matriz de coeficientes.
        self.c_v_var = self.c_v_variables()             # Datos de C y V.
        self.imprimir()
        
    def num_holgura(self):
        """Obtiene el número de variables total."""
        cant = self.num_var + self.num_res              # Para cada restricción tiene una variable de holgura.
        for res in self.restricciones: 
            if res._var_H == [-1, 1]:                   # Método de la gran M.
                cant += 1
        return cant                                     # Retorna el número de variables de holgura.

    def matriz(self):
        """Obtiene una matriz de coeficientes."""
        self.var = ['C', 'V', 'B'] + ['x_'+str(i+1) for i in range(self.num_var)]
        self.identidad = []
        s, r = 1, 1
        is_M = False                                                            # Si se va a desarrollar con el método de la gran M.
        matriz = []                                                             # Matriz.                       
        zj = [0] + self.fo.coef                                                 # Zj.
        self.fo.coef = [0] + self.fo.coef
        count = 0                                                               # Método de la gran M.
        for i in (range(self.num_res)):                                         # Coeficientes de las restricciones.
            aux, hol = self.restricciones[i].holgura(i+count, self.num_hol + 1) # Coeficientes de las variables de Holgura.
            if hol[1] == 1:  
                aux[self.num_var+i+count+1] = 1
                self.var.append('R_'+str(r))
                self.identidad.append('R_'+str(r))
                r+=1
                if not self.fo.max:
                    zj.append(Symbol('M') * -1)                                  # Método de la gran M maximizar.
                    self.fo.coef.append(Symbol('M') * -1)
                else:
                    zj.append(Symbol('M'))                                       # Método de la gran M minimizar.
                    self.fo.coef.append(Symbol('M'))
                if hol[0] == -1:
                    count +=1
                    zj.append(0)  
                    self.fo.coef.append(0)
                    self.var.append('S_'+str(s))
                    s+=1
                    aux[self.num_var+i+count+1] = -1                             # Matriz identidad.
                is_M = True
            else:
                zj.append(0)                                                     # Zj Coeficiente variable de holgura.
                self.fo.coef.append(0)
                self.var.append('S_'+str(s))
                self.identidad.append('S_'+str(s))
                s+=1
                aux[self.num_var+i+count+1] = 1
            matriz.append(aux)                                                   # Matriz identidad.
            self.restricciones[i].coef = aux                                     # Coeficientes más holgura.
        if is_M:
            aux = self.Z()
            for i in range(len(aux)):                             # Coeficientes M despejada.                                              # Igualar los coeficientes a la función objetivo.
                zj[i] += aux[i]
                self.fo.coef[i] = 1 * zj[i]
                print(zj[i], aux[i])
        else:
            for i in range(len(zj)):
                zj[i] *= -1
        return is_M, zj, matriz                                                  # Retorna la matriz de coeficientes.

    def c_v_variables(self, columna=0, fila=0):
        '''Obtiene los valores actualizados de C y V.'''
        _c = []                                                             # Coeficiente.
        _v = []                                                             # Variable.
        if columna == 0:
            for i in self.identidad:
                _v.append(i)                                                # Variables de Holgura.
                _c.append(0)                                                # Coeficientes en 0.
        else:
            _c = self.c_v_var[0]                                            # Coeficientes.
            _v = self.c_v_var[1]                                            # Variables de Holgura.
            _c[fila] = self.fo.coef[columna]                                    # Cambio de coeficiente.
            _v[fila] = self.var[columna + 2]                                       # Cambio de variable.
        return [_c, _v]                                                     # C | V.

    def Cj_Zj(self):
        '''Obtiene los valores de Cj - Zj.'''
        cant = self.num_hol +1                                              # Total de variables.
        zj = [0 for i in range(cant)]                                       # Zj.
        cj = [0 for i in range(cant)]                                       # Cj.
        for i in range(cant):
            for j in range(self.num_res):
                res = self.restricciones[j]
                zj[i] += res.coef[i] * self.c_v_var[0][j]               # Zj (sumatoria de Aj * Cj).
            if self.fo.max:
                cj[i] = zj[i] - self.fo.coef[i]                                  # Cj - Zj.
            else:
                cj[i] = self.fo.coef[i] - zj[i]                                # Cj - Zj.
        self.zj = cj
        return self.zj
    
    def imprimir(self):
        print("#######################################################")
        print("#") 
        print('#', ['Cj'], self.fo.coef)
        for i in range(self.num_res):
            print('#', self.c_v_var[0][i], self.c_v_var[1][i], self.tabla[i])
        print('#', ['Zj'], self.zj)
        print("#") 
        print("#######################################################")

    def operar(self):
        if self.pivote[1] != 0:
            res = self.restricciones[self.pivote[0]]
            if res.coef[self.pivote[1]] > 0:
                res.operDiv(res.coef[self.pivote[1]])
            else:
                return 'error'
            for i in range (self.num_res):
                if i != self.pivote[0]:
                    self.restricciones[i].oper(self.pivote[1], res)
                self.tabla[i] = self.restricciones[i].coef
            self.c_v_variables(columna=self.pivote[1], fila=self.pivote[0])
            if self.Cj_Zj() == True:
                print('No tiene solución.')
            self.pivote[1] = self.columna_pivote()
            self.imprimir()
            return False
        else:
            print("#######################################################")
            print("#") 
            print("# La solución más óptima es:")
            self.z = self.zj[0]
            print("# Z = " + str(self.zj[0]))
            for i in range(self.num_res):
                if self.c_v_var[0][i] != 0:
                    print("# " + self.c_v_var[1][i] + " = " + str(self.tabla[i][0]))
            print("#") 
            print("#######################################################")
            return True

    def columna_pivote(self):
        '''Obtien la posición del valor a Máximizar o Minimizar.'''
        aux = 0
        pos = 0
        if self.is_M:
            for i in range(1, len(self.zj)):
                try:
                    if aux > self.zj[i].coeff(Symbol('M')) and self.zj[i].coeff(Symbol('M')) < 0 and self.fo.max:
                        pos = i                                         # Obtener posición.
                        aux = self.zj[i].coeff(Symbol('M'))             # Método de la Gran M.
                    elif aux < self.zj[i].coeff(Symbol('M')) and self.zj[i].coeff(Symbol('M')) > 0 and not self.fo.max:                        
                        pos = i                                         # Obtener posición.
                        aux = self.zj[i].coeff(Symbol('M'))             # Método de la Gran M.
                except:
                    pass
            if pos == 0:
                for i in range(1, len(self.zj)):
                    try:
                        if aux > self.zj[i] and self.zj[i] < 0 and self.fo.max:
                            pos = i                                         # Obtener posición.
                            aux = self.zj[i]             # Método de la Gran M.
                        elif aux < self.zj[i] and self.zj[i] > 0 and not self.fo.max:                        
                            pos = i                                         # Obtener posición.
                            aux = self.zj[i]             # Método de la Gran M.
                    except:
                        pass
        else:
            for i in range(1, len(self.zj)):                                # Recorre Zj - Cj.
                if aux > self.zj[i] and self.zj[i] < 0 and self.fo.max:
                    pos = i                                         # Obtener posición.
                    aux = self.zj[i]
                elif aux < self.zj[i] and self.zj[i] > 0 and not self.fo.max:
                    pos = i                                         # Obtener posición.
                    aux = self.zj[i]
        self.pivote[1] = pos
        return pos

    def fila_pivote(self):
        '''Obtien la posición del valor a Máximizar o Minimizar.'''
        pos = 0
        for i in range(self.num_res):
            aux = self.restricciones[i].teta(self.pivote[1])
            if aux != 0:
                pos = i
                break
        for i in range(self.num_res):
            teta = self.restricciones[i].teta(self.pivote[1])
            if aux > teta and teta > 0:
                aux = teta
                pos = i
        self.pivote[0] = pos
        return pos

    def Z(self):
        """Obtiene los coeficientes de cada variable con el método de la gran M."""
        cj = [0 for i in range(self.num_hol+1)]                       # Arreglo de coefientes.
        for res in self.restricciones:
            if res.getHolgura()[1] == 1:                            # Método de la gran M.
                for i in range(self.num_hol+1):
                    print(res.coef[i])
                    if self.fo.max:
                        cj[i] += res.coef[i] * Symbol('M') * -1
                    else:
                        cj[i] += res.coef[i] * Symbol('M')
        return cj                                                   # Retorna los coeficientes despejados.