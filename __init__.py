import re
from sympy import Rational
from tkinter import messagebox, ttk, font
from tkinter import *
from FO import FO
from Restriccion import Restriccion
from Tabla import Matriz

class MethodSimplex():
    def __init__(self):
        self.restricciones=[]
        self.__raiz = Tk()
        self.__raiz.title("Method Simplex")
        self.__frame = Frame(self.__raiz,width="1000", height="600", pady=20, padx=100)
        self.vista_1()
        self.__raiz.mainloop()

    def vista_1(self):
        self.clear()
        self.__tabla = None
        self.fo = None
        self.res = None
        self.igu = None
        self.__frame.pack()
        self.numVar = StringVar()
        self.numRes = StringVar()

        Label(self.__frame, text="METHOD SIMPLEX", font=font.Font(family="Helvetica", size=30, weight="bold")).grid(pady=15, columnspan=2)

        Label(self.__frame, text="Ingrese el número de variables:").grid(row=1, column=0)
        Entry(self.__frame, textvariable=self.numVar).grid(row=1, column=1)
        
        Label(self.__frame, text="Ingrese el número de restricciones:  ").grid(row=2, column=0)
        Entry(self.__frame, textvariable=self.numRes).grid(row=2, column=1)

        _btnOk=Button(self.__frame, text="OK", command=self.btnOk)
        _btnOk.grid(pady=15, columnspan=2)
        _btnOk.config(width=15)
        
        return self.numVar, self.numRes

    def vista_2(self):
        self.clear()

        self.fo = [StringVar() for i in range(self.numVar)]
        self.res = [[StringVar() for i in range(self.numVar)] for i in range(self.numRes)]
        self.igu = [[StringVar() for i in range(2)] for i in range(self.numRes)]

        Label(self.__frame, text="METHOD SIMPLEX").grid(pady=15, columnspan=self.numVar*2+2)
        Label(self.__frame, text="Función objetivo:").grid(pady=15, row = 1,columnspan=self.numVar*2+2)

        self.max_min = ttk.Combobox(self.__frame, width=6)
        self.max_min['values'] = ('Máx Z', 'Min Z')
        self.max_min.grid(row=2, column=0)
        self.max_min.current(0)
        Label(self.__frame, text="=").grid(row=2, column=1)

        colum=2
        for i in range(self.numVar):
            aux = Entry(self.__frame, textvariable=self.fo[i])
            aux.grid(row=2, column=colum)
            aux.config(width=10)
            colum+=1
            Label(self.__frame, text="X_"+str(i+1)).grid(row=2, column=colum)
            colum+=1

        Label(self.__frame, text="Sujeto a:").grid(pady=15, row = 3,columnspan=self.numVar*2+2)

        for i in range(self.numRes):
            colum=0
            for j in range(self.numVar):
                aux = Entry(self.__frame, textvariable=self.res[i][j])
                aux.grid(row=i+4, column=colum)
                aux.config(width=10)
                colum+=1
                Label(self.__frame, text="X_"+str(j+1)).grid(row=i+4, column=colum)
                colum+=1
            self.igu[i][0] = ttk.Combobox(self.__frame, width=3)
            self.igu[i][0]['values'] = ('<=', '=', '>=')
            self.igu[i][0].grid(row=i+4, column=colum)
            self.igu[i][0].current(0)
            aux = Entry(self.__frame, textvariable=self.igu[i][1])
            aux.grid(row=i+4, column=colum+1)
            aux.config(width=10)
        
        aux = ""
        for i in range(self.numVar):
            aux += "X_" + str(i+1)
            if i +1 == self.numVar:
                aux += " >= 0"
            else:
                aux += ", "

        Label(self.__frame, text=aux).grid(row=self.numRes+4, columnspan=self.numVar*2+2)

        colum = int((self.numVar*2+2)/2)

        _btnOk=Button(self.__frame, text="Atras", command=self.vista_1)
        _btnOk.grid(pady=15, row=self.numRes+5, columnspan=colum, sticky=E)
        _btnOk.config(width=8)

        _btnOk=Button(self.__frame, text="Siguiente", command=self.btnV_2)
        _btnOk.grid(pady=15, row=self.numRes+5, columnspan=colum, column=colum,sticky=W)
        _btnOk.config(width=8)

        return self.fo, self.res, self.igu

    def vista_3(self):
        self.clear()

        fon = font.Font(family="Helvetica", size=10, weight="bold")

        pivote = [self.__tabla.columna_pivote(), self.__tabla.fila_pivote()]
        print(pivote)

        Label(self.__frame, font=fon, text="METHOD SIMPLEX").grid(pady=15, columnspan=len(self.__tabla.fo.coef)+2, sticky=NSEW)
        Label(self.__frame, bg="black").grid(row=2, column=0, sticky=NSEW)
        Label(self.__frame, font=fon, relief=RIDGE, text="Cj", bg="red").grid(row=2, column=1, sticky=NSEW)

        for i in range(self.__tabla.num_hol + 1):
            Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.fo.coef[i]), bg="red").grid(row=2, column=i+2, sticky=NSEW)
            self.__raiz.update()
        
        for i in range(self.__tabla.num_hol + 3):
            Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.var[i]), bg="grey").grid(row=3, column=i, sticky=NSEW)
            self.__raiz.update()

        for i in range(self.__tabla.num_res):
            Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.c_v_var[0][i])).grid(row=i+4, column=0, sticky=NSEW)
            Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.c_v_var[1][i]), bg="yellow").grid(row=i+4, column=1, sticky=NSEW)
            for j in range(len(self.__tabla.tabla[i])):
                if j == pivote[0] or i == pivote[1] and pivote[0] != 0:
                    Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.tabla[i][j]), bg="blue").grid(row=i+4, column=j+2, sticky=NSEW)
                else:
                    Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.tabla[i][j])).grid(row=i+4, column=j+2, sticky=NSEW)
                self.__raiz.update()
        
        Label(self.__frame, bg="black").grid(row=i+5, column=0, sticky=NSEW)
        Label(self.__frame, font=fon, relief=RIDGE, text="Zj", bg="green").grid(row=i+5, column=1, sticky=NSEW)
        for j in range(len(self.__tabla.zj)):
            if pivote[0] == 0 and j == 0:
                Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.zj[j]), bg="blue").grid(row=i+5, column=j+2, sticky=NSEW)
            elif j == pivote[0]:
                Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.zj[j]), bg="blue").grid(row=i+5, column=j+2, sticky=NSEW)
            else:
                Label(self.__frame, font=fon, relief=RIDGE, text=str(self.__tabla.zj[j])).grid(row=i+5, column=j+2, sticky=NSEW)
            self.__raiz.update()

        colum = int((len(self.__tabla.fo.coef)+3)/2)

        _btnOk=Button(self.__frame, text="Atras", command=self.vista_1)
        _btnOk.grid(pady=15, row=self.numRes+5, columnspan=colum, sticky=E)
        _btnOk.config(width=8)

        _btnOk=Button(self.__frame, text="Siguiente", command=self.operar)
        _btnOk.grid(pady=15, row=self.numRes+5, columnspan=colum, column=colum,sticky=W)
        _btnOk.config(width=8)

    def btnV_2(self):
        sw = True
        expresion = r'^(([0-9]+)|([0-9]+[\.][0-9]+)|([0-9]+[\/][1-9][0-9]*))$'
        for aux in self.fo:
            if False == re.search(expresion, aux.get()):
               sw = False 
               break
        for i in self.res:
            for j in i:
                if not re.search(expresion, j.get()):
                    sw = False 
                    break
        for aux in range(self.numRes):
            if not re.search(expresion, self.igu[aux][1].get()):
               sw = False 
               break

        self.restricciones = []
        if sw:
            fo = []

            for aux in self.fo:
                fo.append(Rational(aux.get()))

            for i in range(self.numRes):
                res = []
                for aux in self.res[i]:
                    res.append(Rational(aux.get()))
                self.restricciones.append(Restriccion(res, self.igu[i][0].get(), Rational(self.igu[i][1].get())))

            if self.max_min.current() == 0:
                sw=True
            else:
                sw=False

            self.fo = None
            self.res = None
            self.igu = None
            self.__FO = FO(fo, sw)
            self.__tabla = Matriz(self.__FO, self.restricciones)
            self.vista_3()
        else:
            messagebox.showinfo("ERROR", "Hay algún número no valido.")

    def operar(self):
        sw = self.__tabla.operar()
        if not sw:
            self.vista_3()
        elif sw:
            aux = "\nZ = "+str(self.__tabla.z)
            for i in range(self.__tabla.num_res):
                if 'x_' in self.__tabla.c_v_var[1][i]:
                    aux += "\n" + self.__tabla.c_v_var[1][i] + " = " + str(self.__tabla.tabla[i][0])
            messagebox.showinfo("Solución", "La solución más óptima es:" + aux)
        else:
            messagebox.showinfo("Solución", "Este problema no tiene solución.")

    def btnOk(self):
        if re.search(r'^[0-9]+$', self.numVar.get()) and re.search(r'^[0-9]+$', self.numRes.get()):
            self.numVar, self.numRes = int(self.numVar.get()), int(self.numRes.get())
            self.vista_2()
        else:
            messagebox.showinfo("ERROR", "La expresión debe ser solo números.")

    def clear(self):
        for a in self.__frame.grid_slaves():
            a.destroy()
        for a in self.__raiz.grid_slaves():
            a.destroy()

a = MethodSimplex()