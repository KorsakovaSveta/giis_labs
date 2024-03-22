import tkinter as tk
import numpy as np
import math, cmath
import numpy as np
class GraphicEditor:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.window = tk.Tk()

        self.canvas = tk.Canvas(self.window, width=self.width, height=self.height)
        
        self.mainmenu = tk.Menu(self.window)
        self.window.config(menu=self.mainmenu)
        self.line_menu = tk.Menu(self.mainmenu, tearoff=0)
        self.line_menu.add_command(label="ЦДА", command=self.activate_canvas_cda)
        self.line_menu.add_command(label="Брезенхем", command=self.activate_canvas_brezenhem)
        self.line_menu.add_command(label="Ву", command=self.activate_canvas_wu)
        self.mainmenu.add_cascade(label="Отрезки", menu=self.line_menu)
        self.mainmenu.add_cascade(label="Окружность", command=self.activate_canvas_circle)
        self.mainmenu.add_cascade(label="Эллипс", command=self.activate_canvas_ellipse)
        self.mainmenu.add_cascade(label="Гипербола", command=self.activate_canvas_giperbola)
        self.mainmenu.add_cascade(label="Парабола", command=self.activate_canvas_parabola)
        self.mainmenu.add_cascade(label="Кривая Безье", command=self.activate_canvas_curve_Bezie)
        self.mainmenu.add_cascade(label="Кривая Эрмита", command=self.activate_canvas_curve_Hermit)
        self.mainmenu.add_cascade(label="В-сплайн", command=self.activate_canvas_B_splain)
        
        #self.cda = tk.Button(self.window, text="Отрезки", command=self.activate_canvas)
        self.mainmenu.add_cascade(label="Отладка", command=self.debug_mode_toggle)
        self.mainmenu.add_cascade(label="Остановить отладку", command=self.delete_grid)
        #self.canvas.create_line(0, 2, self.window.winfo_screenwidth(), 2, fill="gray")
        #self.cda.pack()   
        #self.canvas.bind("<Button-1>", self.on_mouse_click)     
  
        
        self.canvas.pack()
       
        self.points = []
        self.derivatives = []
        self.debug_mode = False
        self.grid_size = 10


    def activate_canvas_cda(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_cda)

    def activate_canvas_brezenhem(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_brezenhem)
    
    def activate_canvas_wu(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_wu)

    def activate_canvas_circle(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_circle)

    def activate_canvas_ellipse(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_ellipse)
    
    def activate_canvas_giperbola(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_giperbola)

    def activate_canvas_parabola(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_parabola)

    def activate_canvas_curve_Bezie(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_curve_Bezie)
        
    def activate_canvas_curve_Hermit(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_curve_Hermit)

    def activate_canvas_B_splain(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_B_splain)

    def debug_mode_toggle(self):
        self.canvas.delete("all")
        self.debug_mode = not self.debug_mode
        self.draw_grid()
       
       

    def on_mouse_click_cda(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:           
            self.draw_line_cda()           
            self.points = []
    
    def on_mouse_click_brezenhem(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:           
            self.draw_line_brezenhem()           
            self.points = []
    
    def on_mouse_click_wu(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:           
            self.draw_line_wu()           
            self.points = []

    def on_mouse_click_circle(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:           
            self.draw_circle()           
            self.points = []

    def on_mouse_click_ellipse(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:           
            self.draw_ellipsis()           
            self.points = []

    def on_mouse_click_giperbola(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:           
            self.draw_giperbola()           
            self.points = []
    
    def on_mouse_click_parabola(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:           
            self.draw_parabola()           
            self.points = []
    
    def on_mouse_click_curve_Bezie(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 4:           
            self.draw_curve_Bezie()           
            self.points = []

    def on_mouse_click_curve_Hermit(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) ==2:   
            #self.calculate_derivatives()        
            self.draw_curve_Hermit()           
            self.points = []

    def on_mouse_click_B_splain(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) ==4:   
            #self.calculate_derivatives()        
            self.draw_B_splain()           
            self.points = []
    
            
    def draw_line_cda(self):
        #self.canvas.create_line(0, 2, self.window.winfo_screenwidth(), 2, fill="gray")
        #self.canvas.bind("<Button-1>", self.on_mouse_click)
        
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]

        length = max(abs(x2-x1), abs(y2-y1))
        
        dx = (x2 - x1)/length if length!= 0 else 0
        dy = (y2 - y1)/length if length!= 0 else 0
        sign_dx = np.sign(dx)
        sign_dy = np.sign(dy)

        x=x1 + 0.5*sign_dx
        y=y1 +0.5*sign_dy
        prev_x_grid=prev_y_grid=0
        for i in range(int(length+1)):
            if self.debug_mode:
                
                 # Округляем значения x и y до ближайшего кратного grid_size
                x_grid = x // self.grid_size
                y_grid = y // self.grid_size
                
                # Отображаем текущие значения x и y на сетке
                if x_grid!= prev_x_grid or y_grid!= prev_y_grid:
                    self.canvas.create_rectangle(x_grid*self.grid_size, y_grid*self.grid_size, (x_grid+1)*self.grid_size, (y_grid+1)*self.grid_size, fill='lightgray')
                    self.window.update()
                prev_x_grid = x_grid
                prev_y_grid = y_grid
                self.canvas.create_rectangle(x, y, x, y)
                self.window.update()
                #self.window.after(50)
            else:
                self.canvas.create_rectangle(x, y, x, y)
            x += dx
            y += dy
    
    


    def draw_line_brezenhem(self):
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]
        dx = abs(x2-x1)
        dy = abs(y2-y1)
        e = 2*dy - dx
        x=x1
        y=y1
        prev_x_grid=prev_y_grid=0
        sx = -1 if x1 > x2 else 1
        sy = -1 if y1 > y2 else 1
        while x!= x2 or y!= y2:
            if self.debug_mode:
                
                # Округляем значения x и y до ближайшего кратного grid_size
                x_grid = x // self.grid_size
                y_grid = y // self.grid_size
                
                # Отображаем текущие значения x и y на сетке
                if x_grid!= prev_x_grid or y_grid!= prev_y_grid:
                    self.canvas.create_rectangle(x_grid*self.grid_size, y_grid*self.grid_size, (x_grid+1)*self.grid_size, (y_grid+1)*self.grid_size, fill='lightgray')
                    self.window.update()
                prev_x_grid = x_grid
                prev_y_grid = y_grid
                self.canvas.create_rectangle(x, y, x, y)
                self.window.update()
                #self.window.after(50)
            else:
                self.canvas.create_rectangle(x, y, x, y)
            if e>=0:
                y+=sy
                e-=2*dx
            else:
                x+=sx
                e+=2*dy
            

    
    def draw_line_wu(self):
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]

        dx = x2 - x1
        dy = y2 - y1 
        x, y = x1, y1
        prev_x_grid=prev_y_grid=0
        # Определяем направление рисования
        if abs(dx) > abs(dy):
            steps = abs(dx)
        else:
            steps = abs(dy)

        # Вычисляем значения интенсивности
        if steps != 0:
            xi = dx / steps
            yi = dy / steps
        else:
            xi = yi = 0

        # Рисуем отрезок
        for _ in range(int(steps)):
            # Округляем координаты
            x_int = int(x)
            y_int = int(y)
            
            # Вычисляем дробную часть координат
            x_frac = x - x_int
            y_frac = y - y_int
            
            # Вычисляем яркость пикселя
            brightness = 1 - (x_frac + y_frac) / 2
            
            # Устанавливаем цвет пикселя
            color = "#%02x%02x%02x" % (int(brightness * 255), int(brightness * 255), int(brightness * 255))
            
            # Рисуем пиксель
            if self.debug_mode:
                
                # Округляем значения x и y до ближайшего кратного grid_size
                x_grid = x // self.grid_size
                y_grid = y // self.grid_size
                
                # Отображаем текущие значения x и y на сетке
                if x_grid!= prev_x_grid or y_grid!= prev_y_grid:
                    self.canvas.create_rectangle(x_grid*self.grid_size, y_grid*self.grid_size, (x_grid+1)*self.grid_size, (y_grid+1)*self.grid_size, fill=color)
                    self.window.update()
                prev_x_grid = x_grid
                prev_y_grid = y_grid
                self.canvas.create_rectangle(x_int, y_int, x_int, y_int, fill="black")
                if steps == abs(dx):
                    self.canvas.create_rectangle(x_int, y_int+1, x_int, y_int+1, fill=color, outline=color)
                else:
                    self.canvas.create_rectangle(x_int+1, y_int, x_int+1, y_int, fill=color, outline=color)
                self.window.update()
                #self.window.after(50)
            else:
                self.canvas.create_rectangle(x_int, y_int, x_int, y_int, fill="black")
                if steps == abs(dx):
                    self.canvas.create_rectangle(x_int, y_int+1, x_int, y_int+1, fill=color, outline=color)
                else:
                    self.canvas.create_rectangle(x_int+1, y_int, x_int+1, y_int, fill=color, outline=color)
            # Переходим к следующему пикселю
            x += xi
            y += yi
       

    def draw_circle(self):
        x_center, y_center = self.points[0]
        x2, y2 = self.points[1]
        radius = math.sqrt((x_center - x2) ** 2 + (y_center - y2) ** 2)
        x = 0
        y = radius
        delta = 2 - 2 * radius
        error = 0
        while y > 0:
            if self.debug_mode:
               
            
                self.canvas.create_rectangle(x_center + x, y_center + y, x_center + x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center + y, x_center - x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center + x, y_center - y, x_center + x, y_center - y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center - y, x_center - x, y_center - y, outline='black')
                self.window.update()
                self.window.after(10)
            else:
                self.canvas.create_rectangle(x_center + x, y_center + y, x_center + x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center + y, x_center - x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center + x, y_center - y, x_center + x, y_center - y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center - y, x_center - x, y_center - y, outline='black')
            if delta < 0:
                error = 2 * (delta + y) - 1
                if error <= 0:
                    x += 1
                    delta += 2 * x + 1
                    continue

           
            if delta > 0:
                error = 2 * (delta - x) - 1
                if error > 0:
                    y -= 1
                    delta += 1 - 2 * y
                    continue

            x += 1
            delta += 2 * (x - y) + 2
            y -= 1  

    def draw_ellipsis(self):
        x_center, y_center = self.points[0]
        x2, y2 = self.points[1]
        #radius = math.sqrt((x_center - x2) ** 2 + (y_center - y2) ** 2)
        a=abs(x2-x_center)
        b=abs(y2-y_center)
        x = 0
        y = b
        delta = a**2 + b**2 - 2*a**2*b
       
        error = 0

        while y > 0:
            if self.debug_mode:
               
            
                self.canvas.create_rectangle(x_center + x, y_center + y, x_center + x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center + y, x_center - x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center + x, y_center - y, x_center + x, y_center - y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center - y, x_center - x, y_center - y, outline='black')
                self.window.update()
                self.window.after(10)
            else:
                self.canvas.create_rectangle(x_center + x, y_center + y, x_center + x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center + y, x_center - x, y_center + y, outline='black')
                self.canvas.create_rectangle(x_center + x, y_center - y, x_center + x, y_center - y, outline='black')
                self.canvas.create_rectangle(x_center - x, y_center - y, x_center - x, y_center - y, outline='black')
            if delta < 0:
                error = 2 * (delta + a**2*y) - 1
                if error <= 0:
                    x += 1
                    delta += b**2*(2 * x + 1)
                    continue

           
            if delta > 0:
                error = 2 * (delta - b**2*x) - 1
                if error > 0:
                    y -= 1
                    delta += a**2*(1 - 2 * y)
                    continue

            x += 1
            delta += b**2*(2*x+1)+a**2*(1-2*y)
            y -= 1  


        
    def draw_giperbola(self):
        x0, y0 = self.points[0]
        x1, y1 = self.points[1]
        c = (x1*y1 - x0*y0)/(y0 - y1)
        k = y0*y1*(x1 - x0)/(y0 - y1)
        x = 0
        prev_x_grid=prev_y_grid=-1
        while x < 800:
            x = x + 0.02 
            y = k / (x + c)
            if self.debug_mode:
                
                # Округляем значения x и y до ближайшего кратного grid_size
                x_grid = x // self.grid_size
                y_grid = y // self.grid_size
                
                # Отображаем текущие значения x и y на сетке
                if x_grid!= prev_x_grid or y_grid!= prev_y_grid:
                    self.canvas.create_rectangle(x_grid*self.grid_size, y_grid*self.grid_size, (x_grid+1)*self.grid_size, (y_grid+1)*self.grid_size, fill='lightgray')
                    
                    self.window.update()
                prev_x_grid = x_grid
                prev_y_grid = y_grid
                self.canvas.create_rectangle(x, y, x, y, outline='black')
            else:
                self.canvas.create_rectangle(x, y, x, y, outline='black')
            
            
    def draw_parabola(self):
        x, y = self.points[0]
        x2, y2 = self.points[1]
        xc,yc=x,y
        self.canvas.create_oval(x-5, y-5, x+5, y+5, fill="red")  # Отметить клик мыши красным кружком
      
        #self.canvas.create_line(x, 0, x, self.height, fill="blue")  # Вертикальная линия через вершину
        #self.canvas.create_line(0, y, self.width, y, fill="blue")  # Горизонтальная линия через вершину
        
        # Генерация параболы
        a = (y2 - y) / ((x2 - x) ** 2)
        b = -2 * a * x
        c = y - a * (x ** 2) - b * x
        x = 0
        y = c
        prev_x_grid=prev_y_grid=-1
        while x <= 800:
            y = a * (x ** 2) + b * x + c
            if self.debug_mode:
                
                # Округляем значения x и y до ближайшего кратного grid_size
                x_grid = x // self.grid_size
                y_grid = y // self.grid_size
                
                # Отображаем текущие значения x и y на сетке
                if x_grid!= prev_x_grid or y_grid!= prev_y_grid:
                    self.canvas.create_rectangle(x_grid*self.grid_size, y_grid*self.grid_size, (x_grid+1)*self.grid_size, (y_grid+1)*self.grid_size, fill='lightgray')
                    self.window.update()
                    
                prev_x_grid = x_grid
                prev_y_grid = y_grid
                self.canvas.create_rectangle(x, y, x, y, outline='black')
                self.window.update()
                #self.window.after(1)
            else:
                self.canvas.create_rectangle(x, y, x, y, fill="black")
            x += 1
        self.canvas.create_oval(xc-5, yc-5, xc+5, yc+5, fill="red")

    def draw_curve_Bezie(self):
        x0, y0 = self.points[0]
        x1, y1 = self.points[1]
        x2, y2 = self.points[2]
        x3, y3 = self.points[3]
       
        #t_arr = np.array([pow(t,3), t*t, t, 1])
        arr=np.array([[-1,3,-3,1],
                      [3,-6,3,0],
                      [-3,3,0,0],
                      [1,0,0,0]
                    ])
        p_arr=np.array([[x0,y0],
                        [x1,y1],
                        [x2,y2],
                        [x3,y3]
                    ])
        p=np.dot(arr,p_arr)
        for t in range(0, 1001, 1):
            t = t / 1000.0
            t_arr = np.array([pow(t,3), t*t, t, 1])
            dot = np.dot(t_arr,p)
            self.canvas.create_rectangle(dot[0], dot[1], dot[0], dot[1], fill="black")

    # def calculate_derivatives(self):
    #     self.derivatives = []
    #     for i in range(2):
    #         dx = self.points[i+2][0] - self.points[i][0]
    #         dy = self.points[i+2][1] - self.points[i][1]
    #         self.derivatives.append((dx, dy))

    def draw_curve_Hermit(self):
        p0= self.points[0]
        p1 = self.points[1]
        v0 = [0, 600]
        v1 = [800, 0]
        
        for i in range(0, 1001):
            t = i / 1000.0
            x = (2*t**3 - 3*t**2 + 1) * self.points[0][0] + (-2*t**3 + 3*t**2) * self.points[1][0] + (t**3 - 2*t**2 + t) * v0[0] + (t**3 - t**2) * v1[0]
            y = (2*t**3 - 3*t**2 + 1) * self.points[0][1] + (-2*t**3 + 3*t**2) * self.points[1][1] + (t**3 - 2*t**2 + t) * v0[1] + (t**3 - t**2) * v1[1]
            self.canvas.create_rectangle(x, y, x, y)

    def draw_B_splain(self):
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]
        x3, y3 = self.points[2]
        x4, y4 = self.points[3]
       
        #t_arr = np.array([pow(t,3), t*t, t, 1])
        arr=np.array([[-1,3,-3,1],
                      [3,-6,3,0],
                      [-3,0,3,0],
                      [1,4,1,0]
                    ])
        p_arr=np.array([[x1,y1],
                        [x2,y2],
                        [x3,y3],
                        [x4,y4]
                    ])
        for i in range(3):
            
            p=np.dot(arr,p_arr)
            self.canvas.create_rectangle(x1, y1, x1+5, y1+5, fill="red")
            self.canvas.create_rectangle(x2, y2, x2+5, y2+5, fill="red")
            self.canvas.create_rectangle(x3, y3, x3+5, y3+5, fill="red")
            self.canvas.create_rectangle(x4, y4, x4+5, y4+5, fill="red")
            
            for t in range(0, 1001, 1):
                t = t / 1000.0
                t_arr = np.array([pow(t,3), t*t, t, 1])
                t_arr=1/6*t_arr
                dot = np.dot(t_arr,p)
                # x = (1-t)**3/6 * x1 + (3*t**3-6*t**2+4)/6 * x2 + (-3*t**3+3*t**2+3*t+1)/6 * x3 + t**3/6 * x4
                # y = (1-t)**3/6 * y1 + (3*t**3-6*t**2+4)/6 * y2 + (-3*t**3+3*t**2+3*t+1)/6 * y3 + t**3/6 * y4
                self.canvas.create_rectangle(dot[0], dot[1], dot[0], dot[1], fill="black")
            p_arr=np.roll(p_arr,-1,axis=0)


    def draw_grid(self):
        for x in range(0, self.window.winfo_screenwidth(), self.grid_size):
            self.canvas.create_line(x, 0, x, self.window.winfo_screenheight(), fill="gray")
            
        for y in range(0, self.window.winfo_screenheight(), self.grid_size):
            self.canvas.create_line(0, y, self.window.winfo_screenwidth(), y, fill="gray")
        
    def delete_grid(self):        
        self.debug_mode = False
        self.canvas.delete("all")
        self.canvas.create_line(0, 2, self.window.winfo_screenwidth(), 2, fill="gray")
        #self.cda = tk.Button(text="Отрезки" , state="disabled")
        #self.debug = tk.Button(text="Отладка", state="disabled")

    def run(self):
        self.window.mainloop()
        

editor = GraphicEditor(800,600)
editor.run()

# def calculate_derivatives(self):
#         self.derivatives = []
#         for i in range(2):
#             dx = self.points[i+2][0] - self.points[i][0]
#             dy = self.points[i+2][1] - self.points[i][1]
#             self.derivatives.append((dx, dy))

# def draw_curve_Hermit(self):
#     for i in range(0, 500):
#         t = i / 500.0
#         x = (2*t**3 - 3*t**2 + 1) * self.points[0][0] + (-2*t**3 + 3*t**2) * self.points[2][0] + (t**3 - 2*t**2 + t) * self.derivatives[0][0] + (t**3 - t**2) * self.derivatives[1][0]
#         y = (2*t**3 - 3*t**2 + 1) * self.points[0][1] + (-2*t**3 + 3*t**2) * self.points[2][1] + (t**3 - 2*t**2 + t) * self.derivatives[0][1] + (t**3 - t**2) * self.derivatives[1][1]
#         self.canvas.create_rectangle(x, y, x, y,)