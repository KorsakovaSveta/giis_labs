import tkinter as tk
import numpy as np
import math, cmath
import numpy as np
from scipy.spatial import Delaunay, Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
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
        self.mainmenu.add_cascade(label="Куб", command=self.activate_canvas_draw_cube)
        self.mainmenu.add_cascade(label="Перемещение", command=self.activate_canvas_translate_cube)
        self.mainmenu.add_cascade(label="Поворот", command=self.activate_canvas_rotate_cube)
        self.mainmenu.add_cascade(label="Скалирование", command=self.activate_canvas_scale_cube)
        self.mainmenu.add_cascade(label="Отображение", command=self.activate_canvas_reflect_cube)
        self.mainmenu.add_cascade(label="Отладка", command=self.debug_mode_toggle)
        self.mainmenu.add_cascade(label="Остановить отладку", command=self.delete_grid)
        self.mainmenu.add_cascade(label="Корректировка", command=self.activate_canvas_editing)
        self.mainmenu.add_cascade(label="Остановить корректировку", command=self.end_editing)
        self.mainmenu.add_cascade(label="Очистить экран", command=self.clear_window)
        self.mainmenu.add_cascade(label="Полигон", command=self.activate_canvas_polygon)
        self.mainmenu.add_cascade(label="Грэхем", command=self.on_mouse_click_graham)
        self.mainmenu.add_cascade(label="Джарвис", command=self.on_mouse_click_jarvis)
        self.mainmenu.add_cascade(label="Пересечение", command=self.activate_canvas_draw_intersection_line)
        self.mainmenu.add_cascade(label="Принадлежность точки", command=self.activate_canvas_points_inside_polygon)
        self.mainmenu.add_cascade(label="Растровая развёртка", command=self.scanline_fill)
        self.mainmenu.add_cascade(label="Активные рёбра", command=self.scanline_active_fill)
        self.mainmenu.add_cascade(label="Затравка", command=self.activate_canvas_flood_fill)
        self.mainmenu.add_cascade(label="Построчная затравка", command=self.activate_canvas_string_flood_fill)
        self.mainmenu.add_cascade(label="Триангуляция", command=self.triangulate)
        self.mainmenu.add_cascade(label="Диаграмма Вороного", command=self.voronoi)
        self.points_label = tk.Label(self.window, text="Координаты точек (x1 y1, x2 y2, ...):")
        self.points_label.pack()
        self.points_entry = tk.Entry(self.window)
        self.points_entry.pack()
        self.center = np.array([self.canvas.winfo_width() / 2, self.canvas.winfo_height() / 2, 0, 1])
        self.new_point = None
        self.changing_and_new_point = []
        self.dragging = False
        self.canvas.pack()
        self.Bezie_editing = False
        self.Hermit_editing = False
        self.B_splain_editing = False
        self.points = []
        self.derivatives = []
        self.debug_mode = False
        self.grid_size = 10
        self.selected_point = None
        self.points_for_change = []
        self.all_points=[]
        self.for_all_points = []
        self.translateCube = False
        self.rotateCube = False
        self.cube_size = 0
        self.rotation_angle = 0
        self.scaleCube = False
        self.reflectCube = False
        self.perspectiveCube = False
        self.exit = False
        self.current_polygon = None
        self.graham = False
        self.jarvis = False
        self.graham_points = []
        self.jarvis_points = []
        self.line = []
        self.points_for_intersection =[]
        self.scanline_fill_points, self.edges, self.rasterized_points,self.flood_point = [], [], [], []

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
        self.Bezie_editing = True
        self.canvas.bind('<Button-1>', self.on_mouse_click_curve_Bezie)
        
    def activate_canvas_curve_Hermit(self):
        self.Hermit_editing = True
        self.canvas.bind('<Button-1>', self.on_mouse_click_curve_Hermit)
        
    def activate_canvas_B_splain(self):
        self.B_splain_editing = True
        self.canvas.bind('<Button-1>', self.on_mouse_click_B_splain)

    def activate_canvas_draw_cube(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_draw_cube)

    def activate_canvas_translate_cube(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_translate_cube)

    def activate_canvas_rotate_cube(self):
        self.window.bind('<KeyPress>', self.on_keybord_click_rotate_cube)

    def activate_canvas_scale_cube(self):
        self.window.bind('<KeyPress>', self.on_keybord_click_scale_cube)
    
    def activate_canvas_reflect_cube(self):
        self.window.bind('<KeyPress>', self.on_keybord_click_reflect_cube)

    def activate_canvas_polygon(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_polygon)
    
    def activate_canvas_draw_intersection_line(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_draw_intersection_line)
   
    def activate_canvas_points_inside_polygon(self):
        self.canvas.bind('<Button-1>', self.point_inside_polygon)

    def activate_canvas_editing(self):
        self.canvas.bind('<Button-1>', self.start_editing)

    def activate_canvas_flood_fill(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_flood_fill)

    def activate_canvas_string_flood_fill(self):
        self.canvas.bind('<Button-1>', self.on_mouse_click_string_flood_fill)

    def debug_mode_toggle(self):
        self.canvas.delete("all")
        self.debug_mode = not self.debug_mode
        self.draw_grid()

    def clear_window(self):
        self.canvas.delete("all")
        self.all_points.clear()

    def start_editing(self, event):
        self.dragging = True  
        self.changing_and_new_point.append((event.x, event.y)) 
        for point in self.points_for_change:
            self.canvas.create_rectangle(point[0]-4, point[1]-4, point[0]+4, point[1]+4, fill="red")
        if len(self.changing_and_new_point) == 2:
            
            changing_point = self.changing_and_new_point[0]
            for point in self.points_for_change:
                if point[0]-3<=changing_point[0]<=point[0]+3 or\
                    point[1]-3<=changing_point[1]<=point[1]+3:  
                    new_point = self.changing_and_new_point[1]
                    index_changed_point = self.points_for_change.index(point)
                    self.points_for_change[index_changed_point] = new_point
                    if self.Bezie_editing:
                        
                        self.canvas.delete("all")
                        for all in self.all_points:
                            self.points=all
                            self.draw_curve_Bezie()
                        self.points = self.points_for_change
                        self.draw_curve_Bezie()
                        self.changing_and_new_point.clear()

                    if self.Hermit_editing:
                        
                        self.canvas.delete("all")
                        for all in self.all_points:
                            self.points=all
                            self.draw_curve_Hermit()
                        self.points = self.points_for_change
                        self.draw_curve_Hermit()
                        self.changing_and_new_point.clear()
                    if self.B_splain_editing:
                        
                        self.canvas.delete("all")
                        for all in self.all_points:
                            self.points=all
                            self.draw_B_splain()
                        self.points = self.points_for_change
                        self.draw_B_splain()
                        self.changing_and_new_point.clear()

    def end_editing(self):
        self.dragging  = False
        self.Hermit_editing= self.Bezie_editing= self.B_splain_editing=False
        for i in range(len(self.points)):
            self.for_all_points.append(self.points[i])   
        self.all_points.append(self.for_all_points)  
        self.for_all_points = []
        self.points.clear()
        self.changing_and_new_point.clear()    

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
            if self.dragging == False:
                self.points_for_change = self.points
                self.points = []
            
    def on_mouse_click_curve_Hermit(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) ==2:   
      
            self.draw_curve_Hermit()   
            if self.dragging == False:
                self.points_for_change = self.points
                self.points = []        
            

    def on_mouse_click_B_splain(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) ==4:         
            self.draw_B_splain()           
            if self.dragging == False:
                self.points_for_change = self.points
                self.points = []
    
    def on_mouse_click_draw_cube(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) ==2:         
            self.draw_cube() 
            self.points_rotate = self.points     
            self.points = []
    
    def on_mouse_click_translate_cube(self, event):  
        self.points.append((event.x, event.y))     
        self.translate_cube()   
        self.points_rotate[0] = self.points[0]        
        self.points = []

    def on_keybord_click_rotate_cube(self, event):
        if event.keysym == "Left":
            self.rotation_angle += 10
        elif event.keysym == "Right":
            self.rotation_angle -= 10
        self.rotate_cube() 

    def on_keybord_click_scale_cube(self, event):
        if event.keysym == "Up":
            self.cube_size += 10
        elif event.keysym == "Down":
            self.cube_size -= 10     
        self.scale_cube() 

    def on_keybord_click_reflect_cube(self, event):
        if event.keysym == "r":
            self.reflect_cube()
        elif event.keysym == "e":
            self.exit = True
            self.points = self.points_rotate
            self.canvas.delete("all")            
            self.draw_cube()
            self.points = []
            self.exit = False
    
    def on_mouse_click_graham(self):
        self.canvas.delete(*self.line)
        self.graham = True
        self.draw_convex_hull_graham(self.graham_points)
        self.graham = False

    def on_mouse_click_jarvis(self):
        self.canvas.delete(*self.line)
        self.jarvis = True
        
        self.draw_convex_hull_jarvis(self.jarvis_points)
        self.jarvis = False

    def on_mouse_click_draw_intersection_line(self, event):
        self.points.append((event.x, event.y))
        if len(self.points) == 2:
            self.canvas.create_line(self.points[0], self.points[1], width=2)
            self.line_segment_intersection()
            self.points = []


    def on_mouse_click_polygon(self, event):
        self.graham_points = []
        self.jarvis_points = []
        self.points.append((event.x, event.y))
        if len(self.points) == 5:            
            self.canvas.delete(self.current_polygon)
            self.current_polygon = None
            if self.is_convex(self.points):
                self.canvas.create_polygon(self.points, fill="", outline="green", tags="polygon")
                self.calculate_normals(self.points)
            else:
                self.canvas.create_polygon(self.points, fill="", outline="red", tags="polygon")
                self.calculate_normals(self.points)
            for point in self.points:
                self.canvas.create_oval(point[0]-3, point[1]-3, point[0]+3, point[1]+3, fill="black")
            self.graham_points = self.points
            self.jarvis_points = self.points
            self.points_for_intersection = self.points
            self.scanline_fill_points = self.points
            self.points = []   

    def on_mouse_click_flood_fill(self, event):  
        self.flood_point.append((event.x, event.y))     
        self.flood_fill()         
        self.flood_point = []

    def on_mouse_click_string_flood_fill(self, event):  
        self.flood_point.append((event.x, event.y))     
        self.string_flood_fill()         
        self.flood_point = []

    def is_convex(self, points):
        cross_product = []
        for i in range(5):
            p1 = points[i]
            p2 = points[(i+1)%5]
            p3 = points[(i+2)%5]            
            x1, y1 = p1
            x2, y2 = p2
            x3, y3 = p3            
            cross_product.append((x2-x1)*(y3-y2) - (y2-y1)*(x3-x2))        
        signs = [1 if cp > 0 else -1 if cp < 0 else 0 for cp in cross_product]        
        if all(sign > 0 for sign in signs) or all(sign < 0 for sign in signs):
            return True
        else:
            return False
    
    def calculate_normals(self, points):
        mid_points = []
        for i in range(5):
            x1, y1 = points[i]
            x2, y2 = points[(i+1)%5]            
            mid_x = (x1 + x2) / 2
            mid_y = (y1 + y2) / 2            
            mid_points.append((mid_x, mid_y))            
        normals = []
        for i in range(5):
            x1, y1 = mid_points[i]
            x2, y2 = mid_points[(i+1)%5]            
            normal_x = y2 - y1
            normal_y = x1 - x2            
            normals.append((normal_x, normal_y))
            self.line.append(self.canvas.create_line(x1, y1, x1 + normal_x, y1 + normal_y, arrow=tk.LAST))            
        return normals
    
    def draw_convex_hull_graham(self, points):        
        hull = []
        start = min(points, key=lambda x: (x[1], x[0]))

        # Sort points by polar angle with the start point
        points.sort(key=lambda p: math.atan2(p[1]-start[1], p[0]-start[0]))

        # Initialize the hull with the start point
        hull.append(start)

        # Add points to the hull
        for point in points[1:]:
            while len(hull) > 1 and not self.is_left_turn(hull[-2], hull[-1], point):
                hull.pop()
            hull.append(point)
        self.canvas.create_polygon(hull, fill="white", outline="blue")
        for point in points:
            self.canvas.create_oval(point[0]-3, point[1]-3, point[0]+3, point[1]+3, fill="black")

    def is_left_turn(self, p1, p2, p3):
        return (p2[0]-p1[0]) * (p3[1]-p2[1]) - (p2[1]-p1[1]) * (p3[0]-p2[0]) > 0

    def draw_convex_hull_jarvis(self, points):
        self.hull_points = []
        p = min(points, key=lambda x: x[1])
        self.hull_points.append(p)
        while True:
            q = None
            for r in points:
                if r != p:
                    if q is None or self.orientation(p, q, r) > 0:
                        q = r
            p = q
            if p == self.hull_points[0]:
                break
            self.hull_points.append(p)
        self.canvas.create_polygon(self.hull_points, fill="white", outline="orange")
        for point in points:
            self.canvas.create_oval(point[0]-2, point[1]-2, point[0]+2, point[1]+2, fill="black")
            
    def orientation(self, p, q, r):
        return (q[1]-p[1]) * (r[0]-q[0]) - (q[0]-p[0]) * (r[1]-q[1])
        
    def line_segment_intersection(self):
        for i in range(len(self.points_for_intersection)):
            side_start = self.points_for_intersection[i]
            side_end = self.points_for_intersection[(i + 1) % len(self.points_for_intersection)]
            x1, y1 = self.points[0]
            x2, y2 = self.points[1]
            x3, y3 = side_start
            x4, y4 = side_end
            denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
            if denominator != 0:
                t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denominator
                u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denominator
                if 0 <= t <= 1 and 0 <= u <= 1:
                    intersection_x = x1 + t * (x2 - x1)
                    intersection_y = y1 + t * (y2 - y1)
                    self.canvas.create_oval(intersection_x-3, intersection_y-3, intersection_x+3, intersection_y+3, fill="red")
                   
    def point_inside_polygon(self, event):
        x, y = event.x, event.y
        overlapping_objects = self.canvas.find_overlapping(x, y, x, y)
        if len(overlapping_objects) > 0:
            self.canvas.create_oval(x-2, y-2, x+2, y+2, fill="blue")
        else:
            self.canvas.create_oval(x-2, y-2, x+2, y+2, fill="black")
        self.points=[]

    def scanline_fill(self):
        self.canvas.delete("all")
        self.canvas.create_polygon(self.scanline_fill_points, fill="", outline="black", tags="polygon")
                
        for point in self.scanline_fill_points:
            self.canvas.create_oval(point[0]-3, point[1]-3, point[0]+3, point[1]+3, fill="black")
        edges = []
        for i in range(5):
            x1, y1 = self.scanline_fill_points[i]
            x2, y2 = self.scanline_fill_points[(i + 1) % 5]
            edges.append((x1, y1, x2, y2))
        ymin = min(point[1] for point in self.scanline_fill_points)
        ymax = max(point[1] for point in self.scanline_fill_points)

        active_edges = []
        for y in range(ymin, ymax+1):
            # if y==ymin:
            #     x = [point[0] for point in self.scanline_fill_points if point[1]==y]
            #     prev_x_start = prev_x_end = x[0]
            #     y_prev = y
            for edge in edges:
                x1, y1, x2, y2 = edge
                if y1 < y <= y2 or y2 < y <= y1:
                    if y2 != y1:
                        x_intersection = int(x1 + (y - y1) * (x2 - x1) / (y2 - y1))
                        active_edges.append(x_intersection)

            active_edges.sort()
            for i in range(0, len(active_edges), 2):
                self.canvas.create_line(active_edges[i], y, active_edges[i+1], y, fill="blue")
                # if y == ymax:
                #     self.canvas.create_polygon(prev_x_start, y_prev, prev_x_end, y_prev, active_edges[i+1], ymax, active_edges[i], ymax, fill="blue")
                # else:
                #     self.canvas.create_polygon(prev_x_start, y_prev, prev_x_end, y_prev, active_edges[i+1], y, active_edges[i], y, fill="blue")
                # prev_x_start = active_edges[i]
                # prev_x_end = active_edges[i+1]
                # y_prev = y
                
            active_edges = [x for x in active_edges if x > active_edges[-1]]
        
    def scanline_active_fill(self):
        self.canvas.delete("all")
        self.canvas.create_polygon(self.scanline_fill_points, fill="", outline="black", tags="polygon")
                
        for point in self.scanline_fill_points:
            self.canvas.create_oval(point[0]-3, point[1]-3, point[0]+3, point[1]+3, fill="black")
        polygon_points = self.scanline_fill_points
        min_y = min(point[1] for point in polygon_points)
        max_y = max(point[1] for point in polygon_points)
        edge_table = [[] for _ in range(max_y - min_y + 1)]
        
        for i in range(len(polygon_points)):
            x1, y1 = polygon_points[i]
            x2, y2 = polygon_points[(i + 1) % len(polygon_points)]
            
            if y1 != y2:
                if y1 > y2:
                    y1, y2 = y2, y1
                    x1, x2 = x2, x1
                
                m = (x2 - x1) / (y2 - y1)
                edge_table[y1 - min_y].append((y2, x1, m))
                
        active_edges = []
        
        for y in range(min_y, max_y + 1):
            for edge in edge_table[y - min_y]:
                active_edges.append(edge)
                
            active_edges.sort()
            
            for i in range(len(active_edges) - 1, -1, -1):
                if y == active_edges[i][0]:
                    active_edges.pop(i)
                else:
                    active_edges[i] = (active_edges[i][0], active_edges[i][1] + active_edges[i][2], active_edges[i][2])
                    
            for i in range(0, len(active_edges), 2):
                x1 = int(active_edges[i][1])
                x2 = int(active_edges[i + 1][1])
                self.canvas.create_line(x1, y, x2, y, fill="green")

    def flood_fill(self):
        self.canvas.delete("all")
        self.canvas.create_polygon(self.scanline_fill_points, fill="", outline="black", tags="polygon")
        stack = []
        visited = set()
        image = np.zeros((self.canvas.winfo_height(), self.canvas.winfo_width()))

        # Помещаем первую вершину в стек
        stack.append(self.flood_point[0])

        # Пока стек не пуст
        while stack:
            x, y = stack.pop()

            # Закрашиваем пиксель и помечаем его как посещенный
            self.canvas.create_rectangle(x, y, x, y, fill="purple", outline="")
            

            # Проверяем соседние пиксели
            for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                nx, ny = x + dx, y + dy

                # Если пиксель не выходит за границы и не посещен
                
                if self.canvas.itemcget(self.canvas.find_closest(nx, ny), 'fill') != "purple" and self.canvas.itemcget(self.canvas.find_closest(nx, ny), 'fill') != "black" and (nx, ny) not in visited:
                    stack.append((nx, ny))

                    # Помечаем пиксель как посещенный
                    visited.add((nx, ny))

    def get_pixel_color(self, x, y):
            item_id = self.canvas.find_closest(x, y)[0]
            color = self.canvas.itemcget(item_id, "fill")
            return color

    def string_flood_fill(self):
        self.canvas.delete("all")
        self.canvas.create_polygon(self.scanline_fill_points, fill="", outline="black", tags="polygon")
        # Создаем стек и помещаем затравочный пиксел в него
       
        min_y = min(point[1] for point in self.scanline_fill_points)
        max_y = max(point[1] for point in self.scanline_fill_points)
        
        for y in range(min_y, max_y+1):
            intersections = []
            for i in range(len(self.scanline_fill_points)):
                x1, y1 = self.scanline_fill_points[i]
                x2, y2 = self.scanline_fill_points[(i+1) % len(self.scanline_fill_points)]
                if (y1 <= y < y2) or (y2 <= y < y1):
                    x = int(x1 + (y - y1) * (x2 - x1) / (y2 - y1))
                    intersections.append(x)
            
            intersections.sort()
            for i in range(0, len(intersections), 2):
                self.canvas.create_line(intersections[i], y, intersections[i+1], y, fill="red")

    
    def triangulate(self):
        # Получение координат точек из поля ввода
        points_str = self.points_entry.get()
        points_list = points_str.split(",")
        points = np.array([tuple(map(float, point.split())) for point in points_list])
        
        # Выполнение триангуляции Делоне
        triangulation = Delaunay(points)
        
        # Визуализация триангуляции Делоне
        plt.triplot(points[:,0], points[:,1], triangulation.simplices)
        plt.plot(points[:,0], points[:,1], 'o')
        plt.title('Триангуляция Делоне')
        plt.show()

    def voronoi(self):
        # Получение координат точек из поля ввода
        points_str = self.points_entry.get()
        points_list = points_str.split(",")
        points = np.array([tuple(map(float, point.split())) for point in points_list])
        
        # Выполнение построения диаграммы Вороного
        vor = Voronoi(points)
    
        # Визуализация диаграммы Вороного
        voronoi_plot_2d(vor)
        plt.plot(points[:,0], points[:,1], 'o')
        plt.title('Диаграмма Вороного')
        plt.show()
                                
    def draw_line_cda(self):
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
        for point in self.points:
            self.canvas.create_rectangle(point[0]-4, point[1]-4, point[0]+4, point[1]+4, fill="red")
        x0, y0 = self.points[0]
        x1, y1 = self.points[1]
        x2, y2 = self.points[2]
        x3, y3 = self.points[3]
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


    def draw_curve_Hermit(self):
        for point in self.points:
            self.canvas.create_rectangle(point[0]-4, point[1]-4, point[0]+4, point[1]+4, fill="red")
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
        for point in self.points:
            self.canvas.create_rectangle(point[0]-4, point[1]-4, point[0]+4, point[1]+4, fill="red")
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]
        x3, y3 = self.points[2]
        x4, y4 = self.points[3]
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
            for t in range(0, 1001, 1):
                t = t / 1000.0
                t_arr = np.array([pow(t,3), t*t, t, 1])
                t_arr=1/6*t_arr
                dot = np.dot(t_arr,p)
       
                self.canvas.create_rectangle(dot[0], dot[1], dot[0], dot[1], fill="black")
            p_arr=np.roll(p_arr,-1,axis=0)

    def draw_cube(self):
        if self.translateCube == False and\
            self.rotateCube == False and\
                self.scaleCube == False and\
                    self.reflectCube == False:
            x, y = self.points[0]
            if self.exit == False:
                self.cube_size = abs(self.points[1][1]-self.points[0][1])
            x2, y2 = self.points[1]
        elif self.translateCube == True:
            x, y = self.points[0]
        elif self.rotateCube == True:
            x,y = self.points_rotate[0]
        elif self.scaleCube == True: 
            x,y = self.points_rotate[0]
        if self.reflectCube == True:
            x,y = self.points_rotate[0]
            vertices = [
            (x - self.cube_size, y - self.cube_size),
            (x + self.cube_size, y - self.cube_size),
            (x + self.cube_size, y + self.cube_size),
            (x - self.cube_size, y + self.cube_size),
            (x - self.cube_size + 50, y - self.cube_size - 50),
            (x + self.cube_size + 50, y - self.cube_size - 50),
            (x + self.cube_size + 50, y + self.cube_size - 50),
            (x - self.cube_size + 50, y + self.cube_size - 50)
        ]    
        else:
        # Вершины куба
            vertices = [
            (x - self.cube_size, y - self.cube_size),
            (x + self.cube_size, y - self.cube_size),
            (x + self.cube_size, y + self.cube_size),
            (x - self.cube_size, y + self.cube_size),
            (x - self.cube_size - 50, y - self.cube_size - 50),
            (x + self.cube_size - 50, y - self.cube_size - 50),
            (x + self.cube_size - 50, y + self.cube_size - 50),
            (x - self.cube_size - 50, y + self.cube_size - 50)
        ]
    # Задаем грани куба
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)
        ]
        if self.rotateCube == True:
            rotation_matrix = [
                [math.cos(math.radians(self.rotation_angle)), -math.sin(math.radians(self.rotation_angle))],
                [math.sin(math.radians(self.rotation_angle)), math.cos(math.radians(self.rotation_angle))]
            ]
        # Поворот вершин куба
            for i in range(len(vertices)):
                x = vertices[i][0] - self.points_rotate[0][0]
                y = vertices[i][1] - self.points_rotate[0][1]
                new_x = x * rotation_matrix[0][0] + y * rotation_matrix[0][1]
                new_y = x * rotation_matrix[1][0] + y * rotation_matrix[1][1]
                vertices[i] = (new_x + self.points_rotate[0][0], new_y + self.points_rotate[0][1])
        # Рисуем куб 
        if self.perspectiveCube == False:
            for edge in edges:
                x1, y1 = vertices[edge[0]]
                x2, y2 = vertices[edge[1]]
                self.canvas.create_line(x1, y1, x2, y2)

    def translate_cube(self):
        self.translateCube = True
        self.canvas.delete("all")
        self.draw_cube()
        self.translateCube = False

    def rotate_cube(self):
        self.rotateCube = True
        self.canvas.delete("all")
        self.draw_cube()
        self.rotateCube = False

    def scale_cube(self):
        self.scaleCube = True
        self.canvas.delete("all")
        self.draw_cube()
        self.scaleCube = False

    def reflect_cube(self):
        self.reflectCube = True
        self.canvas.delete("all")
        self.draw_cube()
        self.reflectCube = False

    def draw_grid(self):
        for x in range(0, self.window.winfo_screenwidth(), self.grid_size):
            self.canvas.create_line(x, 0, x, self.window.winfo_screenheight(), fill="gray")
            
        for y in range(0, self.window.winfo_screenheight(), self.grid_size):
            self.canvas.create_line(0, y, self.window.winfo_screenwidth(), y, fill="gray")
        
    def delete_grid(self):        
        self.debug_mode = False
        self.canvas.delete("all")
        self.canvas.create_line(0, 2, self.window.winfo_screenwidth(), 2, fill="gray")
  
    def run(self):
        self.window.mainloop()
     
editor = GraphicEditor(800,600)
editor.run()
