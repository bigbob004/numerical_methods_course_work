import numpy as np
import matplotlib.pyplot as plt
from math import ceil

#Исходная функция, заданная в декартовых коор-ах
def f(x, y):
    return x + y


#Процедура для вычисления по полярным коор-ам функции в декартовых коор-ах
def polar_f(alpha, r):
    x = r * np.cos(alpha)
    y = r * np.sin(alpha)
    return f(x, y)


def test_cylindrical_func(alpha, r):
    return alpha + r


#Преобразование точки в декатовоую
def to_decart(point):
    phi, r = point
    return (r * np.cos(phi), r * np.sin(phi))


#Преобразование точки в полярную
def to_polar(point):
    x, y = point
    alp = np.arctan(np.abs(y) / np.abs(x)) if x != 0 else 0
    if x * y < 0:
        alp = np.pi - alp
    if y < 0:
        alp = np.pi + alp
    r = (x**2 + y**2)**0.5
    return (alp, r)


#Процедура вычисления билинейного многочлена
def bi_linear_polynom(val_00, val_01, val_10, val_11, point_0, point_1):
    tX_0, tY_0 = point_0
    tX_1, tY_1 = point_1
    def worker(x, y, X_0 = tX_0, X_1 = tX_1, Y_0 = tY_0, Y_1 = tY_1):
        return (val_00 * (X_1 - x) * (Y_1 - y) \
                + val_10 * (x - X_0) * (Y_1 - y) \
                + val_01 * (X_1 - x) * (y - Y_0) \
                + val_11 * (x - X_0) * (y - Y_0)) / ((X_1 - X_0) * (Y_1 - Y_0))
    return worker


#Определим класс элементарной области круга -
#элементарная область круга - сектор, полученного разбиениея,
#переводим такой сектор в полярные коор-ты - получаем прямоугольник

#Cам класс элементарной области
class Elementary_area:
    #Процедура инициализации элементарной области
    #left_up - верхняя левая точка сектора
    #right_down - правая нижняя точка сектора
    #area_func - функция(билинейный многочлен), с помощью которой можно вычислять зн-ие интерполированной функции в точке, находящейся в данной области
    def __init__(self, left_up, right_down, area_func):
        self.left_up = left_up
        self.right_down = right_down
        self.area_func = area_func

#Класс интерполяционного многочлена
class Interpolated_func:
    #Процедура инициализации интерполяционного многочлена
    #list_func_areas - элементарные области круга(см. выше об элементарной области круга)
    #dphi - элементарный полярный угол
    #dr - элементарный полярный радиус
    def __init__(self, list_func_areas):
        self.func_areas = list_func_areas


    #Процедура поиска элементарной области, в которой находится заданная точка.
    #(чтобы вычислять зн-ие интерполяционного многочлена в точке, нужно понять, а в какой элементарной области лежит эта точка)
    #Перебираем все области, которые мы получили с помощью разбиения.
    #Элементарную область ограничивают 2 радиуса и 2 угла.
    #Чтобы точка попала в область, нужно чтобы её полярные радиус находился в отрезке от меньшего и до большего полярного радиуса области(с углом аналогично)
    #Как только найдём такую область, то как результат вернём билинейный многочлен для данной области(да, Python в качестве результат процедуры может возвращать процедуру
    def search(self, phi, r):
        for line in self.func_areas:
            cur_area = line[0]
            if cur_area.right_down[1] <= r <= cur_area.left_up[1]:
                for area in line:
                    if area.right_down[0] <= phi <= area.left_up[0]:
                        return area.area_func


    #Процедура вычисления зн-ия интерполированной функции
    #Процедура принимает на вход саму точку и билиненый многочлен, с помощью которого можно вычислить зн-ие в данной точке(его мы нашли раннее в процедуре search
    def evaluate(self, point, is_decart=True):
        phi, r = point
        if is_decart:
            phi, r = to_polar(point)
        func = self.search(phi, r)
        return func(phi, r)


    #Процедура получения разбиения в полярных коор-ах
    #Делаем разбиение с помощью linspace(разбиение угла и радиуса(с работой linspace можно ознакомится в официальной документации)
    #С заданными шагами обойдём весь круг, обходить будем так:
    #Для данного радиуса переберём все углы(так мы обойдём весь круг радиуса i * dr). После перебора всех углов, увеличим радиус, чтобы перейти на следующий круг
    #Полученные точки будем складировать в список 2 на 2, где в пределах одной строки собарны точки с одним и тем же радиусом, соот-но, в пределах
    #одного столбца собраны точки с одним и тем же углом
def get_points(cnt_dphi, cnt_dr, R):
    points = []
    for cur_r in np.linspace(0, R, cnt_dr, endpoint=True):
        point_with_cur_r = []
        for cur_phi in np.linspace(0, 2 * np.pi, cnt_dphi, endpoint=True):
            point_with_cur_r.append((cur_phi, cur_r))
        points.append(point_with_cur_r)
    return points


#Формируем таблично-заданную функцию f, где f - интерполируемая функция
def get_func_values(f, points):
    values = []
    for list in points:
        buf = []
        for point in list:
            x, y = point
            buf.append(f(x, y))
        values.append(buf)
    return values


#Процедура получения элементарных областей и билинейного многочлена к соот-ей области
#Принцип хранения элементарных областей такой же, как и точек разбиения - список 2 на 2, где в пределах одной строки - области с одинаковым радиусом
#в пределах одного столбца - область с одним и тем же углом
def get_areas(cnt_dphi, cnt_dr, points, f_polar_values, bi_linear_polynom):
    areas = []
    #Проходимся по всем точкам разбиения(процедура get_points)
    for i in range(cnt_dr - 1):
        areas_with_phi = []
        for j in range(cnt_dphi - 1):
            #Получаем правую нижнюю и левую верхнюю точку разбиения(с помощью этих точек зададим область)
            right_down = points[i][j]
            left_up = points[i+1][j+1]
            #Получаем зн-ие интерполируемой функции в узловых точках
            val_00, val_01 = f_polar_values[i][j], f_polar_values[i+1][j]
            val_10, val_11 = f_polar_values[i][j+1], f_polar_values[i+1][j+1]
            #Строим билинейный многочлен по узловым точкам прямоугольника и зн-ий интерполируемой функции в узлах этого прямоугольника
            #(прямоугольник, потому что в полярных коор-ах элементарный сектор преобразуется в прямоугольник)
            func_area = bi_linear_polynom(val_00, val_01, val_10, val_11, (right_down[0], right_down[1]), (left_up[0], left_up[1]))
            areas_with_phi.append(Elementary_area(left_up, right_down, func_area))
        areas.append(areas_with_phi)
    return areas


#Перевод 2-мерного списка полярных точек в 2-ый список декартовых точек
def points_to_decart(points):
    return [[to_decart(point) for point in lst_points] for lst_points in points]


#Перевод 2-мерного списка декартовых точек в 2-ый список полярных точек
def points_to_polar(points):
    return [[to_polar(point) for point in lst_points] for lst_points in points]


#Построение графика поверхности
def my_plot_surf(ax_3d, points, f_values):
    x, y = [], []
    for lst_points in points:
        x_buf, y_buf = [], []
        for point in lst_points:
            x_buf.append(point[0])
            y_buf.append(point[1])
        x.append(x_buf)
        y.append(y_buf)

    ax_3d.set_xlabel('x')
    ax_3d.set_ylabel('y')
    ax_3d.set_zlabel('z')
    ax_3d.plot_surface(np.array(x), np.array(y), np.array(f_values), cmap='gnuplot')


#Построение поточечного 3-мерного графика
def my_plot_scat(ax_3d, points, f_values, clr = 'red', s = 3):
    x, y, z = [], [], []
    for lst_points in points:
        for point in lst_points:
            x.append(point[0])
            y.append(point[1])
    for lst_points in f_values:
        for point in lst_points:
            z.append(point)
    ax_3d.set_xlabel('x')
    ax_3d.set_ylabel('y')
    ax_3d.set_zlabel('z')
    ax_3d.scatter(x, y, z, s=s, c=clr)


#Таблица погрешностей
def errors(points, f_values, interpolated_values):
    max = 0
    lst_errors = []
    print("Точка (x, y)       | Зн-ие реальной функции       | Зн-ин интерпол. фун-ии       | Погрешность")
    for i in range(len(f_values)):
        buf_errors = []
        for j in range(len(f_values[i])):
            real = f_values[i][j]
            interpolated = interpolated_values[i][j]
            error = np.abs(real - interpolated)
            point = points[i][j]
            if error > max: max = error
            print("(%.4f, %.4f)                %.4f                %.4f                 %.4f" % (point[0], point[1], real, interpolated, error))
            buf_errors.append(error)
        lst_errors.append(buf_errors)
    print("Максимум погрешности: ", max)
    return lst_errors



#Таблица точек и зн-ий функции в точке
def print_table_func(points, f_values):
    print("Точка (x, y)       | Зн-ие реальной функции")
    for i in range(1, len(f_values)):
        for j in range(len(f_values[i]) - 1):
            real = f_values[i][j]
            point = points[i][j]
            print("(%.4f, %.4f)                %.10f" % (point[0], point[1], real))


#Построение графика интерполяции по узловым точкам(полярные точки)
def interpolated_func_plot(cnt_dphi, cnt_dr, R, interpolated_f):
    x = np.linspace(0, 2 * np.pi, cnt_dphi, endpoint=True)
    y = np.linspace(0, R, cnt_dr, endpoint=True)

    points = []
    values = []

    for y_point in y:
        buf_point = []
        buf_val = []
        buf_f_val = []
        for x_point in x:
            buf_val.append(interpolated_f.evaluate((x_point, y_point), False))
            buf_point.append((x_point, y_point))
        points.append(buf_point)
        values.append(buf_val)

    decart_points = points_to_decart(points)

    fig = plt.figure(figsize=(15, 7))
    ax_3d = fig.add_subplot(projection='3d')
    ax_3d.set_xlabel('x')
    ax_3d.set_ylabel('y')
    ax_3d.set_zlabel('z')
    my_plot_surf(ax_3d, decart_points, values)
    plt.show()


#Зададим круг радиуса 10
R = 10
cnt_dr, cnt_dphi = 32, 15

#Исходный график
x = y = np.linspace(-10, 10, 10, endpoint=True)
x_grid, y_grid = np.meshgrid(x, y)
z_grid = f(x_grid, y_grid)

# #Функция для объяснения тестовых примеров
# x = np.linspace(0, 2 * np.pi, 10, endpoint=True)
# y = np.linspace(0, R, 10, endpoint=True)
# x_grid, y_grid = np.meshgrid(x, y)
# z_grid = test_cylindrical_func(x_grid, y_grid)

fig = plt.figure(figsize=(30, 20))
ax_3d = fig.add_subplot(projection='3d')
ax_3d.plot_surface(x_grid, y_grid, z_grid)
ax_3d.set_xlabel('x')
ax_3d.set_ylabel('y')
ax_3d.set_zlabel('z')
plt.show()


polar_points = get_points(cnt_dphi, cnt_dr, R)
decart_points = points_to_decart(polar_points)

x, y = [], []
for line in decart_points:
    for point in line:
        x.append(point[0])
        y.append(point[1])

plt.rcParams.update({'font.size': 15})
plt.figure(figsize=(5, 5))
plt.title("Точки разбиения")
plt.scatter(x, y, s=1)
#plt.show()

func_values = get_func_values(polar_f, polar_points)
print_table_func(decart_points, func_values)

func_areas = get_areas(cnt_dphi, cnt_dr, polar_points, func_values, bi_linear_polynom)


interpolated = Interpolated_func(func_areas)

#Построим график интерполяции по узлам интерполяции
interpolated_func_plot(100, 100, R, interpolated)

points = []
values = []
f_values = []

x = np.linspace(-R, R, cnt_dphi, endpoint=True)
y = np.linspace(-R, R, cnt_dr, endpoint=True)

for y_point in y:
    buf_point = []
    buf_val = []
    buf_f_val = []
    for x_point in x:
        if x_point ** 2 + y_point ** 2 <= R ** 2:
            buf_val.append(interpolated.evaluate((x_point, y_point)))
            buf_point.append((x_point, y_point))
            buf_f_val.append(f(x_point, y_point))
    if len(buf_point) != 0 and len(buf_val) != 0:
        points.append(buf_point)
        values.append(buf_val)
        f_values.append(buf_f_val)

decart_points = points_to_decart(points)

lst_errors = errors(decart_points, f_values, values)


fig = plt.figure(figsize=(15, 7))
ax_3d = fig.add_subplot(projection='3d')
my_plot_scat(ax_3d, decart_points, lst_errors)
plt.show()


