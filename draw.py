from display import *
from matrix import *
from gmath import *

def draw_scanline(x0, z0, x1, z1, y, screen, zbuffer, color):
    if x0 > x1:
        tx = x0
        tz = z0
        x0 = x1
        z0 = z1
        x1 = tx
        z1 = tz

    x = x0
    z = z0
    delta_z = (z1 - z0) / (x1 - x0 + 1) if (x1 - x0 + 1) != 0 else 0

    while x <= x1:
        plot(screen, zbuffer, color, x, y, z)
        x+= 1
        z+= delta_z

def scanline_convert(polygons, i, screen, zbuffer, color):
    flip = False
    BOT = 0
    TOP = 2
    MID = 1

    points = [ (polygons[i][0], polygons[i][1], polygons[i][2]),
               (polygons[i+1][0], polygons[i+1][1], polygons[i+1][2]),
               (polygons[i+2][0], polygons[i+2][1], polygons[i+2][2]) ]

    points.sort(key = lambda x: x[1])
    x0 = points[BOT][0]
    z0 = points[BOT][2]
    x1 = points[BOT][0]
    z1 = points[BOT][2]
    y = int(points[BOT][1])

    distance0 = int(points[TOP][1]) - y * 1.0 + 1
    distance1 = int(points[MID][1]) - y * 1.0 + 1
    distance2 = int(points[TOP][1]) - int(points[MID][1]) * 1.0 + 1

    dx0 = (points[TOP][0] - points[BOT][0]) / distance0 if distance0 != 0 else 0
    dz0 = (points[TOP][2] - points[BOT][2]) / distance0 if distance0 != 0 else 0
    dx1 = (points[MID][0] - points[BOT][0]) / distance1 if distance1 != 0 else 0
    dz1 = (points[MID][2] - points[BOT][2]) / distance1 if distance1 != 0 else 0

    while y <= int(points[TOP][1]):
        if ( not flip and y >= int(points[MID][1])):
            flip = True

            dx1 = (points[TOP][0] - points[MID][0]) / distance2 if distance2 != 0 else 0
            dz1 = (points[TOP][2] - points[MID][2]) / distance2 if distance2 != 0 else 0
            x1 = points[MID][0]
            z1 = points[MID][2]

        draw_scanline(int(x0), z0, int(x1), z1, y, screen, zbuffer, color)
        x0+= dx0
        z0+= dz0
        x1+= dx1
        z1+= dz1
        y+= 1

def add_mesh(polygons, filename):
    file = open(filename + '.obj', 'r')
    lines = file.read().split("\n")
    points = ["placeholder"]
    for line in lines:
        x = line.split()
        if (len(x) == 0):
            continue
        if (x[0] == "f"):
            vertices = []
            for y in x[1:]:
                f = y.split("/")
                vertices.append(int(f[0]))
            p_0 = points[vertices[0]]
            p_1 = points[vertices[1]]
            p_2 = points[vertices[2]]
            if (len(vertices) == 4):
                p_3 = points[vertices[3]]
                add_polygon(polygons, p_0[0], p_0[1], p_0[2], p_1[0], p_1[1], p_1[2], p_2[0], p_2[1], p_2[2])
                add_polygon(polygons, p_0[0], p_0[1], p_0[2], p_2[0], p_2[1], p_2[2], p_3[0], p_3[1], p_3[2])
            if (len(vertices) == 3):
                add_polygon(polygons, p_0[0], p_0[1], p_0[2], p_1[0], p_1[1], p_1[2], p_2[0], p_2[1], p_2[2])
        if (x[0] == "v"):
            coords = [SCALING * float(coord) for coord in x[1:4]]
            points.append(coords)

def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0)
    add_point(polygons, x1, y1, z1)
    add_point(polygons, x2, y2, z2)

def draw_scanline_phong(x0, z0, x1, z1, y, screen, zbuffer, view, ambient, light, symbols, reflect, l_0, l_1):
    if x0 > x1:
        tx = x0
        tz = z0
        x0 = x1
        z0 = z1
        x1 = tx
        z1 = tz
        temp = l_0
        l_0 = l_1
        l_1 = temp

    x = x0
    z = z0
    dz = (z1 - z0) / (x1 - x0 + 1) if (x1 - x0 + 1) != 0 else 0

    distance = x1 - x0 + 1
    x_color = l_0[0]
    y_color = l_0[1]
    z_color = l_0[2]
    if(distance == 0):
        dx_color = 0
        dy_color = 0
        dz_color = 0
    else:
        dx_color = (l_1[0] - l_0[0]) / distance
        dy_color = (l_1[1] - l_0[1]) / distance
        dz_color = (l_1[2] - l_0[2]) / distance

    while x <= x1:
        pcolor = [ x_color, y_color, z_color]
        color = get_lighting(pcolor,view,ambient,light,symbols,reflect)
        plot(screen, zbuffer, color, x, y, z)
        x+= 1
        z+= dz
        x_color += dx_color
        y_color += dy_color
        z_color += dz_color

def scanline_phong(polygons, i, screen, zbuffer, colors, view, ambient, light, symbols, reflect):
    flip = False
    BOT = 0
    TOP = 2
    MID = 1

    points = [ (polygons[i][0], polygons[i][1], polygons[i][2], colors[0]),
               (polygons[i+1][0], polygons[i+1][1], polygons[i+1][2], colors[1]),
               (polygons[i+2][0], polygons[i+2][1], polygons[i+2][2], colors[2]) ]

    points.sort(key = lambda x: x[1])
    x0 = points[BOT][0]
    z0 = points[BOT][2]
    x1 = points[BOT][0]
    z1 = points[BOT][2]
    y = int(points[BOT][1])

    distance0 = int(points[TOP][1]) - y * 1.0 + 1
    distance1 = int(points[MID][1]) - y * 1.0 + 1
    distance2 = int(points[TOP][1]) - int(points[MID][1]) * 1.0 + 1

    dx0 = (points[TOP][0] - points[BOT][0]) / distance0 if distance0 != 0 else 0
    dz0 = (points[TOP][2] - points[BOT][2]) / distance0 if distance0 != 0 else 0
    dx1 = (points[MID][0] - points[BOT][0]) / distance1 if distance1 != 0 else 0
    dz1 = (points[MID][2] - points[BOT][2]) / distance1 if distance1 != 0 else 0


    x0_color = points[BOT][3][0]
    y0_color = points[BOT][3][1]
    z0_color = points[BOT][3][2]
    x1_color = points[BOT][3][0]
    y1_color = points[BOT][3][1]
    z1_color = points[BOT][3][2]
    dx0_color = (points[TOP][3][0] - points[BOT][3][0]) / distance0 if distance0 != 0 else 0
    dy0_color = (points[TOP][3][1] - points[BOT][3][1]) / distance0 if distance0 != 0 else 0
    dz0_color = (points[TOP][3][2] - points[BOT][3][2]) / distance0 if distance0 != 0 else 0
    dx1_color = (points[MID][3][0] - points[BOT][3][0]) / distance1 if distance1 != 0 else 0
    dy1_color = (points[MID][3][1] - points[BOT][3][1]) / distance1 if distance1 != 0 else 0
    dz1_color = (points[MID][3][2] - points[BOT][3][2]) / distance1 if distance1 != 0 else 0
    while y <= int(points[TOP][1]):
        if ( not flip and y >= int(points[MID][1])):
            flip = True
            dx1 = (points[TOP][0] - points[MID][0]) / distance2 if distance2 != 0 else 0
            dz1 = (points[TOP][2] - points[MID][2]) / distance2 if distance2 != 0 else 0
            x1 = points[MID][0]
            z1 = points[MID][2]
            dx1_color = (points[TOP][3][0] - points[MID][3][0]) / distance2 if distance2 != 0 else 0
            dy1_color = (points[TOP][3][1] - points[MID][3][1]) / distance2 if distance2 != 0 else 0
            dz1_color = (points[TOP][3][2] - points[MID][3][2]) / distance2 if distance2 != 0 else 0
            x1_color = points[MID][3][0]
            y1_color = points[MID][3][1]
            z1_color = points[MID][3][2]
        l_0 = [x0_color, y0_color, z0_color]
        l_1 = [x1_color, y1_color, z1_color]
        draw_scanline_phong(int(x0), z0, int(x1), z1, y, screen, zbuffer, view, ambient, light, symbols, reflect,l_0,l_1)
        x0+= dx0
        z0+= dz0
        x1+= dx1
        z1+= dz1
        y+= 1
        x0_color += dx0_color
        y0_color += dy0_color
        z0_color += dz0_color
        x1_color += dx1_color
        y1_color += dy1_color
        z1_color += dz1_color

def draw_scanline_gouraud(x0, z0, x1, z1, y, screen, zbuffer, l_0, l_1):
    if x0 > x1:
        tx = x0
        tz = z0
        x0 = x1
        z0 = z1
        x1 = tx
        z1 = tz
        temp = l_0
        l_0 = l_1
        l_1 = temp

    x = x0
    z = z0
    dz = (z1 - z0) / (x1 - x0 + 1) if (x1 - x0 + 1) != 0 else 0

    distance = x1 - x0 + 1
    x_color = l_0[0]
    y_color = l_0[1]
    z_color = l_0[2]
    if(distance == 0):
        dx_color = 0
        dy_color = 0
        dz_color = 0
    else:
        dx_color = (l_1[0] - l_0[0]) / distance
        dy_color = (l_1[1] - l_0[1]) / distance
        dz_color = (l_1[2] - l_0[2]) / distance
    while (x <= x1):
        color = [ int(x_color), int(y_color), int(z_color)]
        plot(screen, zbuffer, color, x, y, z)
        x+= 1
        z+= dz
        x_color += dx_color
        y_color += dy_color
        z_color += dz_color

def scanline_gouraud(polygons, i, screen, zbuffer, colors):
    flip = False
    BOT = 0
    TOP = 2
    MID = 1
    points = [ (polygons[i][0], polygons[i][1], polygons[i][2], colors[0]),
               (polygons[i+1][0], polygons[i+1][1], polygons[i+1][2], colors[1]),
               (polygons[i+2][0], polygons[i+2][1], polygons[i+2][2], colors[2]) ]
    points.sort(key = lambda x: x[1])
    x0 = points[BOT][0]
    z0 = points[BOT][2]
    x1 = points[BOT][0]
    z1 = points[BOT][2]
    y = int(points[BOT][1])

    distance0 = int(points[TOP][1]) - y * 1.0 + 1
    distance1 = int(points[MID][1]) - y * 1.0 + 1
    distance2 = int(points[TOP][1]) - int(points[MID][1]) * 1.0 + 1

    dx0 = (points[TOP][0] - points[BOT][0]) / distance0 if distance0 != 0 else 0
    dz0 = (points[TOP][2] - points[BOT][2]) / distance0 if distance0 != 0 else 0
    dx1 = (points[MID][0] - points[BOT][0]) / distance1 if distance1 != 0 else 0
    dz1 = (points[MID][2] - points[BOT][2]) / distance1 if distance1 != 0 else 0
    x0_color = points[BOT][3][0]
    y0_color = points[BOT][3][1]
    z0_color = points[BOT][3][2]
    x1_color = points[BOT][3][0]
    y1_color = points[BOT][3][1]
    z1_color = points[BOT][3][2]
    dx0_color = (points[TOP][3][0] - points[BOT][3][0]) / distance0 if distance0 != 0 else 0
    dy0_color = (points[TOP][3][1] - points[BOT][3][1]) / distance0 if distance0 != 0 else 0
    dz0_color = (points[TOP][3][2] - points[BOT][3][2]) / distance0 if distance0 != 0 else 0
    dx1_color = (points[MID][3][0] - points[BOT][3][0]) / distance1 if distance1 != 0 else 0
    dy1_color = (points[MID][3][1] - points[BOT][3][1]) / distance1 if distance1 != 0 else 0
    dz1_color = (points[MID][3][2] - points[BOT][3][2]) / distance1 if distance1 != 0 else 0

    while y <= int(points[TOP][1]):
        if ( not flip and y >= int(points[MID][1])):
            flip = True
            dx1 = (points[TOP][0] - points[MID][0]) / distance2 if distance2 != 0 else 0
            dz1 = (points[TOP][2] - points[MID][2]) / distance2 if distance2 != 0 else 0
            x1 = points[MID][0]
            z1 = points[MID][2]
            dx1_color = (points[TOP][3][0] - points[MID][3][0]) / distance2 if distance2 != 0 else 0
            dy1_color = (points[TOP][3][1] - points[MID][3][1]) / distance2 if distance2 != 0 else 0
            dz1_color = (points[TOP][3][2] - points[MID][3][2]) / distance2 if distance2 != 0 else 0
            x1_color = points[MID][3][0]
            y1_color = points[MID][3][1]
            z1_color = points[MID][3][2]
        l_0 = [x0_color, y0_color, z0_color]
        l_1 = [x1_color, y1_color, z1_color]
        draw_scanline_gouraud(int(x0), z0, int(x1), z1, y, screen, zbuffer, l_0, l_1)
        x0+= dx0
        z0+= dz0
        x1+= dx1
        z1+= dz1
        y+= 1
        x0_color += dx0_color
        y0_color += dy0_color
        z0_color += dz0_color
        x1_color += dx1_color
        y1_color += dy1_color
        z1_color += dz1_color

def draw_polygons( polygons, screen, zbuffer, view, ambient, lights, symbols, reflect, shading):
    light = [ lights[0][0], lights[0][1] ]
    if len(polygons) < 2:
        print 'invalid number of points'
        return

    avm = {}
    if (shading == 'gouraud' or shading == 'phong'):
        avm = calculate_vertex_norms( polygons )

    point = 0
    while point < len(polygons) - 2:

        n = calculate_normal(polygons, point)[:]

        if n[2] > 0:
            if shading == 'flat':
                color = [0,0,0]
                for light in lights:
                    c =  get_lighting(n, view, ambient, light, symbols, reflect )
                    color[0] += c[0]
                    color[1] += c[1]
                    color[2] += c[2]
                limit_color(color)
                scanline_convert(polygons, point, screen, zbuffer, color)
            elif shading == 'gouraud':
                n_0 = avm[ str(estimate(polygons[point])) ]
                l_0 = get_lighting( n_0, view, ambient, light, symbols, reflect  )
                n_1 = avm[ str(estimate(polygons[point+1])) ]
                l_1 = get_lighting( n_1, view, ambient, light, symbols, reflect  )
                n_2 = avm[ str(estimate(polygons[point+2])) ]
                l_2 = get_lighting( n_2, view, ambient, light, symbols, reflect  )
                colors = [ l_0, l_1, l_2 ]
                scanline_gouraud(polygons, point, screen, zbuffer, colors)
            elif shading == 'phong':
                n_0 = avm[ str(estimate(polygons[point])) ]
                n_1 = avm[ str(estimate(polygons[point+1])) ]
                n_2 = avm[ str(estimate(polygons[point+2])) ]
                colors = [n_0,n_1,n_2]
                scanline_phong(polygons, point, screen, zbuffer, colors, view, ambient, light, symbols, reflect)
        point+= 3


def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z)
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z)

    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1)
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1)

    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1)
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1)
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z)
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z)

    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1)
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z)
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z)
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1)

def add_sphere(polygons, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)

    lat_start = 0
    lat_stop = step
    longt_start = 0
    longt_stop = step

    step+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * step + longt
            p1 = p0+1
            p2 = (p1+step) % (step * (step-1))
            p3 = (p0+step) % (step * (step-1))

            if longt != step - 2:
                add_polygon( polygons, points[p0][0],
                             points[p0][1],
                             points[p0][2],
                             points[p1][0],
                             points[p1][1],
                             points[p1][2],
                             points[p2][0],
                             points[p2][1],
                             points[p2][2])
            if longt != 0:
                add_polygon( polygons, points[p0][0],
                             points[p0][1],
                             points[p0][2],
                             points[p2][0],
                             points[p2][1],
                             points[p2][2],
                             points[p3][0],
                             points[p3][1],
                             points[p3][2])


def generate_sphere( cx, cy, cz, r, step ):
    points = []

    rot_start = 0
    rot_stop = step
    circ_start = 0
    circ_stop = step

    for rotation in range(rot_start, rot_stop):
        rot = rotation/float(step)
        for circle in range(circ_start, circ_stop+1):
            circ = circle/float(step)

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points

def add_torus(polygons, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)

    lat_start = 0
    lat_stop = step
    longt_start = 0
    longt_stop = step

    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * step + longt;
            if (longt == (step - 1)):
                p1 = p0 - longt;
            else:
                p1 = p0 + 1;
            p2 = (p1 + step) % (step * step);
            p3 = (p0 + step) % (step * step);

            add_polygon(polygons,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(polygons,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )


def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    rot_start = 0
    rot_stop = step
    circ_start = 0
    circ_stop = step

    for rotation in range(rot_start, rot_stop):
        rot = rotation/float(step)
        for circle in range(circ_start, circ_stop):
            circ = circle/float(step)

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points


def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    i = 1

    while i <= step:
        t = float(i)/step
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;

        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        i+= 1

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    i = 1
    while i <= step:
        t = float(i)/step
        x = t * (t * (xcoefs[0] * t + xcoefs[1]) + xcoefs[2]) + xcoefs[3]
        y = t * (t * (ycoefs[0] * t + ycoefs[1]) + ycoefs[2]) + ycoefs[3]
        #x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        #y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]

        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        i+= 1


def draw_lines( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 2 points to draw'
        return

    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   matrix[point][2],
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   matrix[point+1][2],
                   screen, zbuffer, color)
        point+= 2

def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)

def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )



def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):

    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x + 1
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y) + 1
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    dz = (z1 - z0) / distance if distance != 0 else 0

    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):

            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east
        z+= dz
        loop_start+= 1
    plot( screen, zbuffer, color, x, y, z )
