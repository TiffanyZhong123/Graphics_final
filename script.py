import mdl
import sys
from display import *
from matrix import *
from draw import *

def first_pass( commands ):

    name = 'base'
    num_frames = 1

    frames_exist = 0
    vary_exist = 0
    base_exist = 0

    for command in commands:
        # print command
        c = command['op']
        args = command['args']
        if c == 'frames':
            num_frames = int(args[0])
            frames_exist = 1
        elif c == 'basename':
            name = args[0]
            base_exist = 1
        elif c == 'vary':
            vary_exist = 1

    if frames_exist == 0 and vary_exist == 1:
        print 'frame not found, exiting'
        sys.exit()
    if frames_exist == 1 and base_exist == 0:
        print 'using "base"'

    return (name, num_frames)

def second_pass( commands, num_frames ):
    frames = [ {} for i in range(num_frames) ]
    for command in commands:
        if command['op'] == 'vary':
            start_frame = int( command['args'][0] )
            end_frame = int( command['args'][1] )
            start_value = command['args'][2]
            end_value = command['args'][3]
            delta = (end_value - start_value) / ( end_frame - start_frame )
            for i in range( end_frame - start_frame + 1):
                frames[start_frame + i][command['knob']] = start_value + delta * i
    return frames


def run(filename):
    p = mdl.parseFile(filename)
    if p:
        (commands, symbols) = p
    else:
        print "Parsing failed."
        return

    view = [0,
            0,
            1];
    ambient = [50,
               50,
               50]
    lights = []
    color = [0, 0, 0]
    symbols['.white'] = ['constants',
                         {'red': [0.2, 0.5, 0.5],
                        'green': [0.2, 0.5, 0.5],
                          'blue': [0.2, 0.5, 0.5]}]
    reflect = '.white'

    [name, num_frames] = first_pass(commands)
    knobs = second_pass(commands, num_frames)

    '''Set shading type, default is flat'''
    shade_type = 'flat'
    if 'shading' in symbols:
        shade_type = symbols['shading'][1]
    for i in range(int(num_frames)):
        tmp = new_matrix()
        ident( tmp )

        stack = [ [x[:] for x in tmp] ]
        screen = new_screen()
        zbuffer = new_zbuffer()
        tmp = []
        step_3d = 20
        consts = ''
        coords = []
        coords1 = []

        for knob in knobs[i]:
            symbols[knob][1] = knobs[i][knob]

        for command in commands:
            c = command['op']
            args = command['args']
            knob_value = 1

            if c == 'light':
                s = symbols[command['light']]
                sample_color = s[1]['color'][:]
                sample_location = s[1]['location'][:]
                if command['knob']:
                    knob1_value = symbols[command["knob"][0]][1]
                    knob2_value = symbols[command["knob"][1]][1]
                    if command["knob"][0] == "k0":
                        s[1]['color'][0] = min(s[1]['color'][0] * knob1_value, 255)
                        s[1]['color'][1] = min(s[1]['color'][1] * knob1_value, 255)
                        s[1]['color'][2] = min(s[1]['color'][2] * knob1_value, 255)
                    if command["knob"][1] == "k1":
                        s[1]['location'][0] = min(s[1]['location'][0] * knob2_value, 1)
                to_remove = -1
                for j in range(len(lights)):
                    sym = lights[0][2]
                    if sym == command['light']:
                        to_remove = j
                if to_remove >= 0:
                    lights.pop(j)
                lights.append([s[1]['location'], s[1]['color'], command['light']])
                s[1]['color'] = sample_color
                s[1]['location'] = sample_location
            if c == 'mesh':
                if command['constants'] and command['constants'] != ":":
                    reflect = command['constants']
                add_mesh(tmp, args[0])
                matrix_mult( stack[-1], tmp )
                draw_polygons(tmp, screen, zbuffer, view, ambient, lights, symbols, reflect, shade_type)
                tmp = []
                reflect = '.white'
            if c == 'box':
                if command['constants']:
                    reflect = command['constants']
                add_box(tmp,
                        args[0], args[1], args[2],
                        args[3], args[4], args[5])
                matrix_mult( stack[-1], tmp )
                draw_polygons(tmp, screen, zbuffer, view, ambient, lights, symbols, reflect, shade_type)
                tmp = []
                reflect = '.white'
            elif c == 'sphere':
                if command['constants']:
                    reflect = command['constants']
                add_sphere(tmp,
                           args[0], args[1], args[2], args[3], step_3d)
                matrix_mult( stack[-1], tmp )
                draw_polygons(tmp, screen, zbuffer, view, ambient, lights, symbols, reflect, shade_type)
                tmp = []
                reflect = '.white'
            elif c == 'torus':
                if command['constants']:
                    reflect = command['constants']
                add_torus(tmp,
                          args[0], args[1], args[2], args[3], args[4], step_3d)
                matrix_mult( stack[-1], tmp )
                draw_polygons(tmp, screen, zbuffer, view, ambient, lights, symbols, reflect, shade_type)
                tmp = []
                reflect = '.white'
            elif c == 'line':
                add_edge(tmp,
                         args[0], args[1], args[2], args[3], args[4], args[5])
                matrix_mult( stack[-1], tmp )
                draw_lines(tmp, screen, zbuffer, color)
                tmp = []
            elif c == 'move':
                if command["knob"]:
                    knob_value = symbols[command["knob"]][1]
                tmp = make_translate(args[0] * knob_value, args[1] * knob_value, args[2] * knob_value)
                matrix_mult(stack[-1], tmp)
                stack[-1] = [x[:] for x in tmp]
                tmp = []
            elif c == 'scale':
                if command["knob"]:
                    knob_value = symbols[command["knob"]][1]
                tmp = make_scale(args[0] * knob_value, args[1]* knob_value, args[2]* knob_value)
                matrix_mult(stack[-1], tmp)
                stack[-1] = [x[:] for x in tmp]
                tmp = []
            elif c == 'rotate':
                if command["knob"]:
                    knob_value = symbols[command["knob"]][1]
                theta = args[1] * (math.pi/180) * knob_value
                if args[0] == 'x':
                    tmp = make_rotX(theta)
                elif args[0] == 'y':
                    tmp = make_rotY(theta)
                else:
                    tmp = make_rotZ(theta)
                matrix_mult( stack[-1], tmp )
                stack[-1] = [ x[:] for x in tmp]
                tmp = []
            elif c == 'push':
                stack.append([x[:] for x in stack[-1]] )
            elif c == 'pop':
                stack.pop()
            elif c == 'display':
                display(screen)
            elif c == 'save':
                save_extension(screen, args[0])
        save_extension(screen,'anim/' + name + ('%03d' %int(i)))

    if num_frames > 1:
        make_animation(name)
