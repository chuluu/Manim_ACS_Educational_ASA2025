# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 11:54:21 2024

@author: mbluu
"""
from manim import *
CONFIG = {"include_numbers": False,
          "include_tip": False}

MASS_RADIUS = 1.0*0.5
SPRING = 1
MASS = 1
MASS_RADIUS_OSC = 0.5*0.5

def distance_formula(start,end):
    return np.sqrt( (end[0] - start[0])**2 + (end[1] - start[1])**2  + (end[2] - start[2])**2 )

def generate_mode_shape_funcs(mode,omega,damp=1.0):
    #return [(lambda t,bb=mode[ll],cc=2*(ll+1): (bb*np.exp(-1.2*t/2)*np.exp(1j*omega*t)+cc)) for ll in range(len(mode))]
    return [(lambda t,bb=mode[ll],cc=2*(ll+1): (bb*np.exp(-damp*t/2)*np.exp(1j*omega*t))) for ll in range(len(mode))]


def dof1mass(funcs,omega,colors=None,mass_radius = 2*MASS_RADIUS):
    total_masses = len(funcs)
    max_bnd = (total_masses+1)*2
    min_bnd = 0
    axes = (
        Axes(
            x_range = [min_bnd,max_bnd,2],
            x_length = 12,
            y_range = [0,20,4],
            y_length = 2.4,
            axis_config=CONFIG,
            ).set_color(GREY).shift(1.5*UP)
        )
    grid_labels_OG = axes.get_axis_labels(
        Tex("x").scale(0.7), Tex("y").scale(0.7)
    )
    
    " Set up masses section"
    x_mass = ValueTracker(0)
    masses=VGroup()
    final_func = []
    
    for ii in range(total_masses):
        func = axes.plot(funcs[ii], x_range = [0,10], color = YELLOW)
        final_func.append(func)
        

    wall1 = Line(axes.c2p(min_bnd,0),axes.c2p(min_bnd,20))
    masses = VGroup(*[always_redraw( lambda ii=ii: DrawMass(x_mass,final_func[ii],RED)) for ii in range(total_masses)])


    

    springs = VGroup()
    spring2 = always_redraw(lambda ii=ii: DrawSpring(wall1.get_center(), masses[0][0].get_start() - [mass_radius,0,0]))
    springs.add(spring2)
    for ii in range(len(masses)-1):
        if ii == total_masses-1:
            spring2 = always_redraw(lambda ii=ii: DrawSpring(masses[ii][0].get_start(), wall2.get_center()))

        else:
            spring2 = always_redraw(lambda ii=ii: DrawSpring(masses[ii][0].get_end(), masses[ii+1][0].get_end()- [mass_radius,0,0]))
            
        springs.add(spring2)

    masses.remove(masses[-1])

    masses_walls = VGroup(masses,wall1,wall2)
    
    return x_mass,masses_walls, springs

class Spring(VMobject):
    def __init__(self, start=ORIGIN, end=ORIGIN+RIGHT, bumps=8):
        self.end = end
        self.length = distance_formula(start,end)
        self.empty = self.length/(bumps+2)
        self.step = self.length/(bumps+2)/2
        self.bump = 0.18
        super().__init__(color=WHITE)
        vertices = np.array(
            [
                start,
                start + [self.empty, 0, 0],
                start + [self.empty + self.step, self.bump, 0],
                *[start + 
                    [
                        self.empty + self.step + self.step * 2 * i,
                        self.bump * (1 - (i % 2) * 2),
                        0,
                    ]
                    for i in range(1, bumps)
                ],
                start + [self.empty + self.step * 2 * bumps, 0, 0],
                end,
            ]
        )
        vertices = vertices 


        self.start_new_path(np.array(start))
        self.add_points_as_corners(
            [*(np.array(vertex) for vertex in vertices)])

class SpringVert(VMobject):
    def __init__(self, start=ORIGIN, end=ORIGIN+RIGHT, bumps=8):
        self.end = end
        self.length = distance_formula(start,end)
        self.empty = self.length/(bumps+2)
        self.step = self.length/(bumps+2)/2
        self.bump = 0.18
        super().__init__(color=WHITE)
        vertices = np.array(
            [
                start,
                start + [0, self.empty, 0],
                start + [self.bump, self.empty + self.step, 0],
                *[start + 
                    [
                        self.bump * (1 - (i % 2) * 2),
                        self.empty + self.step + self.step * 2 * i,
                        0,
                    ]
                    for i in range(1, bumps)
                ],
                start + [ 0,self.empty + self.step * 2 * bumps, 0],
                end,
            ]
        )
        vertices = vertices 


        self.start_new_path(np.array(start))
        self.add_points_as_corners(
            [*(np.array(vertex) for vertex in vertices)])

def DrawMass(x_mass,func,mass_radius,axxx,color=RED):
    pt_move = axxx.c2p(func.underlying_function(    x_mass.get_value()   ))
    pt_move2 = pt_move @ np.array([[np.cos(np.pi/2), -np.sin(np.pi/2), 0   ],
                                   [np.sin(np.pi/2), np.cos(np.pi/2), 0   ],
                                   [0, 0, 1   ]]) 
    circ = Square(side_length = mass_radius, color = color, fill_opacity=1).scale(0.5).move_to(pt_move2)
    
    text = MathTex(r"m", color = BLACK).move_to(circ.get_center())
    return VGroup(circ,text)

def DrawSpring(p1,p2):
    sprng = SpringVert(p1,p2) 
    text = MathTex(r"k", color = WHITE).move_to(sprng.get_center()).shift(LEFT/2)
    
    return VGroup(sprng,text)

def DrawDashpot(p1,p2):
    dasher_sizer = 0.5
    p1 = np.array(p1)
    p2 = np.array(p2)
    p2_l1 = p2.copy()
    p2_l1[1] = (p1[1]+p2[1])/2 - dasher_sizer/2
    
    
    p1_l2 = p2_l1.copy()
    p1_l2[1] = p1_l2[1] + dasher_sizer
   
    dasher_center = p2_l1.copy()
    dasher_center[1] = dasher_center[1]  + dasher_sizer/2
    
    dasher = Square(dasher_sizer).move_to(dasher_center)
    
    dasher_line_pt1 = p2_l1.copy()
    dasher_line_pt1[0] = dasher_line_pt1[0] - dasher_sizer/2
    dasher_line_pt3 = dasher_line_pt1.copy()
    dasher_line_pt3[1] = dasher_line_pt3[1] + dasher_sizer + dasher_sizer/2
    
    
    dasher_line_pt2 = p2_l1.copy()
    dasher_line_pt2[0] = dasher_line_pt2[0] + dasher_sizer/2
    dasher_line_pt4 = dasher_line_pt2.copy()
    dasher_line_pt4[1] = dasher_line_pt4[1] + dasher_sizer + dasher_sizer/2
    
    dasher_l1 = Line(dasher_line_pt1,dasher_line_pt3)
    dasher_l2 = Line(dasher_line_pt2,dasher_line_pt4)

    
    line1 = Line(p1,p2_l1) 
    line2 = Line(p1_l2,p2) 

    text = MathTex(r"c", color = WHITE).move_to(dasher.get_center()).shift(RIGHT/2)
    
    return VGroup(line1,line2,dasher,dasher_l1,dasher_l2,text)#@Line(p1,p2)

class mass_damped(Scene):
    def setup_math_scene(self):
        text1 = Text('Suppose we had a mass-spring-damping system')
        self.play(Write(text1))
        self.wait(2)
        self.play(text1.animate.to_edge(UP))
        
        return text1
    
    def scene_sys_setup(self, w, period, factor, start_time,damper, time):
        mass_radius = 3*MASS_RADIUS
        axes = (
            Axes(
                x_range = [-10,10,5],
                x_length = 6,
                y_range = [-10,10,5],
                y_length = 6,
                axis_config=CONFIG,
                ).set_color(GREY)
            )
   
        func_1dof    = generate_mode_shape_funcs([-5],w,damp=damper.get_value())
        func_sine    = generate_mode_shape_funcs([5],w,damp=damper.get_value())  
        mass_motion = axes.plot(func_1dof[0], x_range = [-10,10], color = YELLOW)

        axes_timedomain = always_redraw(lambda: Axes(
                x_range = [start_time,start_time*2+period*factor,period],
                x_length = 5,
                y_range = [-10,10,5],
                y_length = 6,
                axis_config=CONFIG,
                ).set_color(GREY).next_to(axes,RIGHT).scale(1).shift(2*LEFT))
        
        equibline = DashedLine(axes.c2p(-5,0),axes.c2p(5,0), dash_length=0.2, color=BLUE)
        sineplot = always_redraw(lambda: axes_timedomain.plot(func_sine[0], x_range = [start_time,time.get_value()], color = YELLOW))
        wall1 = Line(axes.c2p(-5,-10),axes.c2p(5,-10))
        mass = always_redraw( lambda: DrawMass(time,mass_motion,mass_radius,axes))
        spring  = always_redraw(lambda: DrawSpring(wall1.get_center() - [mass_radius/4,0,0], mass[0].get_start() - [mass_radius/2,mass_radius/2,0]))
        dashpot = always_redraw(lambda: DrawDashpot(wall1.get_center() + [mass_radius/4,0,0], mass[0].get_start() - [0,mass_radius/2,0]))
        asd  = VGroup(wall1,mass,spring,dashpot)
        
        return axes, axes_timedomain, sineplot, asd, equibline
    
    def construct(self):
        f = 1
        w = 2*np.pi*f
        period = 1/f
        factor = 6
        start_time = period/4
        damper = ValueTracker(0)
        time = ValueTracker(start_time)
        [axes, axes_timedomain, sineplot, asd, equibline] = self.scene_sys_setup(w, period, factor, start_time,damper,time)
        
        damper2 = ValueTracker(0.7)
        time2 = ValueTracker(start_time)
        [axes2, axes_timedomain2, sineplot2, asd2, equibline2] = self.scene_sys_setup(w, period, factor, start_time,damper2,time2)

        damper3 = ValueTracker(2.0)
        time3 = ValueTracker(start_time)
        [axes3, axes_timedomain3, sineplot3, asd3, equibline3] = self.scene_sys_setup(w, period, factor, start_time,damper3,time3)


        damper4 = ValueTracker(-0.2)
        time4 = ValueTracker(start_time)
        [axes4, axes_timedomain4, sineplot4, asd4, equibline4] = self.scene_sys_setup(w, period, factor, start_time,damper4,time4)

        
        force_vector = Arrow(start=UP, end=DOWN, color=RED).next_to(asd[1],UP)
        force_external = MathTex(r'F_{ext}', color=RED).next_to(force_vector,RIGHT)
        force_vec = VGroup(force_vector,force_external)



        stiff_var = Variable(1, MathTex(r"k"), num_decimal_places=1).next_to(axes,LEFT)
        mass_var = Variable(1, MathTex(r"m"), num_decimal_places=1).next_to(stiff_var,DOWN)   
        
        damping_var = Variable(damper.get_value(), MathTex(r"c"), num_decimal_places=1).next_to(mass_var,DOWN)
        var_tracker = damping_var.tracker

        
        text1 = self.setup_math_scene()
        asd_math = asd.copy().shift(LEFT*4)
        force_vec_2 = VGroup(force_vector.copy(),force_external.copy()).next_to(asd_math,UP).shift(RIGHT/2)
        
        write_sum = Text('Sum Force').next_to(text1,DOWN)
        math_sum = MathTex(r'\sum F = ').next_to(asd_math)
        
        write_mass = Text('Newtons Second Law').next_to(text1,DOWN)
        mass_term = MathTex(r'm \ddot{x}').next_to(math_sum,RIGHT)
        
        
        write_spring = Text('Hookes Law').next_to(text1,DOWN)
        spring_term = MathTex(r'-kx').next_to(mass_term,RIGHT)
        
        write_damp = Text('Viscous Damping').next_to(text1,DOWN)
        damp_term = MathTex(r'-c \dot{x}').next_to(spring_term,RIGHT)
        
        write_F = Text('External Force').next_to(text1,DOWN)
        F_term = MathTex(r'+F_{ext}').next_to(damp_term,RIGHT)
        
        big_var = VGroup(stiff_var,mass_var,damping_var)
        
        math_eqn = VGroup(math_sum,mass_term,spring_term,damp_term,F_term)
        
        self.play(FadeIn(asd_math))
        self.wait(2)
        self.play(Write(math_sum),Write(write_sum))
        self.wait(1)
        self.play(FadeOut(write_sum))
        self.wait(1)
        self.play(Indicate(asd_math[1]),Write(write_mass))
        self.play(Write(mass_term))
        self.wait(1)
        self.play(FadeOut(write_mass))
        self.wait(1)
        self.play(Indicate(asd_math[2]),Write(write_spring))
        self.play(Write(spring_term))
        self.wait(1)
        self.play(FadeOut(write_spring))
        self.wait(1)
        self.play(Indicate(asd_math[3]),Write(write_damp))
        self.play(Write(damp_term))
        self.wait(1)
        self.play(FadeOut(write_damp))
        self.wait(1)
        self.play(FadeIn(force_vec_2),Write(write_F))
        self.wait(1)
        self.play(Write(F_term)) 
        self.wait(2)
        
        #undamped
        text_undamped = Text('Undamped system').to_edge(UP)
        self.play(ReplacementTransform(text1,text_undamped))
        self.play(FadeOut(write_F),FadeOut(force_vec_2),FadeIn(equibline),FadeIn(axes_timedomain),
                  FadeIn(sineplot),FadeIn(big_var),FadeOut(math_eqn),ReplacementTransform(asd_math,asd))
        self.wait(1)
        self.play(Write(force_vec))
        self.wait(1)
        self.play(time.animate(rate_func=rate_functions.linear,run_time=1).increment_value(period/4),
                  FadeOut(force_vec))
        self.wait(2)
        self.play(time.animate(rate_func=rate_functions.linear,run_time=4).increment_value(period*factor))
        
        #
        text_underdamped = Text('Underdamped system c < 2.0').to_edge(UP)
        self.wait(2)
        self.play(ReplacementTransform(text_undamped,text_underdamped))
        self.play(var_tracker.animate.set_value(0.7))
        self.play(FadeOut(sineplot),FadeIn(sineplot2),ReplacementTransform(axes_timedomain,axes_timedomain2),
                  ReplacementTransform(asd,asd2))
        self.wait(1)
        self.play(Write(force_vec))
        self.wait(1)
        self.play(time2.animate(rate_func=rate_functions.linear,run_time=1).increment_value(period/4),
                  FadeOut(force_vec))
        self.wait(2)
        self.play(time2.animate(rate_func=rate_functions.linear,run_time=4).increment_value(period*factor))
        
        #
        text_critical = Text('Critically damped system c = 2.0').to_edge(UP)
        self.wait(2)
        self.play(ReplacementTransform(text_underdamped,text_critical))
        self.play(var_tracker.animate.set_value(2.0))
        self.play(FadeOut(sineplot2),FadeIn(sineplot3),ReplacementTransform(axes_timedomain2,axes_timedomain3),
                  ReplacementTransform(asd2,asd3))
        self.wait(1)
        self.play(Write(force_vec))
        self.wait(1)
        self.play(time3.animate(rate_func=rate_functions.linear,run_time=1).increment_value(period/4),
                  FadeOut(force_vec))
        self.wait(2)
        self.play(time3.animate(rate_func=rate_functions.linear,run_time=4).increment_value(period*factor))
        self.wait(1)
        
        text_unstable = Text('System unstable negative damping').to_edge(UP)
        self.wait(2)
        self.play(ReplacementTransform(text_critical,text_unstable))
        self.play(var_tracker.animate.set_value(-0.2))
        self.play(FadeOut(sineplot3),FadeIn(sineplot4),ReplacementTransform(axes_timedomain3,axes_timedomain4),
                  ReplacementTransform(asd3,asd4))
        self.wait(1)
        self.play(Write(force_vec))
        self.wait(1)
        self.play(time4.animate(rate_func=rate_functions.linear,run_time=1).increment_value(period/4),
                  FadeOut(force_vec))
        self.wait(2)
        self.play(time4.animate(rate_func=rate_functions.linear,run_time=4).increment_value(period*factor))
        self.wait(1)
        
        
        
        