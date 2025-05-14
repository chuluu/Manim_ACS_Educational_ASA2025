# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 13:10:34 2023

@author: mbluu
"""
from manim import *
import numpy as np

CONFIG = {"include_numbers": False,
          "include_tip": False}

MASS_RADIUS = 1.0*0.5
SPRING = 1
MASS = 1
MASS_RADIUS_OSC = 0.5*0.5


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


def emphasis_box(term,color=YELLOW):
    box = SurroundingRectangle(term, corner_radius=0.2,color=color)
    emphasized_term = VGroup(term,box)
    
    return emphasized_term

def seteqnmotioncolors(eqnmotion_tot,st_arr,colores_arr,colors,initial):
    for jj in range(len(eqnmotion_tot)):
        eqnmotion_tot[jj][0][initial[0]:initial[1]].set_color(colors[jj])
    
    for jj in range(len(eqnmotion_tot)):
        for st,color in zip(st_arr[jj],colores_arr[jj]):
            eqnmotion_tot[jj][1][st:st+2].set_color(color)



def equation_generation(colors,eqn_seq = 1):
    if eqn_seq == 1:
        eqnmotion1 = MathTex(r"m\ddot{x}_1", r" = -k x_1 - k \left(x_1 - x_2\right)").shift(DOWN)
        eqnmotion2 = MathTex(r"m\ddot{x}_2", r" = -k \left(x_2 - x_1\right) - k \left(x_2 - x_3\right)").next_to(eqnmotion1,DOWN)
        eqnmotion3 = MathTex(r"m\ddot{x}_3", r" = -k \left(x_3 - x_2\right) - k x_3").next_to(eqnmotion2,DOWN)
        initial = [1,5]
        st_arr = [[3,8,11],[4,7,13,16],[4,7,12]]


    elif eqn_seq == 2:
        eqnmotion1 = MathTex(r"-m\omega^2 a_1 e^{j\omega t}", r" = (-k a_1 - k \left(a_1 - a_2\right)) e^{j\omega t}").shift(DOWN)
        eqnmotion2 = MathTex(r"-m\omega^2 a_2 e^{j\omega t}", r" = \left(-k \left(a_2 - a_1\right) - k \left(a_2 - a_3\right) \right)e^{j\omega t}").next_to(eqnmotion1,DOWN)
        eqnmotion3 = MathTex(r"-m\omega^2 a_3 e^{j\omega t}", r" = \left(-k \left(a_3 - a_2\right) - k a_3\right)e^{j\omega t}").next_to(eqnmotion2,DOWN)
        initial = [4,6]
        st_arr = [[4,9,12],[5,8,14,17],[5,8,13]]

    
    
    elif eqn_seq == 3:
        eqnmotion1 = MathTex(r"-m\omega^2 a_1}", r" = (-k a_1 - k \left(a_1 - a_2\right))").shift(DOWN)
        eqnmotion2 = MathTex(r"-m\omega^2 a_2}", r" = \left(-k \left(a_2 - a_1\right) - k \left(a_2 - a_3\right) \right)").next_to(eqnmotion1,DOWN)
        eqnmotion3 = MathTex(r"-m\omega^2 a_3}", r" = \left(-k \left(a_3 - a_2\right) - k a_3\right)").next_to(eqnmotion2,DOWN)
        initial = [4,6]
        st_arr = [[4,9,12],[5,8,14,17],[5,8,13]]
  
    elif eqn_seq == 4:
        eqnmotion1 = MathTex(r"-\omega^2 a_1}", r" = -\omega_o^2 a_1 - \omega_o^2\left(a_1 - a_2\right)").shift(DOWN)
        eqnmotion2 = MathTex(r"-\omega^2 a_2}", r" = \omega_o^2\left(a_2 - a_1\right) - \omega_o^2\left(a_2 - a_3\right)").next_to(eqnmotion1,DOWN)
        eqnmotion3 = MathTex(r"-\omega^2 a_3}", r" = \omega_o^2\left(a_3 - a_2\right) - \omega_o^2 a_3").next_to(eqnmotion2,DOWN)

        initial = [3,5]
        st_arr = [[5,12,15],[5,8,16,19],[5,8,15]]



    colores_arr = [[colors[0],colors[0],colors[1]],  
          [colors[1],colors[0],colors[1],colors[2]],
           [colors[2],colors[1],colors[2]]]
    
    eqnmotion2[0].align_to(eqnmotion1[0],LEFT) # LEft aligning trick
    eqnmotion2[1].next_to(eqnmotion2[0],RIGHT)
    eqnmotion3[0].align_to(eqnmotion2[0],LEFT) # LEft aligning trick
    eqnmotion3[1].next_to(eqnmotion3[0],RIGHT)
    
    eqnmotion_tot = VGroup(eqnmotion1,eqnmotion2,eqnmotion3)
    

    seteqnmotioncolors(eqnmotion_tot,st_arr,colores_arr,colors,initial)

        
    return eqnmotion_tot

def finaleqnmotion(colors):
    eqnmotion1 = MathTex(r"(-\omega^2 + 2 \omega_o^2) a_1 - \omega_o^2 a_2 + 0 a_3",r" = 0").shift(DOWN)
    eqnmotion2 = MathTex(r" - \omega_o^2 a_1 + (-\omega^2 + 2 \omega_o^2) a_2 - \omega_o^2 a_3 ", r" = 0").next_to(eqnmotion1,DOWN)
    eqnmotion3 = MathTex(r"0 a_1 + -\omega_o^2 a_2 - (-\omega^2 + 2 \omega_o^2) a_3",r"= 0 ").next_to(eqnmotion2,DOWN)
    
    eqnmotion2[1].align_to(eqnmotion1[1],LEFT) # LEft aligning trick
    eqnmotion2[0].next_to(eqnmotion2[1],LEFT)
    eqnmotion3[1].align_to(eqnmotion2[1],LEFT) # LEft aligning trick
    eqnmotion3[0].next_to(eqnmotion3[1],LEFT)
    eqnmotion_tot = VGroup(eqnmotion1,eqnmotion2,eqnmotion3)
    
    
    st_arr = [[10,16,20],[4,17,23],[1,8,21]]
    colores_arr = [[colors[0],colors[1],colors[2]],  
                   [colors[0],colors[1],colors[2]],
                   [colors[0],colors[1],colors[2]]]
    
    for jj in range(len(eqnmotion_tot)):
        for st,color in zip(st_arr[jj],colores_arr[jj]):
            eqnmotion_tot[jj][0][st:st+2].set_color(color)
    

    return eqnmotion_tot

def cancel_line(text,color):
    pre_coord_dl = text.get_corner(DL)
    pre_coord_ur = text.get_corner(UR)
    reference_line = Line(pre_coord_dl,pre_coord_ur,color=color)
    return reference_line

def distance_formula(start,end):
    return np.sqrt( (end[0] - start[0])**2 + (end[1] - start[1])**2  + (end[2] - start[2])**2 )


def generate_mode_shape_funcs(mode,omega):
    return [(lambda t,bb=mode[ll],cc=2*(ll+1): (bb*np.exp(1j*omega*t)+cc)) for ll in range(len(mode))]

def generate_mode_shape_funcs_osc_graph(mode,omega):
    return [lambda t,bb=mode[ll]: bb*np.exp(1j*omega*t) for ll in range(len(mode))]

def Infinite_String(funcs,omega):
    " Setup Axes for the masses "
    total_masses = len(funcs)
    max_bnd = total_masses + 1
    min_bnd = 0
    axes = (
        Axes(
            x_range = [0,max_bnd,1],
            x_length = 12,
            y_range = [-1,1,2],
            y_length = 2.5,
            axis_config=CONFIG,
            ).set_color(GREY)
        ).shift(DOWN)
    grid_labels_OG = axes.get_axis_labels(
        Tex("x").scale(0.7), Tex("y").scale(0.7)
    )

    " Set up masses section"
    x_mass = ValueTracker(0)
    masses=VGroup()
    final_func = []
    def DrawMass(idx,x_mass,func,color=WHITE):
        circ = Dot(radius=0.01, color = color, fill_opacity=1).scale(0.5).move_to(
                 axes.c2p(idx, func.underlying_function(x_mass.get_value())))        
        return circ
    
    for ii in range(total_masses):
        func = axes.plot(funcs[ii], x_range = [0,10], color = YELLOW)
        final_func.append(func)
        
        
    wall1 = Line(axes.c2p(min_bnd,-1),axes.c2p(min_bnd,1))
    wall2 = Line(axes.c2p(max_bnd,-1),axes.c2p(max_bnd,1))
    masses = VGroup(*[always_redraw( lambda ii=ii: DrawMass(ii+1,x_mass,final_func[ii],WHITE)) for ii in range(total_masses)])
    
    masses.add(wall2)

    springs = VGroup()
    spring2 = always_redraw(lambda ii=ii: Line(wall1.get_center(), masses[0][0].get_center()))
    springs.add(spring2)
    for ii in range(len(masses)-1):
        if ii == total_masses-1:
            spring2 = always_redraw(lambda ii=ii: Line(masses[ii][0].get_center(), wall2.get_center()))

        else:
            spring2 = always_redraw(lambda ii=ii: Line(masses[ii][0].get_center(), masses[ii+1][0].get_center()))
            
        springs.add(spring2)
        
    masses.remove(masses[-1])

    masses.add(wall1)
    masses.add(wall2)

    masses_new  = VGroup(wall1,wall2)

    return axes,x_mass,masses,springs

def setup_3DOF_Oscillation_graph(funcs,omega,colors=None,mass_radius=2*MASS_RADIUS_OSC):
    " Setup Axes for the masses "
    total_masses = len(funcs)
    max_bnd = total_masses + 1
    min_bnd = 0
    axes = (
        Axes(
            x_range = [0,max_bnd,1],
            x_length = 12,
            y_range = [-1,1,2],
            y_length = 2.5,
            axis_config=CONFIG,
            ).set_color(GREY)
        ).shift(DOWN)
    grid_labels_OG = axes.get_axis_labels(
        Tex("x").scale(0.7), Tex("y").scale(0.7)
    )

    " Set up masses section"
    x_mass = ValueTracker(0)
    masses=VGroup()
    final_func = []
    def DrawMass(idx,x_mass,func,color=RED):
        circ = Circle(radius = mass_radius, color = color, fill_opacity=1).scale(0.5).move_to(
                 axes.c2p(idx, func.underlying_function(x_mass.get_value())))
        text = MathTex(r"m", color = BLACK).scale(0.5).move_to(circ.get_center())
        
        return VGroup(circ,text)
    
    for ii in range(total_masses):
        func = axes.plot(funcs[ii], x_range = [0,10], color = YELLOW)
        final_func.append(func)
        
        
    wall1 = Line(axes.c2p(min_bnd,-1),axes.c2p(min_bnd,1))
    wall2 = Line(axes.c2p(max_bnd,-1),axes.c2p(max_bnd,1))
    if colors == None:
        masses = VGroup(*[always_redraw( lambda ii=ii: DrawMass(ii+1,x_mass,final_func[ii],RED)) for ii in range(total_masses)])
    else:
        masses = VGroup(*[always_redraw( lambda ii=ii: DrawMass(ii+1,x_mass,final_func[ii],colors[ii])) for ii in range(total_masses)])

    
    masses.add(wall2)

    springs = VGroup()
    spring2 = always_redraw(lambda ii=ii: Line(wall1.get_center(), masses[0][0].get_start() - [mass_radius,0,0]))
    springs.add(spring2)
    for ii in range(len(masses)-1):
        if ii == total_masses-1:
            spring2 = always_redraw(lambda ii=ii: Line(masses[ii][0].get_start(), wall2.get_center()))

        else:
            spring2 = always_redraw(lambda ii=ii: Line(masses[ii][0].get_end(), masses[ii+1][0].get_end()- [mass_radius,0,0]))
            
        springs.add(spring2)
        
    masses.remove(masses[-1])

    masses.add(wall1)
    masses.add(wall2)

    return axes,x_mass,masses,springs

def setup_3DOF_spring_mass_system(funcs,omega,colors=None,mass_radius = 2*MASS_RADIUS):
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
    def DrawMass(x_mass,func,color=RED):
        circ = Circle(radius = mass_radius, color = color, fill_opacity=1).scale(0.5).move_to(
                 axes.c2p(func.underlying_function(x_mass.get_value()),10))
        text = MathTex(r"m", color = BLACK).move_to(circ.get_center())
        
        return VGroup(circ,text)
    
    for ii in range(total_masses):
        func = axes.plot(funcs[ii], x_range = [0,10], color = YELLOW)
        final_func.append(func)
        

    wall1 = Line(axes.c2p(min_bnd,0),axes.c2p(min_bnd,20))
    wall2 = Line(axes.c2p(max_bnd,0),axes.c2p(max_bnd,20))
    
    if colors == None:
        masses = VGroup(*[always_redraw( lambda ii=ii: DrawMass(x_mass,final_func[ii],RED)) for ii in range(total_masses)])
    else:
        masses = VGroup(*[always_redraw( lambda ii=ii: DrawMass(x_mass,final_func[ii],colors[ii])) for ii in range(total_masses)])
    masses.add(wall2)


    def DrawSpring(p1,p2):
        sprng = Spring(p1,p2)
        text = MathTex(r"k", color = WHITE).move_to(sprng.get_center()).shift(UP/2)
        
        return VGroup(sprng,text)
    

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

class Speed3dof(Scene):
    def MathScene1(self,masses_walls, springs, player=True):
        colors = self.colors
        eqnmotion = equation_generation(colors,1)
        
        x_arrows = VGroup()
        for ii in range(3):
            x_arrow_p1 = Line(masses_walls[0][ii].get_center(),masses_walls[0][ii].get_center() + [0,MASS_RADIUS/2,0]).shift(MASS_RADIUS*UP)
            x_arrow_p2 = Arrow(buff=0.5,start = x_arrow_p1.get_end(),
                                end = x_arrow_p1.get_end() + [MASS_RADIUS/2,0,0],max_stroke_width_to_length_ratio=15)
            x_arrow_txt = MathTex(r"x_{}".format(ii+1)).scale(0.8).next_to(x_arrow_p2,RIGHT/2)
            x_arrows.add( VGroup(x_arrow_p1,x_arrow_p2,x_arrow_txt).set_color(colors[ii]))
        
        if player == True:
            self.play(*[FadeIn(x_arrows[ii]) for ii in range(3)])
            for ii in range(3):
                self.play(Indicate(masses_walls[0][ii]),
                          Indicate(x_arrows[ii]))
                self.play(Write(eqnmotion[ii][0]))
                self.play(Indicate(masses_walls[0][ii]),
                          Indicate(springs[ii]),
                          Indicate(springs[ii+1]))
                self.play(Write(eqnmotion[ii][1]))

            return eqnmotion,x_arrows
        
        else:
            self.add(eqnmotion,x_arrows)
            return eqnmotion,x_arrows

    def MathScene2(self,eqnmotion,x_arrows,masses_walls,player = True):
        colors = self.colors
        eqnmotion.generate_target()
        eqnmotion.target.shift(LEFT*2.5)
        
        "Equation of motion set 2"
        eqnmotion_g2 = equation_generation(colors,2).shift(LEFT*2.5)
    
        "Equation of motion set 3"
        eqnmotion_g3 = equation_generation(colors,3).shift(LEFT*2.5)
        
        "Equation of motion set 4"
        eqnmotion_g4 = equation_generation(colors,4).shift(LEFT*2.5)
        
        
        "Harmonic Equations"
        harm1 = MathTex(r"x_1",r"=", r"a_1"," e^{j\omega t").shift(DOWN).shift(RIGHT*4)
        harm1[0].set_color(self.colors[0])
        harm1[2].set_color(self.colors[0])

        harm2 = MathTex(r"x_2",r"=", r"a_2", "e^{j\omega t").next_to(harm1,DOWN)
        harm2[0].set_color(self.colors[1])
        harm2[2].set_color(self.colors[1])

        harm3 = MathTex(r"x_3",r"=", r"a_3", " e^{j\omega t").next_to(harm2,DOWN)
        harm3[0].set_color(self.colors[2])
        harm3[2].set_color(self.colors[2])

        harmonic = VGroup(harm1,harm2,harm3)

        cancel_t1 = cancel_line(eqnmotion_g2[0][0][6:10],ORANGE)
        cancel_t2 = cancel_line(eqnmotion_g2[0][1][16:20],ORANGE)
        cancel_eqn1 = VGroup(cancel_t1,cancel_t2)
    
        cancel_t1 = cancel_line(eqnmotion_g2[1][0][6:10],ORANGE)
        cancel_t2 = cancel_line(eqnmotion_g2[1][1][21:25],ORANGE)
        cancel_eqn2 = VGroup(cancel_t1,cancel_t2)
        
        cancel_t1 = cancel_line(eqnmotion_g2[2][0][6:10],ORANGE)
        cancel_t2 = cancel_line(eqnmotion_g2[2][1][16:20],ORANGE)
        cancel_eqn3 = VGroup(cancel_t1,cancel_t2)
        
        cancel_lines = VGroup(cancel_eqn1,cancel_eqn2,cancel_eqn3)

        
        phrase = Text("Let:").scale(0.7).shift(DOWN).shift(RIGHT*4.5)
        state = MathTex(r"w_o",r"=\sqrt{\frac{k}{m}}").next_to(phrase,DOWN)
        wo_condition = VGroup(phrase,state)

            
        eqn_last = finaleqnmotion(colors).shift(LEFT*2.5)

        A_mat = Matrix([["-\omega^2 + 2 \omega_o^2", "- \omega_o^2            ", "0"],
                              ["-\omega_o^2             ", "-\omega^2 + 2 \omega_o^2", " - \omega_o^2"],
                              ["0"                       , "-\omega_o^2             ", "-\omega^2 + 2 \omega_o^2"]],
                             element_alignment_corner=ORIGIN,
                             v_buff = 1.2,
                             h_buff = 2.5)
        b_mat = Matrix([["a_1"],["a_2"],["a_3"]], 
                       element_alignment_corner=ORIGIN,v_buff = 1.2,).next_to(A_mat,RIGHT).set_row_colors(colors[0],colors[1],colors[2])
                        
        equals = MathTex(r"=").next_to(b_mat,RIGHT)
        f_mat = Matrix([["0"],
                      ["0"],
                      ["0"]], 
               element_alignment_corner=ORIGIN,v_buff = 1.2,).next_to(equals,RIGHT)
                
        
                        
        systemofeqn = VGroup(A_mat,b_mat,equals,f_mat).move_to(ORIGIN).shift(2*DOWN)
        
        if player == True:
            self.play(MoveToTarget(eqnmotion),runtime=2)
            self.wait(2)
            self.play(Write(harmonic))
            self.wait(1)
            self.play(Circumscribe(harmonic))
            self.wait(1)
            self.play(Circumscribe(eqnmotion))
            self.wait(1)
            self.play(FadeOut(harmonic,target_position=eqnmotion),
                      ReplacementTransform(eqnmotion,eqnmotion_g2))
            self.wait(1)
            self.play(Write(cancel_lines))
            self.wait(1)
            self.play(FadeOut(cancel_lines),
                  TransformMatchingTex(eqnmotion_g2,eqnmotion_g3))
            self.wait(1)
            self.play(Write(wo_condition))
            self.wait(2)
            self.play(TransformMatchingTex(eqnmotion_g3,eqnmotion_g4))
            self.wait(2)
            self.play(TransformMatchingTex(eqnmotion_g4,eqn_last))
            self.wait(3)
            self.play(FadeOut(wo_condition),ReplacementTransform(eqn_last,systemofeqn))

            return systemofeqn
        else:
            self.play(ReplacementTransform(eqnmotion,systemofeqn))
            return systemofeqn
        
    def MathScene3(self,systemofeqn,x_arrows,masses_walls,springs,player = True):
        colors = self.colors
        systemofeqn.generate_target()
        systemofeqn.target.move_to(ORIGIN)
        
        K_mat = Matrix([["2 \omega_o^2", "- \omega_o^2            ", "0"],
                      ["-\omega_o^2             ", " 2 \omega_o^2", " - \omega_o^2"],
                      ["0"                       , "-\omega_o^2             ", "2 \omega_o^2"]],
                     element_alignment_corner=ORIGIN,
                     v_buff = 1.2,
                     h_buff = 1.4)
        
        subsign = MathTex("-").next_to(K_mat,RIGHT)
        
        M_mat = Matrix([["\omega^2", "0", "0"],
              ["0", "\omega^2", "0"],
              ["0" , "0", "\omega^2"]],
             element_alignment_corner=ORIGIN,
             v_buff = 1.2,
             h_buff = 1.0).next_to(subsign,RIGHT)
        

        sec1 = VGroup(K_mat, subsign, M_mat)
        
        b_mat = Matrix([["a_1"],["a_2"],["a_3"]], 
           element_alignment_corner=ORIGIN,v_buff = 1.2,).set_row_colors(colors[0],colors[1],colors[2])
                
        
        equals = MathTex(r"=").next_to(b_mat,RIGHT)
        f_mat = Matrix([["0"],
                      ["0"],
                      ["0"]], 
               element_alignment_corner=ORIGIN,v_buff = 1.2,).next_to(equals,RIGHT)
                
        sec2 = VGroup(b_mat, equals, f_mat).next_to(sec1,RIGHT)
        
        total = VGroup(sec1,sec2).move_to(ORIGIN)
        
        
        eigval = MathTex(r"\omega^2").next_to(sec1[1],RIGHT)
        M_mat = Matrix([["1", "0", "0"],
              ["0", "1", "0"],
              ["0" , "0", "1"]],
             element_alignment_corner=ORIGIN,
             v_buff = 1.2,
             h_buff = 1.0).next_to(eigval,RIGHT)
        Newmat = VGroup(eigval,M_mat)
        
        sec1 = VGroup(K_mat, subsign, Newmat)
        sec2 = total[1].next_to(Newmat,RIGHT)
        total2 = VGroup(sec1,sec2)

        Eigprobtext = Text("Eigen-Problem").to_edge(UP)
        divider = Line([-4,0,0], [4,0,0]).next_to(Eigprobtext,DOWN)
        title = VGroup(Eigprobtext,divider)
        
        eigprob = MathTex(r"(",r"\textbf{A}","-",r"\lambda \textbf{I} ",r")","v","= 0").scale(2).next_to(total2,DOWN*2.5)
        
        A_term_eqn = emphasis_box(eigprob[1],ORANGE)
        A_term_mat = emphasis_box(total2[0][0],ORANGE)

        M_term_eqn = emphasis_box(eigprob[3],RED)
        M_term_mat = emphasis_box(total2[0][2],RED)

        vec_term_eqn = emphasis_box(eigprob[5],YELLOW)
        vec_term_mat = emphasis_box(total2[1][0],YELLOW)

        emphasis_boxes = VGroup(A_term_eqn[1],A_term_mat[1],M_term_eqn[1],M_term_mat[1],vec_term_eqn[1],vec_term_mat[1])
        Eigprobtext = Text("Solve The Eigen-Problem").to_edge(UP)
        

        
        title2 = VGroup(Eigprobtext,divider)
        eigprob_reformat = MathTex(r"(\textbf{A} - \omega^2 \textbf{I} ) \vec{a}", r"= 0").scale(2).next_to(title2,DOWN)
        
        eigprob_solve = MathTex(r"\text{det}",r"(\textbf{A} - \omega^2 \textbf{I} )", r"= 0").scale(2).next_to(eigprob_reformat,DOWN)


        if player == True:
            self.play(FadeOut(x_arrows),
                      FadeOut(masses_walls),
                      FadeOut(springs),
                      MoveToTarget(systemofeqn))
            self.wait(1)
            self.play(ReplacementTransform(systemofeqn,total))
            self.wait(1)
            self.play(FadeOut(total))
            return title2
        
        else:
            self.play(FadeOut(x_arrows),
                      FadeOut(masses_walls),
                      FadeOut(springs),
                      FadeOut(systemofeqn))
            
            self.play(FadeOut(title2))
            return title2
        
    def MathScene5(self,title,player=True):
        colors = self.colors
        text_soln = Text("Eign-Problem Solution").to_edge(UP)
        divider = Line([-4,0,0], [4,0,0]).next_to(text_soln,DOWN)
        title_soln = VGroup(text_soln,divider)
        final_eig_sys = self.solve_3DOF_problem(SPRING,MASS)
        eigprob_reformat = MathTex(r"(\textbf{A} - \omega^2 \textbf{I} ) \vec{a}", r"= 0").scale(1.5).next_to(title_soln,DOWN)

        
        
        Mode_Freqs = VGroup()
        previous = [eigprob_reformat]
        for ii in range(3):
            first_eigval = MathTex(r"\omega_{}^2 = ".format(ii+1),r"{}".format(final_eig_sys["eigvals"][ii])).next_to(previous[ii],DOWN*2).to_edge(LEFT).shift(RIGHT)
            text_mode_1 = MathTex(r"\vec{a}" + "_{}".format(ii+1)," = ")
            first_mode = Matrix([list(final_eig_sys["modes"][ii])],
                 element_alignment_corner=ORIGIN,
                 v_buff = 1.2,
                 h_buff = 1.3).next_to(text_mode_1,RIGHT).set_column_colors(colors[0],colors[1],colors[2])
            
            First_Mode = VGroup(text_mode_1,first_mode).next_to(first_eigval,RIGHT*4)
            
            first_freq = MathTex(r"f_{} = ".format(ii+1),r"{}".format(final_eig_sys["freqs"][ii]), r"\text{Hz}").next_to(First_Mode,RIGHT*4)

            First_Freq_Mode = VGroup(first_eigval,first_freq,First_Mode)
            Mode_Freqs.add(First_Freq_Mode)
            previous.append(first_eigval)
            
            
        
        Mode_Freqs[0].generate_target()
        Mode_Freqs[0].target.to_edge(DOWN)
        

        if player == True:
            self.play(ReplacementTransform(title,title_soln),Write(eigprob_reformat))
            self.play(Write(Mode_Freqs[0]))
            self.play(Write(Mode_Freqs[1]))
            self.play(Write(Mode_Freqs[2]))
            self.wait(2)
            
        else:
            None

        return Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys
    
    def solve_3DOF_problem(self,k,m):
        import numpy as np
        import scipy
        wo = np.sqrt(k/m)
        
        arr = np.array([[2*wo**2, - wo**2, 0],
                        [-wo**2, 2*wo**2, -wo**2],
                        [0, -wo**2, 2*wo**2]])
        
        [vals, modes] = scipy.linalg.eig(arr)
        idx = np.argsort(vals)
        
        
        
        mode_list = []
        freqs = []
        eigvals = []
        for ii in idx:
            mode_list.append(np.around(modes[:,ii],2))
            freqs.append(np.around(np.real(np.sqrt(vals[ii])/(2*np.pi)),2))
            eigvals.append(np.around(np.real(vals[ii]),2))
            
        final_eig_sys = dict({"eigvals": eigvals, "modes": mode_list, "freqs" : freqs})
        return final_eig_sys
        
    def ApplicationScene1(self,Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys,player=True):

        x_mass = []
        mass_wall = []
        spring_mode = []
        T0_vals = []
        
        x_mass_osc = []
        mass_wall_osc = []
        spring_mode_osc = []
        axes_mode_osc = []
        
        titles = []
        for ii in range(3):
            omega = 2*np.pi*final_eig_sys["freqs"][ii]
            mode_funcs = generate_mode_shape_funcs(final_eig_sys["modes"][ii],omega)
            mode_funcs_osc_graph = generate_mode_shape_funcs_osc_graph(final_eig_sys["modes"][ii],omega)
            axes,x_mass_osc_graph,masses_osc_graph, springs_osc_graph = setup_3DOF_Oscillation_graph(mode_funcs_osc_graph,omega,self.colors)
            x_mass_mode,masses_walls_mode, springs_mode = setup_3DOF_spring_mass_system(mode_funcs,omega,self.colors)
            
            
            x_mass_osc.append(x_mass_osc_graph)
            mass_wall_osc.append(masses_osc_graph)
            spring_mode_osc.append(springs_osc_graph)
            axes_mode_osc.append(axes)
            
            x_mass.append(x_mass_mode)
            mass_wall.append(masses_walls_mode)
            spring_mode.append(springs_mode)
            T0_vals.append(1/omega)
            title_txt = Text("Modal Visualization - Mode: {}".format(ii+1))
            divider = Line([-6,0,0], [6,0,0]).next_to(title_txt,DOWN)
            
            title = VGroup(title_txt,divider).scale(0.9).move_to(ORIGIN).to_edge(UP)
            titles.append(title)
            
            
        if player == True:
            
            self.play(FadeIn(mass_wall[0]),FadeIn(spring_mode[0]),
                      FadeIn(axes_mode_osc[0]),FadeIn(mass_wall_osc[0]),FadeIn(spring_mode_osc[0]),
                      Write(titles[0]))
            
            Mode_Freqs[1].to_edge(DOWN)
            Mode_Freqs[2].to_edge(DOWN)
            self.play(x_mass[0].animate(rate_func=rate_functions.linear, run_time=4).increment_value(16*T0_vals[0]),
                      x_mass_osc[0].animate(rate_func=rate_functions.linear, run_time=4).increment_value(16*T0_vals[0]))


        else:
            None

            
    
    def construct(self):
        self.colors = [RED,BLUE,GREEN]

        " Inputs "
        omega = 0.765
        T0    = 1/omega
        
        mode_funcs = generate_mode_shape_funcs([0,0,0],10)

        
        x_mass,masses_walls, springs = setup_3DOF_spring_mass_system(mode_funcs,omega,self.colors)

        self.add(*masses_walls, springs)
        
        eqnmotion,x_arrows = self.MathScene1(masses_walls, springs,True)
        systemofeqn = self.MathScene2(eqnmotion,x_arrows,masses_walls,False)
        title_eig_solve = self.MathScene3(systemofeqn,x_arrows,masses_walls,springs,True)
        Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys = self.MathScene5(title_eig_solve,False)
        self.ApplicationScene1(Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys,player=True)
        


class Multi_3DOF_Derivation(Scene):
    def MathScene1(self,masses_walls, springs, player=True):
        colors = self.colors
        eqnmotion = equation_generation(colors,1)
        
        x_arrows = VGroup()
        for ii in range(3):
            x_arrow_p1 = Line(masses_walls[0][ii].get_center(),masses_walls[0][ii].get_center() + [0,MASS_RADIUS/2,0]).shift(MASS_RADIUS*UP)
            x_arrow_p2 = Arrow(buff=0.5,start = x_arrow_p1.get_end(),
                                end = x_arrow_p1.get_end() + [MASS_RADIUS/2,0,0],max_stroke_width_to_length_ratio=15)
            x_arrow_txt = MathTex(r"x_{}".format(ii+1)).scale(0.8).next_to(x_arrow_p2,RIGHT/2)
            x_arrows.add( VGroup(x_arrow_p1,x_arrow_p2,x_arrow_txt).set_color(colors[ii]))
        
        if player == True:
            self.play(*[FadeIn(x_arrows[ii]) for ii in range(3)])
            for ii in range(3):
                self.wait(1)
                self.play(Indicate(masses_walls[0][ii]),
                          Indicate(x_arrows[ii]))
                self.play(Write(eqnmotion[ii][0]))
                self.wait(1)
                self.play(Indicate(masses_walls[0][ii]),
                          Indicate(springs[ii]),
                          Indicate(springs[ii+1]))
                self.play(Write(eqnmotion[ii][1]))

            return eqnmotion,x_arrows
        
        else:
            self.add(eqnmotion,x_arrows)
            return eqnmotion,x_arrows

    def MathScene2(self,eqnmotion,x_arrows,masses_walls,player = True):
        colors = self.colors
        eqnmotion.generate_target()
        eqnmotion.target.shift(LEFT*2.5)
        
        "Equation of motion set 2"
        eqnmotion_g2 = equation_generation(colors,2).shift(LEFT*2.5)
    
        "Equation of motion set 3"
        eqnmotion_g3 = equation_generation(colors,3).shift(LEFT*2.5)
        
        "Equation of motion set 4"
        eqnmotion_g4 = equation_generation(colors,4).shift(LEFT*2.5)
        
        
        "Harmonic Equations"
        harm1 = MathTex(r"x_1",r"=", r"a_1"," e^{j\omega t").shift(DOWN).shift(RIGHT*4)
        harm1[0].set_color(self.colors[0])
        harm1[2].set_color(self.colors[0])

        harm2 = MathTex(r"x_2",r"=", r"a_2", "e^{j\omega t").next_to(harm1,DOWN)
        harm2[0].set_color(self.colors[1])
        harm2[2].set_color(self.colors[1])

        harm3 = MathTex(r"x_3",r"=", r"a_3", " e^{j\omega t").next_to(harm2,DOWN)
        harm3[0].set_color(self.colors[2])
        harm3[2].set_color(self.colors[2])

        harmonic = VGroup(harm1,harm2,harm3)

        cancel_t1 = cancel_line(eqnmotion_g2[0][0][6:10],ORANGE)
        cancel_t2 = cancel_line(eqnmotion_g2[0][1][16:20],ORANGE)
        cancel_eqn1 = VGroup(cancel_t1,cancel_t2)
    
        cancel_t1 = cancel_line(eqnmotion_g2[1][0][6:10],ORANGE)
        cancel_t2 = cancel_line(eqnmotion_g2[1][1][21:25],ORANGE)
        cancel_eqn2 = VGroup(cancel_t1,cancel_t2)
        
        cancel_t1 = cancel_line(eqnmotion_g2[2][0][6:10],ORANGE)
        cancel_t2 = cancel_line(eqnmotion_g2[2][1][16:20],ORANGE)
        cancel_eqn3 = VGroup(cancel_t1,cancel_t2)
        
        cancel_lines = VGroup(cancel_eqn1,cancel_eqn2,cancel_eqn3)

        
        phrase = Text("Let:").scale(0.7).shift(DOWN).shift(RIGHT*4.5)
        state = MathTex(r"w_o",r"=\sqrt{\frac{k}{m}}").next_to(phrase,DOWN)
        wo_condition = VGroup(phrase,state)

            
        eqn_last = finaleqnmotion(colors).shift(LEFT*2.5)

        A_mat = Matrix([["-\omega^2 + 2 \omega_o^2", "- \omega_o^2            ", "0"],
                              ["-\omega_o^2             ", "-\omega^2 + 2 \omega_o^2", " - \omega_o^2"],
                              ["0"                       , "-\omega_o^2             ", "-\omega^2 + 2 \omega_o^2"]],
                             element_alignment_corner=ORIGIN,
                             v_buff = 1.2,
                             h_buff = 2.5)
        b_mat = Matrix([["a_1"],["a_2"],["a_3"]], 
                       element_alignment_corner=ORIGIN,v_buff = 1.2,).next_to(A_mat,RIGHT).set_row_colors(colors[0],colors[1],colors[2])
                        
        equals = MathTex(r"=").next_to(b_mat,RIGHT)
        f_mat = Matrix([["0"],
                      ["0"],
                      ["0"]], 
               element_alignment_corner=ORIGIN,v_buff = 1.2,).next_to(equals,RIGHT)
                
        
                        
        systemofeqn = VGroup(A_mat,b_mat,equals,f_mat).move_to(ORIGIN).shift(2*DOWN)
        
        if player == True:
            self.play(MoveToTarget(eqnmotion),runtime=2)
            self.wait(2)
            self.play(Write(harmonic))
            self.wait(1)
            self.play(Circumscribe(harmonic))
            self.wait(1)
            self.play(Circumscribe(eqnmotion))
            self.wait(1)
            self.play(FadeOut(harmonic,target_position=eqnmotion),
                      ReplacementTransform(eqnmotion,eqnmotion_g2))
            self.wait(1)
            self.play(Write(cancel_lines))
            self.wait(1)
            self.play(FadeOut(cancel_lines),
                  TransformMatchingTex(eqnmotion_g2,eqnmotion_g3))
            self.wait(1)
            self.play(Write(wo_condition))
            self.wait(2)
            self.play(TransformMatchingTex(eqnmotion_g3,eqnmotion_g4))
            self.wait(2)
            self.play(TransformMatchingTex(eqnmotion_g4,eqn_last))
            self.wait(3)
            self.play(FadeOut(wo_condition),ReplacementTransform(eqn_last,systemofeqn))

            return systemofeqn
        else:
            self.play(ReplacementTransform(eqnmotion,systemofeqn))
            return systemofeqn
        
    def MathScene3(self,systemofeqn,x_arrows,masses_walls,springs,player = True):
        colors = self.colors
        systemofeqn.generate_target()
        systemofeqn.target.move_to(ORIGIN)
        
        K_mat = Matrix([["2 \omega_o^2", "- \omega_o^2            ", "0"],
                      ["-\omega_o^2             ", " 2 \omega_o^2", " - \omega_o^2"],
                      ["0"                       , "-\omega_o^2             ", "2 \omega_o^2"]],
                     element_alignment_corner=ORIGIN,
                     v_buff = 1.2,
                     h_buff = 1.4)
        
        subsign = MathTex("-").next_to(K_mat,RIGHT)
        
        M_mat = Matrix([["\omega^2", "0", "0"],
              ["0", "\omega^2", "0"],
              ["0" , "0", "\omega^2"]],
             element_alignment_corner=ORIGIN,
             v_buff = 1.2,
             h_buff = 1.0).next_to(subsign,RIGHT)
        

        sec1 = VGroup(K_mat, subsign, M_mat)
        
        b_mat = Matrix([["a_1"],["a_2"],["a_3"]], 
           element_alignment_corner=ORIGIN,v_buff = 1.2,).set_row_colors(colors[0],colors[1],colors[2])
                
        
        equals = MathTex(r"=").next_to(b_mat,RIGHT)
        f_mat = Matrix([["0"],
                      ["0"],
                      ["0"]], 
               element_alignment_corner=ORIGIN,v_buff = 1.2,).next_to(equals,RIGHT)
                
        sec2 = VGroup(b_mat, equals, f_mat).next_to(sec1,RIGHT)
        
        total = VGroup(sec1,sec2).move_to(ORIGIN)
        
        
        eigval = MathTex(r"\omega^2").next_to(sec1[1],RIGHT)
        M_mat = Matrix([["1", "0", "0"],
              ["0", "1", "0"],
              ["0" , "0", "1"]],
             element_alignment_corner=ORIGIN,
             v_buff = 1.2,
             h_buff = 1.0).next_to(eigval,RIGHT)
        Newmat = VGroup(eigval,M_mat)
        
        sec1 = VGroup(K_mat, subsign, Newmat)
        sec2 = total[1].next_to(Newmat,RIGHT)
        total2 = VGroup(sec1,sec2)

        Eigprobtext = Text("Eigen-Problem").to_edge(UP)
        divider = Line([-4,0,0], [4,0,0]).next_to(Eigprobtext,DOWN)
        title = VGroup(Eigprobtext,divider)
        
        eigprob = MathTex(r"(",r"\textbf{A}","-",r"\lambda \textbf{I} ",r")","v","= 0").scale(2).next_to(total2,DOWN*2.5)
        
        A_term_eqn = emphasis_box(eigprob[1],ORANGE)
        A_term_mat = emphasis_box(total2[0][0],ORANGE)

        M_term_eqn = emphasis_box(eigprob[3],RED)
        M_term_mat = emphasis_box(total2[0][2],RED)

        vec_term_eqn = emphasis_box(eigprob[5],YELLOW)
        vec_term_mat = emphasis_box(total2[1][0],YELLOW)

        emphasis_boxes = VGroup(A_term_eqn[1],A_term_mat[1],M_term_eqn[1],M_term_mat[1],vec_term_eqn[1],vec_term_mat[1])
        Eigprobtext = Text("Solve The Eigen-Problem").to_edge(UP)
        

        
        title2 = VGroup(Eigprobtext,divider)
        eigprob_reformat = MathTex(r"(\textbf{A} - \omega^2 \textbf{I} ) \vec{a}", r"= 0").scale(2).next_to(title2,DOWN)
        
        eigprob_solve = MathTex(r"\text{det}",r"(\textbf{A} - \omega^2 \textbf{I} )", r"= 0").scale(2).next_to(eigprob_reformat,DOWN)


        if player == True:
            self.play(FadeOut(x_arrows),
                      FadeOut(masses_walls),
                      FadeOut(springs),
                      MoveToTarget(systemofeqn))
            self.wait(1)
            self.play(ReplacementTransform(systemofeqn,total))
            self.wait(2)
            self.play(ReplacementTransform(total,total2))
            self.wait(1)
            self.play(Write(eigprob),Write(title))
            self.wait(1)
            self.play(Create(A_term_eqn[1]),Create(A_term_mat[1]))
            self.wait(1)
            self.play(Create(M_term_eqn[1]),Create(M_term_mat[1]))
            self.wait(1)
            self.play(Create(vec_term_eqn[1]),Create(vec_term_mat[1]))
            self.wait(1)
            self.play(ReplacementTransform(title,title2),
                      TransformMatchingTex(eigprob,eigprob_reformat),
                      FadeOut(emphasis_boxes),
                      FadeOut(total2),
                      FadeOut(Newmat),
                      FadeOut(total))
            self.wait(1)
            self.play(Write(eigprob_solve))
            self.wait(2)
            self.play(FadeOut(eigprob_solve),
                      FadeOut(eigprob_reformat))
            self.wait(1)
            
            return title2
        
        else:
            self.play(FadeOut(x_arrows),
                      FadeOut(masses_walls),
                      FadeOut(springs),
                      FadeOut(systemofeqn))
            
            self.play(FadeOut(title2))
            return title2
        
    def MathScene5(self,title,player=True):
        colors = self.colors
        text_soln = Text("Eign-Problem Solution").to_edge(UP)
        divider = Line([-4,0,0], [4,0,0]).next_to(text_soln,DOWN)
        title_soln = VGroup(text_soln,divider)
        final_eig_sys = self.solve_3DOF_problem(SPRING,MASS)
        eigprob_reformat = MathTex(r"(\textbf{A} - \omega^2 \textbf{I} ) \vec{a}", r"= 0").scale(1.5).next_to(title_soln,DOWN)

        
        
        Mode_Freqs = VGroup()
        previous = [eigprob_reformat]
        for ii in range(3):
            first_eigval = MathTex(r"\omega_{}^2 = ".format(ii+1),r"{}".format(final_eig_sys["eigvals"][ii])).next_to(previous[ii],DOWN*2).to_edge(LEFT).shift(RIGHT)
            text_mode_1 = MathTex(r"\vec{a}" + "_{}".format(ii+1)," = ")
            first_mode = Matrix([list(final_eig_sys["modes"][ii])],
                 element_alignment_corner=ORIGIN,
                 v_buff = 1.2,
                 h_buff = 1.3).next_to(text_mode_1,RIGHT).set_column_colors(colors[0],colors[1],colors[2])
            
            First_Mode = VGroup(text_mode_1,first_mode).next_to(first_eigval,RIGHT*4)
            
            first_freq = MathTex(r"f_{} = ".format(ii+1),r"{}".format(final_eig_sys["freqs"][ii]), r"\text{Hz}").next_to(First_Mode,RIGHT*4)

            First_Freq_Mode = VGroup(first_eigval,first_freq,First_Mode)
            Mode_Freqs.add(First_Freq_Mode)
            previous.append(first_eigval)
            
            
        
        Mode_Freqs[0].generate_target()
        Mode_Freqs[0].target.to_edge(DOWN)
        

        if player == True:
            self.play(ReplacementTransform(title,title_soln),Write(eigprob_reformat))
            self.play(Write(Mode_Freqs[0]))
            self.play(Write(Mode_Freqs[1]))
            self.play(Write(Mode_Freqs[2]))
            self.wait(2)
            
        else:
            None

        return Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys
    
    def solve_3DOF_problem(self,k,m):
        import numpy as np
        import scipy
        wo = np.sqrt(k/m)
        
        arr = np.array([[2*wo**2, - wo**2, 0],
                        [-wo**2, 2*wo**2, -wo**2],
                        [0, -wo**2, 2*wo**2]])
        
        [vals, modes] = scipy.linalg.eig(arr)
        idx = np.argsort(vals)
        
        
        
        mode_list = []
        freqs = []
        eigvals = []
        for ii in idx:
            mode_list.append(np.around(modes[:,ii],2))
            freqs.append(np.around(np.real(np.sqrt(vals[ii])/(2*np.pi)),2))
            eigvals.append(np.around(np.real(vals[ii]),2))
            
        final_eig_sys = dict({"eigvals": eigvals, "modes": mode_list, "freqs" : freqs})
        return final_eig_sys
        
    def ApplicationScene1(self,Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys,player=True):

        x_mass = []
        mass_wall = []
        spring_mode = []
        T0_vals = []
        
        x_mass_osc = []
        mass_wall_osc = []
        spring_mode_osc = []
        axes_mode_osc = []
        
        titles = []
        for ii in range(3):
            omega = 2*np.pi*final_eig_sys["freqs"][ii]
            mode_funcs = generate_mode_shape_funcs(final_eig_sys["modes"][ii],omega)
            mode_funcs_osc_graph = generate_mode_shape_funcs_osc_graph(final_eig_sys["modes"][ii],omega)
            axes,x_mass_osc_graph,masses_osc_graph, springs_osc_graph = setup_3DOF_Oscillation_graph(mode_funcs_osc_graph,omega,self.colors)
            x_mass_mode,masses_walls_mode, springs_mode = setup_3DOF_spring_mass_system(mode_funcs,omega,self.colors)
            
            
            x_mass_osc.append(x_mass_osc_graph)
            mass_wall_osc.append(masses_osc_graph)
            spring_mode_osc.append(springs_osc_graph)
            axes_mode_osc.append(axes)
            
            x_mass.append(x_mass_mode)
            mass_wall.append(masses_walls_mode)
            spring_mode.append(springs_mode)
            T0_vals.append(1/omega)
            title_txt = Text("Modal Visualization - Mode: {}".format(ii+1))
            divider = Line([-6,0,0], [6,0,0]).next_to(title_txt,DOWN)
            
            title = VGroup(title_txt,divider).scale(0.9).move_to(ORIGIN).to_edge(UP)
            titles.append(title)
            
            
        if player == True:
            
            self.play(FadeOut(Mode_Freqs[1]), FadeOut(Mode_Freqs[2]), FadeOut(eigprob_reformat),
                      FadeOut(title_soln), MoveToTarget(Mode_Freqs[0]),
                      FadeIn(mass_wall[0]),FadeIn(spring_mode[0]),
                      FadeIn(axes_mode_osc[0]),FadeIn(mass_wall_osc[0]),FadeIn(spring_mode_osc[0]),
                      Write(titles[0]))
            
            Mode_Freqs[1].to_edge(DOWN)
            Mode_Freqs[2].to_edge(DOWN)
            self.play(x_mass[0].animate(rate_func=rate_functions.linear, run_time=6).increment_value(12*T0_vals[0]),
                      x_mass_osc[0].animate(rate_func=rate_functions.linear, run_time=6).increment_value(12*T0_vals[0]))
            
            self.wait(1)
            for ii in range(2):
                self.play(ReplacementTransform(Mode_Freqs[ii],Mode_Freqs[ii+1]),
                       ReplacementTransform(mass_wall[ii],mass_wall[ii+1]),
                      ReplacementTransform(spring_mode[ii],spring_mode[ii+1]),
                      ReplacementTransform(axes_mode_osc[ii],axes_mode_osc[ii+1]),
                       ReplacementTransform(mass_wall_osc[ii],mass_wall_osc[ii+1]),
                      ReplacementTransform(spring_mode_osc[ii],spring_mode_osc[ii+1]),
                      ReplacementTransform(titles[ii],titles[ii+1]))
                
                self.wait(1)
                self.play(x_mass[ii+1].animate(rate_func=rate_functions.linear, run_time=6).increment_value(12*T0_vals[0]),
                      x_mass_osc[ii+1].animate(rate_func=rate_functions.linear, run_time=6).increment_value(12*T0_vals[0]))
                self.wait(1)

        else:
            None

            
    
    def construct(self):
        self.colors = [RED,BLUE,GREEN]

        " Inputs "
        omega = 0.765
        T0    = 1/omega
        
        mode_funcs = generate_mode_shape_funcs([0,0,0],10)

        
        x_mass,masses_walls, springs = setup_3DOF_spring_mass_system(mode_funcs,omega,self.colors)

        self.add(*masses_walls, springs)
        
        eqnmotion,x_arrows = self.MathScene1(masses_walls, springs,False)
        self.wait(2)
        systemofeqn = self.MathScene2(eqnmotion,x_arrows,masses_walls,True)
        self.wait(2)
        title_eig_solve = self.MathScene3(systemofeqn,x_arrows,masses_walls,springs,True)
        self.wait(2)
        Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys = self.MathScene5(title_eig_solve,False)
        self.wait(2)

        self.ApplicationScene1(Mode_Freqs,eigprob_reformat,title_soln,final_eig_sys,player=False)
        


class Multi_DOF_Short(Scene):
    def solve_3DOF_problem(self,k,m):
        import numpy as np
        import scipy
        wo = np.sqrt(k/m)
        
        arr = np.array([[2*wo**2, - wo**2, 0],
                        [-wo**2, 2*wo**2, -wo**2],
                        [0, -wo**2, 2*wo**2]])
        
        [vals, modes] = scipy.linalg.eig(arr)
        idx = np.argsort(vals)
        
        
        
        mode_list = []
        freqs = []
        eigvals = []
        for ii in idx:
            mode_list.append(np.around(modes[:,ii],2))
            freqs.append(np.around(np.real(np.sqrt(vals[ii])/(2*np.pi)),2))
            eigvals.append(np.around(np.real(vals[ii]),2))
            
        final_eig_sys = dict({"eigvals": eigvals, "modes": mode_list, "freqs" : freqs})
        return final_eig_sys
        
    def ApplicationScene1(self,final_eig_sys,player=True):
        colors = self.colors
        x_mass = []
        mass_wall = []
        spring_mode = []
        T0_vals = []
        
        x_mass_osc = []
        mass_wall_osc = []
        spring_mode_osc = []
        axes_mode_osc = []
        text_mode_1 = MathTex(r"\vec{a}" + "_{}".format(1)," = ")
        first_mode = Matrix([list(final_eig_sys["modes"][0])],
             element_alignment_corner=ORIGIN,
             v_buff = 1.2,
             h_buff = 1.3).next_to(text_mode_1,RIGHT).set_column_colors(colors[0],colors[1],colors[2])
        
        First_Mode = VGroup(text_mode_1,first_mode).move_to(ORIGIN).to_edge(DOWN)
        
        for ii in range(3):
            omega = 2*np.pi*final_eig_sys["freqs"][ii]
            mode_funcs = generate_mode_shape_funcs(final_eig_sys["modes"][ii],omega)
            mode_funcs_osc_graph = generate_mode_shape_funcs_osc_graph(final_eig_sys["modes"][ii],omega)
            axes,x_mass_osc_graph,masses_osc_graph, springs_osc_graph = setup_3DOF_Oscillation_graph(mode_funcs_osc_graph,omega,colors)
            x_mass_mode,masses_walls_mode, springs_mode = setup_3DOF_spring_mass_system(mode_funcs,omega,colors)
            
            
            x_mass_osc.append(x_mass_osc_graph)
            mass_wall_osc.append(masses_osc_graph)
            spring_mode_osc.append(springs_osc_graph)
            axes_mode_osc.append(axes)
            
            x_mass.append(x_mass_mode)
            mass_wall.append(masses_walls_mode)
            spring_mode.append(springs_mode)
            T0_vals.append(1/omega)

    

        
        if player == True:
            self.add(mass_wall[0],spring_mode[0],axes_mode_osc[0],mass_wall_osc[0],spring_mode_osc[0],First_Mode)
            #self.play(x_mass[0].animate(rate_func=rate_functions.linear, run_time=10).increment_value(20*T0_vals[0]),
            #          x_mass_osc[0].animate(rate_func=rate_functions.linear, run_time=10).increment_value(20*T0_vals[0]))
            #self.wait(1)

    def ImgScene1(self,final_eig_sys,player=True):
        colors = self.colors
        x_mass = []
        mass_wall = []
        spring_mode = []
        T0_vals = []
        
        x_mass_osc = []
        mass_wall_osc = []
        spring_mode_osc = []
        axes_mode_osc = []
        text_mode_1 = MathTex(r"\vec{a}" + "_{}".format(1)," = ")
        first_mode = Matrix([list(final_eig_sys["modes"][0])],
             element_alignment_corner=ORIGIN,
             v_buff = 1.2,
             h_buff = 1.3).next_to(text_mode_1,RIGHT).set_column_colors(colors[0],colors[1],colors[2])
        
        First_Mode = VGroup(text_mode_1,first_mode).move_to(ORIGIN).to_edge(DOWN)
        final_eig_sys["modes"][0] = [0]*3

        for ii in range(3):
            omega = 2*np.pi*final_eig_sys["freqs"][ii]
            mode_funcs = generate_mode_shape_funcs(final_eig_sys["modes"][ii],omega)
            mode_funcs_osc_graph = generate_mode_shape_funcs_osc_graph(final_eig_sys["modes"][ii],omega)
            axes,x_mass_osc_graph,masses_osc_graph, springs_osc_graph = setup_3DOF_Oscillation_graph(mode_funcs_osc_graph,omega,colors)
            x_mass_mode,masses_walls_mode, springs_mode = setup_3DOF_spring_mass_system(mode_funcs,omega,colors)
            
            
            x_mass_osc.append(x_mass_osc_graph)
            mass_wall_osc.append(masses_osc_graph)
            spring_mode_osc.append(springs_osc_graph)
            axes_mode_osc.append(axes)
            
            x_mass.append(x_mass_mode)
            mass_wall.append(masses_walls_mode)
            spring_mode.append(springs_mode)
            T0_vals.append(1/omega)

    
        eqn_last = finaleqnmotion(colors)

        
        if player == True:
            self.add(mass_wall[0],spring_mode[0],eqn_last)
            #self.add(mass_wall[0].move_to(ORIGIN),spring_mode[0].move_to(ORIGIN))
            #self.play(x_mass[0].animate(rate_func=rate_functions.linear, run_time=10).increment_value(20*T0_vals[0]),
            #          x_mass_osc[0].animate(rate_func=rate_functions.linear, run_time=10).increment_value(20*T0_vals[0]))
            #self.wait(1)


    def construct(self):
        self.colors = [RED,BLUE,GREEN]
        k = 1
        m = 1
        final_eig_sys = self.solve_3DOF_problem(k,m)
        #self.ApplicationScene1(final_eig_sys,player=True)
        self.ImgScene1(final_eig_sys,player=True)

class Multi_DOF_LongForm(Scene):
    def solve_3DOF_problem(self,k,m,n=6):
        import numpy as np
        import scipy
        wo = np.sqrt(k/m)
        from scipy.sparse import diags
        import numpy as np
        
        k = [-wo**2*np.ones(n-1),2*wo**2*np.ones(n),- wo**2*np.ones(n-1)]
        offset = [-1,0,1]
        arr = diags(k,offset).toarray()
        [vals, modes] = scipy.linalg.eig(arr)
        idx = np.argsort(vals)
        
        
        
        mode_list = []
        freqs = []
        eigvals = []
        for ii in idx:
            mode_list.append(np.around(modes[:,ii],2))
            freqs.append(np.around(np.real(np.sqrt(vals[ii])/(2*np.pi)),2))
            eigvals.append(np.around(np.real(vals[ii]),2))
            
        final_eig_sys = dict({"eigvals": eigvals, "modes": mode_list, "freqs" : freqs})
        return final_eig_sys
        
    def ApplicationScene1(self,player=True):
        colors = self.colors
        x_mass = []
        mass_wall = []
        spring_mode = []
        T0_vals = []
        
        x_mass_osc = []
        mass_wall_osc = []
        spring_mode_osc = []
        axes_mode_osc = []
        mode_choice = 2
        base = 3
        n = 12
        high_DOF_Soln = self.solve_3DOF_problem(1,1,50)
        omega_high = 2*np.pi*high_DOF_Soln["freqs"][mode_choice]
        print(omega_high)
        mode_shape_High_DOF = high_DOF_Soln["modes"][mode_choice]/max(abs(high_DOF_Soln["modes"][mode_choice]))
        print(mode_shape_High_DOF)
        mode_funcs_High_DOF = generate_mode_shape_funcs_osc_graph(mode_shape_High_DOF,omega_high)
        axes_high,x_mass_osc_graph_high,masses_osc_graph_high, springs_osc_graph_high = Infinite_String(mode_funcs_High_DOF,omega_high)

        for ii in range(base,n):
            final_eig_sys = self.solve_3DOF_problem(1,1,ii)
            mode_shape = final_eig_sys["modes"][mode_choice]/max(abs(final_eig_sys["modes"][mode_choice]))
            
            mode_funcs = generate_mode_shape_funcs(mode_shape,omega_high)
            mode_funcs_osc_graph = generate_mode_shape_funcs_osc_graph(mode_shape,omega_high)
            axes,x_mass_osc_graph,masses_osc_graph, springs_osc_graph = setup_3DOF_Oscillation_graph(mode_funcs_osc_graph,omega_high,mass_radius=MASS_RADIUS_OSC)

            x_mass_mode,masses_walls_mode, springs_mode = setup_3DOF_spring_mass_system(mode_funcs,omega_high,mass_radius=MASS_RADIUS)

            
            x_mass_osc.append(x_mass_osc_graph)
            mass_wall_osc.append(masses_osc_graph)
            spring_mode_osc.append(springs_osc_graph)
            axes_mode_osc.append(axes)
            
            x_mass.append(x_mass_mode)
            mass_wall.append(masses_walls_mode)
            spring_mode.append(springs_mode)

            T0_vals.append(1/(omega_high/(2*np.pi)))


        self.add(mass_wall[0],spring_mode[0],axes_mode_osc[0],mass_wall_osc[0],spring_mode_osc[0])
        
        shape_num = n-base
        if player == True:
            for ii in range(shape_num-2):
                self.play(x_mass_osc[ii].animate(rate_func=rate_functions.linear,run_time=3).increment_value(T0_vals[ii]*3),
                          x_mass[ii].animate(rate_func=rate_functions.linear,run_time=3).increment_value(T0_vals[ii]*3))

                self.play(ReplacementTransform(axes_mode_osc[ii],axes_mode_osc[ii+1]),
                          ReplacementTransform(mass_wall_osc[ii],mass_wall_osc[ii+1]),
                          ReplacementTransform(spring_mode_osc[ii],spring_mode_osc[ii+1]),
                          ReplacementTransform(mass_wall[ii],mass_wall[ii+1]),
                          ReplacementTransform(spring_mode[ii],spring_mode[ii+1]))
                
            ii = ii + 1
            self.play(x_mass_osc[ii].animate(rate_func=rate_functions.linear,run_time=3).increment_value(T0_vals[ii]*3),
                      x_mass[ii].animate(rate_func=rate_functions.linear,run_time=3).increment_value(T0_vals[ii]*3))

            self.play(ReplacementTransform(mass_wall_osc[ii],masses_osc_graph_high),
                ReplacementTransform(spring_mode_osc[ii],springs_osc_graph_high),
                ReplacementTransform(mass_wall[ii],mass_wall[ii+1]),
                ReplacementTransform(spring_mode[ii],spring_mode[ii+1]))
                    
            self.play(x_mass_osc_graph_high.animate(rate_func=rate_functions.linear,run_time=5).increment_value(T0_vals[ii]*5),
                  x_mass[ii+1].animate(rate_func=rate_functions.linear,run_time=5).increment_value(T0_vals[ii]*5))

    def construct(self):
        self.colors = [RED,BLUE,GREEN]
        self.ApplicationScene1(player=True)

