# -*- coding: utf-8 -*-
"""
Created on Mon May 15 17:40:36 2023

@author: mbluu
"""

from manim import *

INITIAL_SHIFT = 2


def cancel_line(text,color):
    pre_coord_dl = text.get_corner(DL)
    pre_coord_ur = text.get_corner(UR)
    reference_line = Line(pre_coord_dl,pre_coord_ur,color=color)
    return reference_line


        
        
class Galerkin_Method(Scene):
    def create_bar_up(self):
        rect = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [3.0,-0.5,0], [3.0,0.5,0],
            fill_color = BLUE, fill_opacity=0.5
        ) # Stationary bar movement

        l  = Line([-3., 1., 0.], [-3., -1., 0.])
        l2 = Line([-3., 1., 0.], [-4., 0.5, 0.])
        l3 = Line([-3., 0.5, 0.], [-4., 0.0, 0.])
        l4 = Line([-3., 0., 0.], [-4., -0.5, 0.])
        l5 = Line([-3., -0.5, 0.], [-4., -1.0, 0.])
        l6 = Line([-3., -1., 0.], [-4., -1.5, 0.])
        BC = VGroup(l,l2,l3,l4,l5,l6)
        
        Bar = VGroup(rect,BC).shift(UP*2)     
        return Bar
    
    def Create_Parsel1(self):
        # Describe Parsel
        Parsel1 = Rectangle(height=1, width=1, stroke_width=2,color=WHITE,fill_opacity=1.0).move_to([2.5,0,0]).shift(INITIAL_SHIFT*UP)
        emphasisline1 = Line([2.0,-0.5,0],[2.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        emphasisline2 = Line([3.0,-0.5,0],[3.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        arrow_wave1 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([1.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        arrow_wave2 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([2.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        xi1 = MathTex(r"\xi",color=BLUE).next_to(arrow_wave1,UP/2)
        xi2 = MathTex(r"\xi + d\xi",color=BLUE).next_to(arrow_wave2,UR/2,aligned_edge=LEFT)
        brace = BraceBetweenPoints([2.0,-0.5,0],[3.0,-0.5,0],color=GREY).shift(INITIAL_SHIFT*UP)
        dx = Text("dx").next_to(brace,DOWN).set_color(GREY)
        Parsel_Group = Group(Parsel1,brace,dx,emphasisline1,emphasisline2,arrow_wave1,arrow_wave2,xi1,xi2)

        return Parsel_Group
    
    def Create_Parsel2(self):
        # Describe Parsel
        Parsel1 = Rectangle(height=1, width=2, stroke_width=2,color=WHITE,fill_opacity=1.0).move_to([3.0,0,0]).shift(INITIAL_SHIFT*UP)
        emphasisline1 = Line([2.0,-0.5,0],[2.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        emphasisline2 = Line([4.0,-0.5,0],[4.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        arrow_wave1 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([1.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        arrow_wave2 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([3.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        xi1 = MathTex(r"\xi",color=BLUE).next_to(arrow_wave1,UP/2)
        xi2 = MathTex(r"\xi + d\xi",color=BLUE).next_to(arrow_wave2,UR/2,aligned_edge=LEFT)
        brace = BraceBetweenPoints([2.0,-0.5,0],[4.0,-0.5,0],color=GREY).shift(INITIAL_SHIFT*UP)
        dx = Text("dx").next_to(brace,DOWN).set_color(GREY)
        
        bracetot = Group(brace,dx)
        displacetot = Group(xi2, arrow_wave2,emphasisline2)
        
        Parsel_Group = Group(Parsel1,brace,dx,emphasisline1,emphasisline2,arrow_wave1,arrow_wave2,xi1,xi2).move_to([4.0,-1.0,0])
        fadeoutgroup = Group(xi1,Parsel1,emphasisline1,arrow_wave1)
        
        return Parsel_Group,displacetot,bracetot,fadeoutgroup
    
    def DemonstrateLongWaveMotionScene1(self,player=True):
        # Mobject Creation Beam Bar Physical Creation
        width = ValueTracker(3.0)

        rect1 = always_redraw(lambda : Polygon(
            [-3,0.5,0], [-3,-0.5,0], [width.get_value(),-0.5,0], [width.get_value(),0.5,0],
            fill_color = BLUE, fill_opacity=0.5
        )) # Longitudinal movement

        
        rectfake = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [2.0,-0.5,0], [2.0,0.5,0],
            fill_color = BLUE, fill_opacity=0.5
        ) # Stationary bar movement
   
        rectfake_forward = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [2.0,-0.5,0], [2.0,0.5,0],
            fill_color = BLUE, fill_opacity=0.1
        ).shift(INITIAL_SHIFT*UP) # Stationary bar movement
        
        rectfake_forward_expand = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [3.0,-0.5,0], [3.0,0.5,0],
            fill_color = RED, fill_opacity=0.5
        ).shift(INITIAL_SHIFT*UP)
        # Boundary Condition fixed
        l  = Line([-3., 1., 0.], [-3., -1., 0.])
        l2 = Line([-3., 1., 0.], [-4., 0.5, 0.])
        l3 = Line([-3., 0.5, 0.], [-4., 0.0, 0.])
        l4 = Line([-3., 0., 0.], [-4., -0.5, 0.])
        l5 = Line([-3., -0.5, 0.], [-4., -1.0, 0.])
        l6 = Line([-3., -1., 0.], [-4., -1.5, 0.])
        BC = Group(l,l2,l3,l4,l5,l6)

        # Forcing function
        arrow_1 = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(rect1,RIGHT)
        Fs = Text("F",color=RED).next_to(arrow_1,RIGHT)
        Force = Group(arrow_1,Fs)
        Bar = Group(BC,rect1)


        # Differential Equations
        MainDiff = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = F")
        MainDiff[3].set_color(RED)
        for ii in range(3):
            MainDiff[ii].set_color(BLUE)

        Bar1 = Group(BC,rectfake)


        # Animation Portion
        if player == True:
            self.add(Bar)
            self.wait(1)
            self.play(FadeIn(Force))
            self.wait(1)
            self.play(FadeOut(Force))
            
            length_of_bar = [2.0,4.0,2.0,4.0,2.0]
            for l_bar in length_of_bar:
                self.play(width.animate.set_value(l_bar),run_time=1.0)

    
            self.play(ReplacementTransform(Bar,Bar1))
            self.play(Bar1.animate.shift(INITIAL_SHIFT*UP))
            
            Question = Text("Can we create an equation to describe this motion? ").scale(0.8)
            NecessaryEquations = Text("Necessary Equations").scale(1)
            ul = Underline(NecessaryEquations)  # Underlining the word
            NecessaryEquations = VGroup(NecessaryEquations,ul)       
            
            
            
            Eqn1 = Text(" 1. ???").scale(1).next_to(NecessaryEquations,DOWN)
            Eqn1_real = Text(" 1. Elastic Motion").scale(1).next_to(NecessaryEquations,DOWN)
            box1_real = SurroundingRectangle(Eqn1_real, corner_radius=0.2)
            Eqn1_real_group = Group(Eqn1_real,box1_real)
            
            Eqn2 = Text("2. ???").scale(1).next_to(Eqn1,DOWN*2)
            Eqn3 = Text("3. ???").scale(1).next_to(Eqn2,DOWN*2)
            
            Eqns = Group(Eqn1,Eqn2,Eqn3)
            box1 = SurroundingRectangle(Eqn1, corner_radius=0.2)
            Eqn1_box = Group(Eqn1,box1)
            self.wait(1)
            self.play(FadeIn(Question))
            self.wait(2)
            self.play(Transform(Question,NecessaryEquations))
            self.wait(1)
            self.play(FadeIn(Eqns))
            self.wait(2)
            self.play(Create(box1))
            self.wait(1)
            self.play(ReplacementTransform(Eqn1_box,Eqn1_real_group))
            Eqns = Group(Eqn1_real_group,Eqn2,Eqn3)
    
            self.wait(3)
            self.play(FadeOut(Question,NecessaryEquations,Eqns,box1,Eqn1_real_group,Eqn1_real))
        else:
            None
            
        return Bar1,BC,rectfake
    

    
    def HookesLawScene2(self,Bar1,BC,rectfake,player=True):
        # Mobject Creation Beam Bar Physical Creation
        width = ValueTracker(3.0)
              
        # Set up problem
        rectfake_forward = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [2.0,-0.5,0], [2.0,0.5,0],
            fill_color = BLUE, fill_opacity=0.1
        ).shift(INITIAL_SHIFT*UP) # Stationary bar movement
        
        rectfake_forward_expand = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [3.0,-0.5,0], [3.0,0.5,0],
            fill_color = RED, fill_opacity=0.5
        ).shift(INITIAL_SHIFT*UP)
        
        arrow_1 = Arrow(start=LEFT, end=RIGHT, color=RED).next_to(rectfake_forward_expand,RIGHT)
        Fs = Text("F",color=RED).next_to(arrow_1,RIGHT)
        Force = Group(arrow_1,Fs)
        
        interanl_force1 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=RED).next_to([2.0,-1.0,0],RIGHT)
        Fs = Text("f",color=RED).next_to(interanl_force1,RIGHT)
        IF1 = Group(interanl_force1,Fs)
        
        interanl_force2 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=RED).next_to([4.0,-1.0,0],RIGHT)
        Fs = Text("f + df",color=RED).next_to(interanl_force2,RIGHT)
        IF2 = Group(interanl_force2,Fs)
        
        
        
        # Describe Parsel
        Parsel_Group = self.Create_Parsel1()
        Parsel_Group2,xi2,dx,fadeoutgroup = self.Create_Parsel2()
        # Describe surface face for area
        surface_face = Rectangle(height=1, width=0.2, stroke_width=2,color=GREEN,fill_opacity=1.0).move_to([3.0,0,0]).shift(INITIAL_SHIFT*UP)
        surface_face2 = Rectangle(height=2, width=2, stroke_width=2,color=GREEN,fill_opacity=1.0).move_to([-4.0,-1.0,0])
        S_text = Text("S",color=GREEN).next_to(surface_face2,DOWN)
        
        if player == True:
            self.play(FadeIn(rectfake_forward))
            self.play(FadeIn(Force))
    
            self.wait(1)
            self.play(ReplacementTransform(rectfake_forward,rectfake_forward_expand),FadeIn(Parsel_Group))
            self.wait(1)
            self.play(FadeIn(surface_face))
            self.wait(2)
            self.play(ReplacementTransform(surface_face,surface_face2),FadeIn(S_text))
            self.wait(1)
            self.play(ReplacementTransform(Parsel_Group,Parsel_Group2))
            self.wait(1)
            self.play(FadeIn(IF1,IF2))
            
            Parsel_Total = Group(Parsel_Group2,IF1,IF2)
            Surface_Total = Group(S_text,surface_face2)
            
            self.play(FadeOut(rectfake,BC,rectfake_forward_expand,Force)) # Fadeout all
            # self.play(FadeOut(rectfake_forward_expand,Force)) # fadeout some
    
            
            Parsel_Total.generate_target()
            Parsel_Total.target.shift(3*UP)
            Surface_Total.generate_target()
            Surface_Total.target.shift(3*UP)
            
            self.play(MoveToTarget(Parsel_Total), MoveToTarget(Surface_Total))
            self.wait(2)
            
            Stress_eqn1 = MathTex(r"\sigma = \frac{f}{S}").move_to([-3.0,-2.0,0]).scale(2)
            Stress_eqn1[0][2].set_color(RED)
            Stress_eqn1[0][4].set_color(GREEN)
            Stress_words = Text("Stress").next_to(Stress_eqn1,UP).scale(1).set_color(RED_E)
            Stress_eqn2 = MathTex(r"").move_to([-3.0,-2.0,0])
    
            Strain_eqn1 = MathTex(r"\epsilon = \frac{d\xi}{dx}").move_to([3.0,-2.0,0]).scale(2)
            Strain_eqn1[0][2:4].set_color(BLUE)
            Strain_eqn1[0][5:7].set_color(GREY)
            Strain_words = Text("Strain").next_to(Strain_eqn1,UP*1.05).scale(1).set_color(BLUE_E)
            Strain_eqn2 = MathTex(r"").move_to([3.0,-2.0,0])
    
            Hookes = Text("= -(Modulus)").next_to(Stress_words,RIGHT).scale(1)
            Hookes_Law = MathTex(r"\frac{f}{S}", " =-E ", r"\frac{d\xi}{dx}").move_to([0.0,-2.0,0]).scale(2)
            Hookes_Law_Test = Text("Hooke's Law").next_to(Hookes,UP*1.1)
            Hookes_Law[0].set_color(RED_E)
            Hookes_Law[2].set_color(BLUE_E)
    
    
            Hookes_words_tot = Group(Strain_words,Stress_words,Hookes,Hookes_Law_Test)
            
            NecessaryEquations = Text("Necessary Equations").to_edge(UL)
            ul = Underline(NecessaryEquations)  # Underlining the word
            NecessaryEquations = VGroup(NecessaryEquations,ul)            
            Eqn1 = Text(" 1. Elastic Motion").scale(1).next_to(NecessaryEquations,DOWN).to_edge(LEFT)
            Eqn1_real = MathTex(r" 1. \frac{f}{S}", " =-E ", r"\frac{d\xi}{dx}").scale(1.4).next_to(NecessaryEquations,DOWN).to_edge(LEFT).set_color(BLUE)
    
            Eqn2 = Text("2. ???").scale(1).next_to(Eqn1,DOWN*6).to_edge(LEFT)
            box2 = SurroundingRectangle(Eqn2, corner_radius=0.2)
            Eqn2_Group = Group(Eqn2,box2)
            
            Eqn2_real = Text("2. Force Relation").scale(1).next_to(Eqn1,DOWN*6).to_edge(LEFT)
            box2_real = SurroundingRectangle(Eqn2_real, corner_radius=0.2)
            Eqn2_real_Group = Group(Eqn2_real,box2_real)
    
            Eqn3 = Text("3. ???").scale(1).next_to(Eqn2,DOWN*5).to_edge(LEFT)
            
            Eqns = Group(Eqn1,Eqn2,Eqn3)
            
            " Show strain and stress eqns "
            self.play(ReplacementTransform(S_text,Stress_eqn1),FadeIn(Stress_words),ReplacementTransform(IF1,Stress_eqn2)) # Stress Eqn
            self.wait(1)
            self.play(ReplacementTransform(xi2,Strain_eqn1),FadeIn(Strain_words),ReplacementTransform(dx,Strain_eqn2)) # Strain Eqn
            self.wait(2)
            self.play(FadeIn(Hookes),FadeIn(Hookes_Law_Test)) # FadeIn that this is hookes law
            self.wait(2)
            self.play(FadeOut(Strain_eqn1, Strain_eqn2, shift=LEFT),FadeOut(Stress_eqn1,Stress_eqn2,shift=RIGHT), FadeIn(Hookes_Law))
            self.wait(2)
            
            " Ending things to transition to next part"
            self.play(FadeOut(fadeoutgroup),FadeOut(surface_face2),FadeOut(IF2))
            self.play(ReplacementTransform(Hookes_words_tot,NecessaryEquations),ReplacementTransform(Hookes_Law,Eqn1),FadeIn(Eqn2,Eqn3))
            self.wait(1)
            self.play(Transform(Eqn1,Eqn1_real))
            self.play(Create(box2))
            self.wait(1)
            self.play(ReplacementTransform(Eqn2_Group,Eqn2_real_Group))
            
            
                    
            Bar = self.create_bar_up()

            
            self.wait(3)
    
            self.play(ReplacementTransform(Eqn2_real_Group,Bar),FadeOut(Eqn1_real,Eqn1,Eqn3,NecessaryEquations))
        else:
                                
            Bar = self.create_bar_up()

        return Bar

    def InternalForceScene3(self,Bar,player=True):
        arrow_1 = Arrow(start=LEFT, end=RIGHT, color=RED).next_to(Bar[0],RIGHT)
        Fs = MathTex(r"F",color=RED).next_to(arrow_1,RIGHT)
        Force = Group(arrow_1,Fs)
        
        
        Parsel1 = Rectangle(height=1, width=3, stroke_width=2,color=RED,fill_opacity=0.1).next_to(Bar[0],RIGHT).shift(LEFT*5)
        emphasis_line1 = Line([0,0,0],[0,2,0],color=RED).next_to(Parsel1,LEFT,buff=0).shift(0.5*UP)
        emphasis_line2 = Line([0,0,0],[0,2,0],color=RED).next_to(Parsel1,RIGHT,buff=0).shift(0.5*UP)
        IF1_arrow = Arrow(start=LEFT, end=RIGHT, color=RED).next_to(emphasis_line1[0],RIGHT,buff=0).shift(UP/2)
        IF1 = MathTex(r"f",color=RED).next_to(IF1_arrow,buff = 0.1)
        IF2_arrow = Arrow(start=LEFT, end=RIGHT, color=RED).next_to(emphasis_line2[0],RIGHT,buff=0).shift(UP/2)
        IF2 = MathTex(r"f + df",color=RED).next_to(IF2_arrow,buff = 0.1)
        
        IF_Group = VGroup(emphasis_line1,emphasis_line2,IF1_arrow,IF1,IF2_arrow,IF2)
        
        x_val   = MathTex("x",color=GRAY).next_to(emphasis_line1,DOWN*1.3)
        dx_val  = MathTex("x + dx",color=GRAY).next_to(emphasis_line2,DOWN)
        brace = BraceBetweenPoints(x_val.get_center() , dx_val.get_center(),color=GREY)
        dx_val_real = MathTex("dx",color=GRAY).next_to(brace,DOWN)
        
        dx_Group = VGroup(x_val,dx_val,brace,dx_val_real)
        
        
        
        Parsel2 = Rectangle(height=1, width=3, stroke_width=2,color=BLUE,fill_opacity=0.1).next_to(Bar[0],RIGHT).shift(LEFT*4)


        
        IF_Text = Text("What are the total forces on this bar to the right?").scale(0.7).next_to(dx_Group[3],DOWN)

        boxT1 = SurroundingRectangle(IF_Group[3], corner_radius=0.2,color=ORANGE)
        T1_Group = VGroup(boxT1,IF_Group[3])
        
        boxT2_1 = SurroundingRectangle(IF_Group[-1], corner_radius=0.2)
        T21_Group = VGroup(boxT2_1,IF_Group[-1])
        boxT2_2 = SurroundingRectangle(dx_Group[1], corner_radius=0.2)
        T22_Group = VGroup(boxT2_2,dx_Group[1])       
        
        Fnet = MathTex(r"F_{net}"," ="," f ","-", r"(f + " , r"\frac{df}{dx}dx",")").scale(1.5).next_to(IF_Text,DOWN)
        
        Fnet[0].set_color(RED)
        Fnet[2].set_color(ORANGE)
        Fnet[4].set_color(YELLOW)
        Fnet[5].set_color(YELLOW)
        Fnet[6].set_color(YELLOW)

        if player == True:
            self.play(FadeIn(Force))

            self.play(FadeIn(Parsel1),FadeIn(Parsel2),FadeIn(IF_Group))
            self.play(FadeIn(dx_Group))
            self.play(Write(IF_Text))
            self.wait(2)
            self.play(Create(boxT1))
            self.wait(1)
            self.play(ReplacementTransform(T1_Group,Fnet[2]))
            self.wait(1)
            self.play(Create(boxT2_1),Create(boxT2_2))
            self.wait(1)
            self.play(ReplacementTransform(T21_Group,Fnet[4]),ReplacementTransform(T22_Group,Fnet[5:7]))
            self.wait(1)
            self.play(FadeIn(Fnet[0]),FadeIn(Fnet[3]),FadeIn(Fnet[1]))
            self.play(FadeOut(Parsel1),FadeOut(Parsel2),FadeOut(dx_Group),FadeOut(IF_Group),FadeOut(Force))
    
    
            Fnet[5].generate_target()
            Fnet[5].target.shift(LEFT*1.3)
            
            Fnet[0].generate_target()
            Fnet[0].target.shift(RIGHT*1.5)
            Fnet[1].generate_target()
            Fnet[1].target.shift(RIGHT*1.5) 
            Negative_Sign = MathTex("-").next_to(Fnet[1],RIGHT,buff = 0.1).shift(RIGHT*1.5).set_color(YELLOW)

            Fnet2 = VGroup(Fnet[0],Fnet[1],Fnet[5],Negative_Sign)
            self.wait(1)
            self.play(FadeOut(Fnet[2],Fnet[3:5],Fnet[-1]))
            self.play(MoveToTarget(Fnet[5]),FadeIn(Negative_Sign),MoveToTarget(Fnet[0]),MoveToTarget(Fnet[1]))
            self.wait(2)
            NecessaryEquations = Text("Necessary Equations").to_edge(UL)
            ul = Underline(NecessaryEquations)  # Underlining the word
            NecessaryEquations = VGroup(NecessaryEquations,ul)            
            Eqn1 = Text(" 1. Elastic Motion").scale(1).next_to(NecessaryEquations,DOWN).to_edge(LEFT)
            Eqn1_real = MathTex(r" 1. \frac{f}{S}", " =-E ", r"\frac{d\xi}{dx}").scale(1.4).next_to(NecessaryEquations,DOWN).to_edge(LEFT).set_color(BLUE)
    
            Eqn2 = Text("2. Forcing Relation").scale(1).next_to(Eqn1,DOWN*8).to_edge(LEFT)
            Eqn2_real = MathTex(r"2. F_{net} = -\frac{df}{dx}dx").scale(1.4).next_to(Eqn1,DOWN*6).to_edge(LEFT).set_color(RED)
    
            Eqn3 = Text("3. ???").scale(1).next_to(Eqn2_real,DOWN*4).to_edge(LEFT)
            box3 = SurroundingRectangle(Eqn3, corner_radius=0.2)
            Eqn3_Group = Group(Eqn3,box3)
                
            Eqn3_real = Text("3. Newton Second Law").scale(1).next_to(Eqn2_real,DOWN*4).to_edge(LEFT)
            box3_real = SurroundingRectangle(Eqn3_real, corner_radius=0.2)
            Eqn3_real_Group = Group(Eqn3_real,box3_real)
            
            Eqns = Group(Eqn1,Eqn2,Eqn3)
            self.play(FadeOut(IF_Text),ReplacementTransform(Fnet2,Eqn2),Transform(Bar,NecessaryEquations),FadeIn(Eqn1_real,Eqn3))
            self.wait(2)
            self.play(ReplacementTransform(Eqn2,Eqn2_real))
            self.play(Create(box3))
            self.wait(1)
            self.play(ReplacementTransform(Eqn3_Group,Eqn3_real_Group))
            
            
            Barnew = self.create_bar_up()
     
            self.wait(2)
            self.play(ReplacementTransform(Eqn3_real_Group,Barnew),FadeOut(Eqn1_real,Eqn2_real,Bar))
        else:

            Barnew = self.create_bar_up()
            
        return Barnew

    def Newton2ndLawScene4(self,Barnew,player=True):
        #self.add(Barnew)
        " Newton Face "
        img1 = ImageMobject('newton.png')
        img1.to_corner(LEFT+DOWN)
        img2 = ImageMobject('google_eyes.png').scale(0.85).next_to(img1,RIGHT).shift(LEFT*2.27).shift(UP)

        " Parsel Group"
        Parsel2 = Rectangle(height=1, width=3, stroke_width=2,color=WHITE,fill_opacity=0.1).next_to(Barnew[0],RIGHT).shift(LEFT*5)
        emphasis_line1 = Line([0,0,0],[0,2,0],color=BLUE).next_to(Parsel2,LEFT,buff=0).shift(0.5*UP)
        emphasis_line2 = Line([0,0,0],[0,2,0],color=BLUE).next_to(Parsel2,RIGHT,buff=0).shift(0.5*UP)
        
        IF1_arrow = Arrow(start=LEFT, end=RIGHT, color=BLUE).next_to(emphasis_line1[0],RIGHT,buff=0).shift(UP/2)
        IF1 = MathTex(r"\xi",color=BLUE).next_to(IF1_arrow,buff = 0.1)
        IF2_arrow = Arrow(start=LEFT, end=RIGHT, color=BLUE).next_to(emphasis_line2[0],RIGHT,buff=0).shift(UP/2)
        IF2 = MathTex(r"\xi + d\xi",color=BLUE).next_to(IF2_arrow,buff = 0.1)
        
        brace = BraceBetweenPoints(emphasis_line1.get_center() , emphasis_line2.get_center(),color=GREY).shift(DOWN)
        dx_val_real = MathTex("dx",color=GRAY).next_to(brace,DOWN)
        brace_group = VGroup(brace,dx_val_real)
        brace_box = SurroundingRectangle(brace_group, corner_radius=0.2)
        brace_box_group = VGroup(brace_group,brace_box)
        
        
        deflection_group = VGroup(emphasis_line1,emphasis_line2,IF1_arrow,IF1,IF2_arrow,IF2)
        deflection_group2 = VGroup(Parsel2,deflection_group)
        newtonlaw = Text("Newton's 2nd Law").next_to(Barnew,DOWN).shift(RIGHT*3)
        newtonlaw = VGroup(newtonlaw,Underline(newtonlaw))
        newtonlaw_math = MathTex(r"F_{net} =  ",r"m", r"a").scale(1.2).next_to(newtonlaw,DOWN)
        
        " Surface Area Discussion Boxing for mass derivation "
        surface_face = Rectangle(height=1, width=0.2, stroke_width=2,color=GREEN,fill_opacity=1.0).next_to(Barnew[0],RIGHT,buff = 0).shift(LEFT/8)
        area = MathTex("S",color=GREEN).next_to(surface_face,RIGHT,buff=0.1)
        area_group = VGroup(surface_face,area)
        boxA = SurroundingRectangle(area_group, corner_radius=0.2)
        Area_Group = VGroup(brace_box_group,area_group,boxA)
        box_parsel = SurroundingRectangle(deflection_group2, corner_radius=0.2)
        parsel_box_group = VGroup(deflection_group2,box_parsel)
        
        " Final math Transform Section "
        mass = MathTex(r'm = \rho S dx',color=GREEN).scale(1.2).next_to(newtonlaw_math,DOWN)
        accel = MathTex(r'a = \frac{d^2\xi}{dt^2}',color=BLUE).scale(1.2).next_to(mass,DOWN)
        newtonlaw_full = MathTex(r"F_{net} =  ",r"\rho S dx", r"\frac{d^2\xi}{dt^2}").scale(1.2).next_to(newtonlaw,DOWN)
        newtonlaw_full[1].set_color(GREEN)
        newtonlaw_full[2].set_color(BLUE)

        eqn_group = VGroup(mass,accel)
        
        if player == True:
            " Animation Section "
            self.play(FadeIn(img1)) # Get image of newton
            self.play(FadeIn(newtonlaw),FadeIn(newtonlaw_math)) # Display Newton meth
            self.wait(1)
            
            " Introduce variables for equatins"
            self.play(FadeIn(Parsel2))
            self.play(FadeIn(deflection_group),FadeIn(brace_group))
            self.play(FadeIn(surface_face),FadeIn(area))
            self.wait(2)
            
            " Mass Derivation"
            self.play(Create(boxA),Create(brace_box))
            self.wait(2)
            newtonlaw_math[1].set_color(GREEN)
            self.play(ReplacementTransform(Area_Group,mass),TransformMatchingTex(newtonlaw_math,newtonlaw_math))        
            self.wait(2)
            
            " Acceleration Derivation"
            self.play(Create(box_parsel))
            self.wait(2)
            newtonlaw_math[2].set_color(BLUE)
            self.play(ReplacementTransform(parsel_box_group,accel),TransformMatchingTex(newtonlaw_math,newtonlaw_math))      
            self.wait(2)
            Transform_Group = VGroup(eqn_group,newtonlaw_math)
            " Final Equations "
            self.play(ReplacementTransform(Transform_Group,newtonlaw_full))
            self.wait(1)
            self.play(FadeIn(img2)) # Newton say what?!?
            self.wait(2)
                    
            " Final Equations for necessary equation review "
            NecessaryEquations = Text("Necessary Equations").to_edge(UL)
            ul = Underline(NecessaryEquations)  # Underlining the word
            NecessaryEquations = VGroup(NecessaryEquations,ul)
            Eqn1_real = MathTex(r" 1. \frac{f}{S}", " =-E ", r"\frac{d\xi}{dx}").scale(1.4).next_to(NecessaryEquations,DOWN).to_edge(LEFT).set_color(BLUE)
            Eqn1 = Text(" 1. Elastic Motion").scale(1).next_to(NecessaryEquations,DOWN).to_edge(LEFT)
    
            Eqn2_real = MathTex(r"2. F_{net} = -\frac{df}{dx}dx").scale(1.4).next_to(Eqn1,DOWN*6).to_edge(LEFT).set_color(RED)
                
            Eqn3 = Text("3. Newton Second Law").scale(1).next_to(Eqn2_real,DOWN*4).to_edge(LEFT)
            Eqn3_real = MathTex(r"3. F_{net} = \rho S dx \frac{d^2\xi}{dt^2}").scale(1.4).next_to(Eqn2_real,DOWN*2.5).to_edge(LEFT).set_color(GREEN)
    
            Eqn_Group = VGroup(NecessaryEquations,Eqn1_real,Eqn2_real,Eqn3_real,newtonlaw_full)
            self.play(FadeOut(newtonlaw),FadeOut(img1),FadeOut(img2),ReplacementTransform(newtonlaw_full,Eqn3),ReplacementTransform(Barnew,Eqn_Group[0]),FadeIn(Eqn_Group[1],Eqn_Group[2]))
            self.wait(2)
            self.play(ReplacementTransform(Eqn3,Eqn_Group[3]))
            
            return Eqn_Group
        
        else:             
            " Final Equations for necessary equation review "
            NecessaryEquations = Text("Necessary Equations").to_edge(UL)
            ul = Underline(NecessaryEquations)  # Underlining the word
            NecessaryEquations = VGroup(NecessaryEquations,ul)
            Eqn1_real = MathTex(r" 1. \frac{f}{S}", " =-E ", r"\frac{d\xi}{dx}").scale(1.4).next_to(NecessaryEquations,DOWN).to_edge(LEFT).set_color(BLUE)
            Eqn1 = Text(" 1. Elastic Motion").scale(1).next_to(NecessaryEquations,DOWN).to_edge(LEFT)
    
            Eqn2_real = MathTex(r"2. F_{net} = -\frac{df}{dx}dx").scale(1.4).next_to(Eqn1,DOWN*6).to_edge(LEFT).set_color(RED)
                
            Eqn3 = Text("3. Newton Second Law").scale(1).next_to(Eqn2_real,DOWN*4).to_edge(LEFT)
            Eqn3_real = MathTex(r"3. F_{net} = \rho S dx \frac{d^2\xi}{dt^2}").scale(1.4).next_to(Eqn2_real,DOWN*2.5).to_edge(LEFT).set_color(GREEN)
    
            Eqn_Group = VGroup(NecessaryEquations,Eqn1_real,Eqn2_real,Eqn3_real)

            return Eqn_Group

    def WaveEquationScene5(self,Eqns,player=True):
                
        " Emphasize old Eqnuations"
        box_eqn1 = SurroundingRectangle(Eqns[1], corner_radius=0.2)
        box_eqn2 = SurroundingRectangle(Eqns[2], corner_radius=0.2)
        box_eqn3 = SurroundingRectangle(Eqns[3], corner_radius=0.2)
        
        Eqn1_group = VGroup(Eqns[1],box_eqn1)
        Eqn2_group = VGroup(Eqns[2],box_eqn2)
        Eqn3_group = VGroup(Eqns[3],box_eqn3)
        " New Equations "
        Eqn1_new = MathTex(r"1. f", " =-SE ", r"\frac{d\xi}{dx}").scale(1.4).next_to(Eqns[0],DOWN).to_edge(LEFT).set_color(BLUE)
        box_eqn_new = SurroundingRectangle(Eqn1_new, corner_radius=0.2)
        Eqn1_group_new = VGroup(Eqn1_new,box_eqn_new)

        Eqn_1_2  = MathTex(r"-\frac{df}{dx}","dx","=",r"\rho S","dx",r"\frac{d^2\xi}{dt^2}").scale(1.4).to_edge(RIGHT).shift(LEFT)
        Eqn_1_2[0:2].set_color(RED)
        Eqn_1_2[3:6].set_color(GREEN)
        
        if player == True:
            self.wait(2)
            self.play(Create(box_eqn2),Create(box_eqn3))
            self.wait(2)
            self.play(ReplacementTransform(Eqn2_group,Eqn_1_2[0:2]))
            self.wait(1)
            self.play(ReplacementTransform(Eqn3_group,Eqn_1_2[2:6]))
            self.wait(1)
            
            " Cancel Terms dx eqn1_2"
            cancel_t1 = cancel_line(Eqn_1_2[1],ORANGE)
            cancel_t2 = cancel_line(Eqn_1_2[4],ORANGE)
            Eqn_1_2[0].generate_target()
            Eqn_1_2[0].target.shift(RIGHT*0.8)
    
            
            Eqn_1_2[5].generate_target()
            Eqn_1_2[5].target.shift(LEFT*0.8)
    
            self.play(FadeIn(cancel_t1),FadeIn(cancel_t2))
            self.wait(1)
            self.play(FadeOut(cancel_t1),FadeOut(Eqn_1_2[1]),FadeOut(cancel_t2),
                      FadeOut(Eqn_1_2[4]),MoveToTarget(Eqn_1_2[0]),MoveToTarget(Eqn_1_2[5]))
            self.wait(2)
            
            " Bring in Equation 1 "
            self.play(Create(box_eqn1))
            self.wait(1)
            self.play(ReplacementTransform(Eqn1_group,Eqn1_group_new))
            self.wait(1)
            box_f_term = SurroundingRectangle(Eqn_1_2[0][2], corner_radius=0.2)
            self.play(Create(box_f_term))
    
            " Transform to singular Equation "
            Transform_To_Final_Eqn_Group = VGroup(Eqn1_group_new,Eqn_1_2[0],Eqn_1_2[2],Eqn_1_2[3],Eqn_1_2[5],box_f_term)
            Eqn_1_2_final = MathTex(r"\frac{1}{dx}","S",r"E \frac{d\xi}{dx}","=",r"\rho", "S",r"\frac{d^2\xi}{dt^2}").scale(1.4).shift(UP*1.5)
            Eqn_1_2_final[0].set_color(RED)
            Eqn_1_2_final[1:3].set_color(BLUE)
            Eqn_1_2_final[4:7].set_color(GREEN)
            
            # All this crap gotta go, just here for canceling out and emphasis
            cancel_t1 = cancel_line(Eqn_1_2_final[1],ORANGE)
            cancel_t2 = cancel_line(Eqn_1_2_final[5],ORANGE)
            box_t1 = SurroundingRectangle(Eqn_1_2_final[0], corner_radius=0.2)
            box_t2 = SurroundingRectangle(Eqn_1_2_final[2], corner_radius=0.2)
    
            
            self.wait(1)
            self.play(FadeOut(Eqns[0]),ReplacementTransform(Transform_To_Final_Eqn_Group,Eqn_1_2_final))
            self.wait(2)
            self.play(FadeIn(cancel_t1),FadeIn(cancel_t2))
            self.wait(1)
            self.play(Create(box_t1),Create(box_t2))
            
            " Final Wave Equation Derivations "
            Question = Text("Can we create an equation to describe this motion? ").scale(0.8).next_to(Eqn_1_2_final,UP)
            Answer = Text("The Wave Equation ").scale(0.8).next_to(Question,DOWN)
            Answer = VGroup(Answer,Underline(Answer))
    
            Wave_eqn = MathTex(r"E \frac{d^2\xi}{dx^2}","=",r"\rho\frac{d^2\xi}{dt^2}").scale(1.4).next_to(Eqn_1_2_final,DOWN*2)
            math_arrow = Arrow(start=Eqn_1_2_final.get_center(), end=Wave_eqn.get_center(), color=WHITE).shift(DOWN*0.8)
            Wave_eqn.shift(DOWN)
            Final_Transform_Group = VGroup(Wave_eqn,math_arrow,Eqn_1_2_final,cancel_t1,cancel_t2,box_t1,box_t2)
            Wave_eqn2 = MathTex(r"E \frac{\partial ^2\xi}{\partial x^2}","=",r"\rho\frac{\partial ^2\xi}{\partial t^2}").scale(1.5)
            
            self.play(FadeIn(math_arrow))
            self.wait(1)
            self.play(FadeIn(Wave_eqn))
            self.wait(1)
            self.play(Write(Question))
            self.wait(3)
            self.play(ReplacementTransform(Final_Transform_Group,Wave_eqn2),Write(Answer))
            
        else:
            None
            
    def construct(self):
        Bar1,BC,rectfake = self.DemonstrateLongWaveMotionScene1(True)
        Bar              = self.HookesLawScene2(Bar1,BC,rectfake,True)
        Barnew           = self.InternalForceScene3(Bar,True)
        Eqns             = self.Newton2ndLawScene4(Barnew,True)
        self.WaveEquationScene5(Eqns,True)
        self.wait(5)

        

class Galerkin_KeyFrame(Scene):
    def Scene_1_Demonstrate_Long_Wave_Movement(self):
        # Mobject Creation Beam Bar Physical Creation
        width = ValueTracker(3.0)

        rect1 = always_redraw(lambda : Polygon(
            [-3,0.5,0], [-3,-0.5,0], [width.get_value(),-0.5,0], [width.get_value(),0.5,0],
            fill_color = BLUE, fill_opacity=0.5
        )) # Longitudinal movement

        
        rectfake = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [2.0,-0.5,0], [2.0,0.5,0],
            fill_color = BLUE, fill_opacity=0.5
        ) # Stationary bar movement
   
        rectfake_forward = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [2.0,-0.5,0], [2.0,0.5,0],
            fill_color = BLUE, fill_opacity=0.1
        ).shift(INITIAL_SHIFT*UP) # Stationary bar movement
        
        rectfake_forward_expand = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [3.0,-0.5,0], [3.0,0.5,0],
            fill_color = RED, fill_opacity=0.5
        ).shift(INITIAL_SHIFT*UP)
        # Boundary Condition fixed
        l  = Line([-3., 1., 0.], [-3., -1., 0.])
        l2 = Line([-3., 1., 0.], [-4., 0.5, 0.])
        l3 = Line([-3., 0.5, 0.], [-4., 0.0, 0.])
        l4 = Line([-3., 0., 0.], [-4., -0.5, 0.])
        l5 = Line([-3., -0.5, 0.], [-4., -1.0, 0.])
        l6 = Line([-3., -1., 0.], [-4., -1.5, 0.])
        BC = Group(l,l2,l3,l4,l5,l6)

        # Forcing function
        arrow_1 = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(rect1,RIGHT)
        Fs = Text("F",color=RED).next_to(arrow_1,RIGHT)
        Force = Group(arrow_1,Fs)
        Bar = Group(BC,rect1)


        # Differential Equations
        MainDiff = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = F")
        MainDiff[3].set_color(RED)
        for ii in range(3):
            MainDiff[ii].set_color(BLUE)

        Bar1 = Group(BC,rectfake)

    
    
        return Bar1,BC,rectfake
    
    def Create_Parsel1(self):
        # Describe Parsel
        Parsel1 = Rectangle(height=1, width=1, stroke_width=2,color=WHITE,fill_opacity=1.0).move_to([2.5,0,0]).shift(INITIAL_SHIFT*UP)
        emphasisline1 = Line([2.0,-0.5,0],[2.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        emphasisline2 = Line([3.0,-0.5,0],[3.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        arrow_wave1 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([1.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        arrow_wave2 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([2.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        xi1 = MathTex(r"\xi",color=BLUE).next_to(arrow_wave1,UP/2)
        xi2 = MathTex(r"\xi + d\xi",color=BLUE).next_to(arrow_wave2,UR/2,aligned_edge=LEFT)
        brace = BraceBetweenPoints([2.0,-0.5,0],[3.0,-0.5,0],color=GREY).shift(INITIAL_SHIFT*UP)
        dx = Text("dx").next_to(brace,DOWN).set_color(GREY)
        Parsel_Group = Group(Parsel1,brace,dx,emphasisline1,emphasisline2,arrow_wave1,arrow_wave2,xi1,xi2)

        return Parsel_Group
    
    def Create_Parsel2(self):
        # Describe Parsel
        Parsel1 = Rectangle(height=1, width=2, stroke_width=2,color=WHITE,fill_opacity=1.0).move_to([3.0,0,0]).shift(INITIAL_SHIFT*UP)
        emphasisline1 = Line([2.0,-0.5,0],[2.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        emphasisline2 = Line([4.0,-0.5,0],[4.0,1.2,0],color=BLUE).shift(INITIAL_SHIFT*UP)
        arrow_wave1 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([1.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        arrow_wave2 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=BLUE).next_to([3.8,1.2,0],RIGHT).shift((INITIAL_SHIFT-0.2)*UP)
        xi1 = MathTex(r"\xi",color=BLUE).next_to(arrow_wave1,UP/2)
        xi2 = MathTex(r"\xi + d\xi",color=BLUE).next_to(arrow_wave2,UR/2,aligned_edge=LEFT)
        brace = BraceBetweenPoints([2.0,-0.5,0],[4.0,-0.5,0],color=GREY).shift(INITIAL_SHIFT*UP)
        dx = Text("dx").next_to(brace,DOWN).set_color(GREY)
        
        bracetot = Group(brace,dx)
        displacetot = Group(xi2, arrow_wave2,emphasisline2)
        
        Parsel_Group = Group(Parsel1,brace,dx,emphasisline1,emphasisline2,arrow_wave1,arrow_wave2,xi1,xi2).move_to([4.0,-1.0,0])

        return Parsel_Group,displacetot,bracetot
    
    def construct(self):
        # Mobject Creation Beam Bar Physical Creation
        width = ValueTracker(3.0)
              
        # Set up problem
        rectfake_forward = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [2.0,-0.5,0], [2.0,0.5,0],
            fill_color = BLUE, fill_opacity=0.1
        ).shift(INITIAL_SHIFT*UP) # Stationary bar movement
        
        rectfake_forward_expand = Polygon(
            [-3,0.5,0], [-3,-0.5,0], [3.0,-0.5,0], [3.0,0.5,0],
            fill_color = RED, fill_opacity=0.5
        ).shift(INITIAL_SHIFT*UP)
        
        arrow_1 = Arrow(start=LEFT, end=RIGHT, color=RED).next_to(rectfake_forward_expand,RIGHT)
        Fs = Text("F",color=RED).next_to(arrow_1,RIGHT)
        Force = Group(arrow_1,Fs)
        
        interanl_force1 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=RED).next_to([2.0,-1.0,0],RIGHT)
        Fs = Text("f",color=RED).next_to(interanl_force1,RIGHT)
        IF1 = Group(interanl_force1,Fs)
        
        interanl_force2 = Arrow(max_tip_length_to_length_ratio=0.2,start=LEFT, end=RIGHT/4, color=RED).next_to([4.0,-1.0,0],RIGHT)
        Fs = Text("f + df",color=RED).next_to(interanl_force2,RIGHT)
        IF2 = Group(interanl_force2,Fs)
        
        
        
        # Describe Parsel
        Parsel_Group = self.Create_Parsel1()
        Parsel_Group2,xi2,dx = self.Create_Parsel2()
        # Describe surface face for area
        surface_face = Rectangle(height=1, width=0.2, stroke_width=2,color=GREEN,fill_opacity=1.0).move_to([3.0,0,0]).shift(INITIAL_SHIFT*UP)
        surface_face2 = Rectangle(height=2, width=2, stroke_width=2,color=GREEN,fill_opacity=1.0).move_to([-4.0,-1.0,0])
        S_text = Text("S",color=GREEN).next_to(surface_face2,DOWN)
        
        Bar1,BC,rectfake = self.Scene_1_Demonstrate_Long_Wave_Movement()
       # Bar1.shift(2*UP)
        
        Parsel_Total = Group(Parsel_Group2,IF1,IF2)
        Surface_Total = Group(S_text,surface_face2)

        
        Parsel_Total.generate_target()
        Parsel_Total.target.shift(3*UP)
        Surface_Total.generate_target()
        Surface_Total.target.shift(3*UP)
        self.add(Bar1)#,Parsel_Total,Surface_Total)





        
        
        
        
        

