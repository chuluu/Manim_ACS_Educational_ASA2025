# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:20:13 2023

@author: mbluu
"""

from manim import *
import numpy as np

CONFIG = {"include_numbers": True,
          "include_tip": False}

CONFIG2 = {"include_numbers": False,
          "include_tip": False}

def function(x):
    return x**2

def inv_func(x):
    return 2*x

def deriv(x):
    return 2*x


def emphasis_box(term,color=YELLOW):
    box = SurroundingRectangle(term, corner_radius=0.2,color=color)
    emphasized_term = VGroup(term,box)
    
    return emphasized_term

def arrow_text(text,arrow_coord,color=YELLOW):
    arrow = Arrow(start=arrow_coord[0], end=arrow_coord[1], color=color)
    txt = Text(text, color=color).scale(0.5).next_to(arrow,LEFT,buff = 0.1)
    
    return VGroup(arrow,txt)

class FE_Difference_Inverse(Scene):
    def MathScene1(self,player = True):
        txt = Text("Using Central Differences").to_edge(UP).set_color(GREEN)
        divider = Line([-4,0,0], [4,0,0]).next_to(txt,DOWN)

        TextExplain = Paragraph('Solve differential equation', 'using central differences',line_spacing=1.2).scale(0.9).next_to(divider,DOWN)
        TextIDA = Paragraph("IDEA: ").next_to(TextExplain[0],LEFT).set_color(RED).scale(0.9)
        Group1Explain = VGroup(TextExplain,TextIDA)

        
        
        Finite_Difference = MathTex(r"\frac{dy}{dx} =", r"\frac{y_{i+1} - y_{i-1}}{2dx}").next_to(Group1Explain,DOWN*1.5).shift(RIGHT)

        diffeqn = MathTex(r"\frac{dy}{dx}","=","2").next_to(Finite_Difference,DOWN)

        fe_difference_eqn   = MathTex(r" {y_{i+1}"," -",r" y_{i-1} \over ","2dx}","=","2")
        fe_difference_eqn_2 = MathTex(r"y_{i+1}","-", r"y_{i-1}","=","4dx")
        fe_difference_eqn_3 = MathTex(r"y_{i+1}", r"=" ,"4dx","+",r" y_{i-1}")
        fe_difference_eqn_3 = MathTex(r"y_{i+1}", r"=" ,"4dx","+",r" y_{i-1}")

        Explain_text_PD = Text("Point difference ").scale(0.6).next_to(fe_difference_eqn_3[2],UP).set_color(ORANGE)
        Explain_text_IC = Text("Boundary conditions").scale(0.6).next_to(fe_difference_eqn_3[-1],DOWN).set_color(YELLOW)
        Explain_text_eqn = arrow_text("Original Function",[LEFT,RIGHT],color=BLUE)
        Explain_text_eqn.next_to(fe_difference_eqn_3[0],LEFT)
        
        diffeqn_boxed =emphasis_box(fe_difference_eqn[-1])
        diffeqn_boxed_PD =emphasis_box(fe_difference_eqn_3[2],ORANGE)
        diffeqn_boxed_IC =emphasis_box(fe_difference_eqn_3[-1],YELLOW)

        func_tot = MathTex(r"y_{i+1}", r"=", r"\frac{dy}{dx}","2dx", "+", r" y_{i-1}").scale(0.8).to_edge(UR).shift(RIGHT/6)
        func_tot[0].set_color(RED)
        func_tot[3].set_color(GRAY)
        func_tot[-1].set_color(MAROON)

        axes = (
            Axes(
                x_range = [0,10,2],
                x_length = 8,
                y_range = [0,20,4],
                y_length = 2.4,
                axis_config=CONFIG,
                ).to_edge(DL).set_color(GREY)
            )
        grid_labels_OG = axes.get_axis_labels(
            Tex("x").scale(0.7), Tex("y").scale(0.7)
        ).shift(DOWN/4).set_color(LIGHT_GREY)
        OG_funcbox = SurroundingRectangle(axes, corner_radius=0.2,buff = 0.4)

        axes_derivative = (
            Axes(
                x_range = [0,10,2],
                x_length = 8,
                y_range = [0,10,2],
                y_length = 2.4,
                axis_config=CONFIG,
                ).to_edge(UL).set_color(GREY)
            )
        grid_labels_deriv = axes_derivative.get_axis_labels(
            Tex("x").scale(0.7), MathTex(r"\frac{dy}{dx}").scale(0.7).set_color(LIGHT_GREY)
        ).shift(DOWN/4)
        
        Derivative_funcbox = SurroundingRectangle(axes_derivative, corner_radius=0.2,buff=0.4)
        if player == True:
            self.play(FadeIn(txt))
            self.play(Write(divider))
            self.play(FadeIn(Group1Explain))
            self.wait(1)
            self.play(Write(Finite_Difference))
            self.wait(2)
            self.play(Write(diffeqn))
            self.wait(3)
            self.play(TransformMatchingTex(diffeqn,fe_difference_eqn),FadeOut(Finite_Difference, target_position=diffeqn))
            self.wait(1)
            self.play(FadeOut(Group1Explain))
            self.wait(1)
            self.play(Create(diffeqn_boxed[1]))
            self.wait(1)
            self.play(FadeOut(diffeqn_boxed[1]), TransformMatchingTex(fe_difference_eqn,fe_difference_eqn_2))
            self.wait(1)
            self.play(TransformMatchingTex(fe_difference_eqn_2,fe_difference_eqn_3))
            self.wait(2)
            self.play(Create(diffeqn_boxed_PD[1]))
            self.wait(1)
            self.play(FadeIn(Explain_text_PD))
            self.wait(2)
            self.play(Create(diffeqn_boxed_IC[1]))
            self.wait(1)
            self.play(FadeIn(Explain_text_IC))
            self.wait(2)
            self.play(FadeIn(Explain_text_eqn))
            self.wait(4)
            self.play(FadeOut(divider, txt))
            self.play(ReplacementTransform(fe_difference_eqn_3,func_tot),
                      FadeOut(diffeqn_boxed_PD[1]),
                      FadeOut(diffeqn_boxed_IC[1]),
                      FadeOut(Explain_text_PD),
                      FadeOut(Explain_text_IC),
                      FadeOut(Explain_text_eqn),
                      FadeIn(axes),
                      FadeIn(grid_labels_OG),
                      FadeIn(axes_derivative),
                      FadeIn(grid_labels_deriv))
            self.play(Create(OG_funcbox),
                      Create(Derivative_funcbox))
            self.wait(2)
        else:
            None
            
        return [func_tot,axes,grid_labels_OG,OG_funcbox,axes_derivative,grid_labels_deriv,Derivative_funcbox]
    
    def construct(self):
        player = True
        T1 = MAROON
        T2 = RED
        DX_COL = BLUE
        
        [func_tot,axes,grid_labels_OG,OG_funcbox,axes_derivative,grid_labels_deriv,Derivative_funcbox] = self.MathScene1(True)
        self.add(func_tot,axes,grid_labels_OG,OG_funcbox,axes_derivative,grid_labels_deriv,Derivative_funcbox)
        diffeqn = MathTex(r"\frac{dy}{dx}","=","2").scale(0.8).next_to(func_tot,DOWN)

        x = [0, 1]
        y = [0, 2]
        line = axes.plot_line_graph(x, y, add_vertex_dots=True)
        x = [0, 1,2,3,4,5,6,7,8,9,10]
        y = [2, 2,2,2,2,2,2,2,2,2,2]
        derivline = axes_derivative.plot_line_graph(x, y, add_vertex_dots=False)


        " Start Creating yo "
        x = ValueTracker(0)
        e = ValueTracker(1)
        dx = ValueTracker(2)
        func = axes.plot(lambda x: 2*x, x_range = [0,10], color = YELLOW)
        func_deriv = axes.plot(lambda x: 2, x_range = [0,10], color = YELLOW)

        Dot_x = always_redraw(
            lambda: Dot(color = T1).scale(1).move_to(axes.c2p(x.get_value(), 
                                                    func.underlying_function(x.get_value()))))
        
        Dot_dx = always_redraw(
            lambda: Dot(color = T2).scale(1).move_to(axes.c2p(x.get_value() + dx.get_value(), 
                                                    func.underlying_function(x.get_value() + dx.get_value()))))
        

        Dot_deriv = always_redraw(
            lambda: Dot(color = WHITE).scale(1.2).move_to(axes_derivative.c2p(e.get_value(), 
                                                    func_deriv.underlying_function(e.get_value()))))


        func_OG = always_redraw(lambda: 
                                axes.plot(lambda x: (2*x),
                                                     x_range = [1,x.get_value()+2], color = YELLOW))

        brace_dx = always_redraw( 
            lambda: BraceBetweenPoints(Dot_x.get_center() , Dot_dx.get_center(),color=GREY,buff=0.1 ))
        
        dx_val_real =always_redraw( 
            lambda:  MathTex("2dx",color=GRAY).scale(0.5).next_to(brace_dx,DOWN,buff=0).shift(RIGHT/8).rotate(PI/10))
        
        brace_dx_group = VGroup(brace_dx,dx_val_real)
        #Dot_OG = always_redraw(
        #    lambda: Dot(color = WHITE).scale(1.2).move_to(axes.c2p(e.get_value(), 
        #                                            func.underlying_function(e.get_value()))))
        
        eqn1 = always_redraw(lambda: MathTex(r"y_{i-1}",color=T1).scale(0.6).next_to(Dot_x,UP,buff=0.15))
        eqn2 = always_redraw(lambda:  MathTex(r"y_{i+1}",color=T2).scale(0.6).next_to(Dot_dx,UP,buff=0.15))
        
        
        tt1 = 0 
        y1_val = Variable(tt1, 'x', num_decimal_places=3)
        on_screen_var_y1 = Variable(inv_func(tt1), MathTex("y_{i-1}"), num_decimal_places=2).next_to(diffeqn,DOWN*2)
        on_screen_var_y1.label.set_color(WHITE)
        on_screen_var_y1.value.set_color(T1)
        on_screen_var_y1.add_updater(lambda a: a.tracker.set_value(inv_func(y1_val.tracker.get_value())))

        tt2 = 2 
        y2_val = Variable(tt2, 'x', num_decimal_places=3)
        on_screen_var_y2 = Variable(inv_func(tt2), MathTex("y_{i+1}"), num_decimal_places=2).next_to(on_screen_var_y1,DOWN*2)
        on_screen_var_y2.label.set_color(WHITE)
        on_screen_var_y2.value.set_color(T2)
        on_screen_var_y2.add_updater(lambda b: b.tracker.set_value(  inv_func(y2_val.tracker.get_value())))

        
        
        
        # self.add(line,Dot_x,Dot_dx,eqn1,eqn2,derivline)
        if player == True:
            self.play(FadeIn(line))
            self.play(FadeIn(derivline),FadeIn(diffeqn))
            self.wait(2)
            self.play(FadeIn(on_screen_var_y1),FadeIn(on_screen_var_y2))
            self.play(FadeIn(Dot_x),FadeIn(Dot_dx),FadeIn(Dot_deriv),FadeIn(brace_dx_group))
            self.play(FadeIn(func_OG),FadeIn(eqn1),FadeIn(eqn2))
            
            for ii in range(8):
                self.play(x.animate.increment_value(1),e.animate.increment_value(1),
                          y1_val.tracker.animate.set_value(ii+1),
                          y2_val.tracker.animate.set_value(ii+3))
                self.wait(1.3)
        else:
            None
            
class FE_Difference(Scene):
    def BackwardDiff(self):
        axes = Axes(
                x_range = [0,4,1],
                x_length = 3,
                y_range = [0,16,4],
                y_length = 1.5,
                axis_config=CONFIG2,
                ).to_edge(DL).shift(RIGHT).set_color(GREY)
        
        dt11 = Dot(axes.coords_to_point(1, 15), color=GREEN)
        dt12 = Dot(axes.coords_to_point(1, 0), color=GREEN)

        line1 = DashedLine(dt11, dt12, color=GREEN)
        
        dt21 = Dot(axes.coords_to_point(2, 12), color=GREEN)
        dt22 = Dot(axes.coords_to_point(2, 0), color=GREEN)

        line2 = DashedLine(dt21, dt22, color=GREEN)
        line3 = Line(dt11, dt21, color=GREEN)

        func_OG = axes.plot(lambda x: 16-x**2, x_range = [0,4], color = WHITE)
        
        val1 = MathTex(r"x-h").next_to(line1,DOWN).scale(0.5).set_color(GREEN).shift(UP*0.08)
        val2 = MathTex(r"x").next_to(line2,DOWN).scale(0.5).set_color(GREEN)

        Total_Group = VGroup(axes,dt11,line1,dt21,line2,func_OG,line3,val1,val2)
        
        return Total_Group
    
    
    def ForwardDiff(self):
        axes = Axes(
                x_range = [0,4,1],
                x_length = 3,
                y_range = [0,16,4],
                y_length = 1.5,
                axis_config=CONFIG2,
                ).to_edge(DL).shift(RIGHT).set_color(GREY)
        
        dt11 = Dot(axes.coords_to_point(2, 12), color=RED)
        dt12 = Dot(axes.coords_to_point(2, 0), color=RED)

        line1 = DashedLine(dt11, dt12, color=RED)
        
        dt21 = Dot(axes.coords_to_point(3, 7), color=RED)
        dt22 = Dot(axes.coords_to_point(3, 0), color=RED)

        line2 = DashedLine(dt21, dt22, color=RED)
        line3 = Line(dt11, dt21, color=RED)

        func_OG = axes.plot(lambda x: 16-x**2, x_range = [0,4], color = WHITE)
        
        val1 = MathTex(r"x").next_to(line1,DOWN).scale(0.5).set_color(RED)
        val2 = MathTex(r"x+h").next_to(line2,DOWN).scale(0.5).set_color(RED).shift(UP*0.08)

        Total_Group = VGroup(axes,dt11,line1,dt21,line2,func_OG,line3,val1,val2)
        
        return Total_Group
    
   
    def CentralDiff(self):
        axes = Axes(
                x_range = [0,4,1],
                x_length = 3,
                y_range = [0,16,4],
                y_length = 1.5,
                axis_config=CONFIG2,
                ).to_edge(DL).shift(RIGHT).set_color(GREY)
        
        dt11 = Dot(axes.coords_to_point(1, 15), color=BLUE)
        dt12 = Dot(axes.coords_to_point(1, 0), color=BLUE)

        line1 = DashedLine(dt11, dt12, color=BLUE)
        
        dt21 = Dot(axes.coords_to_point(3, 7), color=BLUE)
        dt22 = Dot(axes.coords_to_point(3, 0), color=BLUE)

        line2 = DashedLine(dt21, dt22, color=BLUE)
        line3 = Line(dt11, dt21, color=BLUE)

        func_OG = axes.plot(lambda x: 16-x**2, x_range = [0,4], color = WHITE)
        
        val1 = MathTex(r"x-h").next_to(line1,DOWN).scale(0.5).set_color(BLUE).shift(UP*0.08)
        val2 = MathTex(r"x+h").next_to(line2,DOWN).scale(0.5).set_color(BLUE).shift(UP*0.08)

        Total_Group = VGroup(axes,dt11,line1,dt21,line2,func_OG,line3,val1,val2)
        
        return Total_Group
    
    def MathScene1(self,player = True):
        """
        This scene develops the mathematical background for fintie differences, where
        the next scene will dive deeper into what these finite differences application looks like.

        """
        txt = Text("Finite Differences").to_edge(UP).set_color(GREEN)
        divider = Line([-4,0,0], [4,0,0]).next_to(txt,DOWN)

        TextExplain = Paragraph('Approximate derivatives discreetly', 'using difference quotients',line_spacing=1.2).scale(0.9).next_to(divider,DOWN)
        TextIDA = Paragraph("IDEA: ").next_to(TextExplain[0],LEFT).set_color(RED).scale(0.9)
        Group1Explain = VGroup(TextExplain,TextIDA)
        Matheqn1 =  MathTex(r"f'(x)", r" = \lim_{h\to 0} \frac{f(x+h) - f(x)}{h}").next_to(Group1Explain,DOWN*2).shift(RIGHT)

        BackwardDiff =  MathTex(r"f'(x) \approx", r" \frac{f(x) - f(x-h)}{h}").to_corner(DL).scale(0.7).shift(2*UP)
        BackwardDiff[1][9].set_color(ORANGE)
        BackwardDiff[1][12].set_color(ORANGE)

        CentralDiff =  MathTex(r"\approx", r" \frac{f(x+h) - f(x-h)}{2h}").next_to(BackwardDiff,RIGHT).scale(0.7).shift(LEFT*0.3)
        CentralDiff[1][4].set_color(ORANGE)
        CentralDiff[1][11].set_color(ORANGE)
        CentralDiff[1][15].set_color(ORANGE)

        ForwardDiff =  MathTex(r"\approx", r" \frac{f(x+h) - f(x)}{h}").next_to(CentralDiff,RIGHT*1).scale(0.7)
        ForwardDiff[1][4].set_color(ORANGE)
        ForwardDiff[1][12].set_color(ORANGE)

        BDText = Text("Backward Difference").next_to(BackwardDiff,UP).scale(0.7).set_color(GREEN)
        CDText = Text("Central Difference").next_to(CentralDiff,UP).scale(0.7).set_color(BLUE)
        FDText = Text("Foward Difference").next_to(ForwardDiff,UP*1).scale(0.7).shift(RIGHT*0.1).set_color(RED)

        Backward = VGroup(BDText,BackwardDiff)
        Central = VGroup(CDText,CentralDiff)
        Forward = VGroup(FDText,ForwardDiff)

        



        if player == True:
            self.play(FadeIn(txt))
            self.play(Write(divider))
            self.play(FadeIn(Group1Explain))
            self.wait(1)
            self.play(Write(Matheqn1))
            self.wait(2)
            self.play(FadeOut(Group1Explain))
                        
            Matheqn1.generate_target()
            Matheqn1.target.next_to(divider,DOWN)

            self.play(MoveToTarget(Matheqn1))
            self.wait(1)
            
            small =  Text("for small h").next_to(Matheqn1,DOWN).set_color(ORANGE).scale(0.7)

            
            self.play(FadeIn(small))
            self.play(FadeIn(Backward))
            self.play(FadeIn(Central))
            self.play(FadeIn(Forward))
            self.wait(1)
            BackGraph = self.BackwardDiff()
            BackGraph.next_to(Backward,DOWN)

            CentralGraph = self.CentralDiff()
            CentralGraph.next_to(Central,DOWN)

            ForwardGraph = self.ForwardDiff()
            ForwardGraph.next_to(Forward,DOWN)
            
            self.play(FadeIn(BackGraph),FadeIn(ForwardGraph),FadeIn(CentralGraph))

            BackwardTotal = VGroup(Backward,BackGraph)
            CentralTotal = VGroup(Central,CentralGraph)
            ForwardTotal = VGroup(Forward,ForwardGraph)
            centbox = SurroundingRectangle(CentralTotal, corner_radius=0.2)
            
            Matheqn2 =  MathTex(r"f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}").next_to(divider,DOWN)
            txt2 = Text("Central Differences").to_edge(UP).set_color(GREEN)
            
            Group1Fadeout = VGroup(BackwardTotal,ForwardTotal,CentralTotal,centbox,small)

            self.wait(2)
            self.play(Create(centbox))
            self.wait(1)

            self.play(ReplacementTransform(txt,txt2),TransformMatchingTex(Matheqn1,Matheqn2),FadeOut(Group1Fadeout, target_position=Matheqn2))

            return Matheqn2,txt2,divider

        else:
            Matheqn2 =  MathTex(r"f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}").next_to(divider,DOWN)
            txt2 = Text("Central Differences").to_edge(UP).set_color(GREEN)
            self.add(Matheqn2)
            self.add(txt2)
            self.add(divider)
            return Matheqn2,txt2,divider
        
        
    def MathScene2(self,Matheqn2,divider,player = True):
        
        Central_Diff_Discrete = MathTex(r"y'(x) \approx \frac{y_{i+1} - y_{i-1}}{2dx}").next_to(divider,DOWN)
        mattextx = MathTex(r"x = ")
        m0x = Matrix([[0,1,2,3,4]]).next_to(mattextx,RIGHT)
        matgroupx = VGroup(mattextx,m0x).next_to(Matheqn2,DOWN)
        
        mattexty = MathTex(r"y = ")
        m0y = Matrix([[0,1,4,9,16]]).next_to(mattexty,RIGHT)
        matgroupy = VGroup(mattexty,m0y).next_to(matgroupx,DOWN)
        
        
        arrow_1 = Arrow(buff = 0.1,start=DOWN/2, end=UP, color=BLUE).next_to(matgroupy,DOWN/2).shift(LEFT*0.81)
        arrow_2 = Arrow(buff = 0.1,start=DOWN/2, end=UP, color=BLUE).next_to(matgroupy,DOWN/2).shift(RIGHT*1.79)
        arrow_x = Arrow(buff = 0.1,start=LEFT, end=RIGHT/2, color=ORANGE).next_to(matgroupx,RIGHT/2)

        arrow1_math = MathTex(r"y_{i-1}").next_to(arrow_1,DOWN).set_color(BLUE)
        arrow2_math = MathTex(r"y_{i+1}").next_to(arrow_2,DOWN).set_color(BLUE)

        arrowx_math = MathTex(r"h = dx").next_to(arrow_x,RIGHT).set_color(ORANGE)

        arrow_group1 = VGroup(arrow_1,arrow_2,arrow1_math,arrow2_math)
        arrow_group2 = VGroup(arrow_x,arrowx_math)

        matreplacex = MathTex(r"x = x").next_to(Matheqn2,DOWN)
        matreplacey = MathTex(r"y = x^2").next_to(matreplacex,DOWN)


        self.play(FadeIn(matgroupx))
        self.play(FadeIn(matgroupy))
        self.wait(3)

        self.play(FadeIn(arrow_group1),Matheqn2[0][6:12].animate.set_color(BLUE),
                  Matheqn2[0][13:19].animate.set_color(BLUE))

        self.wait(2)
        self.play(FadeIn(arrow_group2),Matheqn2[0][21].animate.set_color(ORANGE))

        self.wait(3)

        self.play(FadeOut(arrow_group2,target_position=Matheqn2),
                  FadeOut(arrow_group1,target_position=Matheqn2),
                  TransformMatchingTex(Matheqn2,Central_Diff_Discrete))
        
        self.wait(2)
        
        self.play(ReplacementTransform(matgroupx,matreplacex),
          ReplacementTransform(matgroupy,matreplacey))
                
        
        self.play(ReplacementTransform(matgroupx,matreplacex),
          ReplacementTransform(matgroupy,matreplacey))
        
        return matreplacex,matreplacey,Central_Diff_Discrete
        
    def ApplicationScene1(self,MathScene2Group_Fadeout,MathScene2Group_Transform,player = True):
        T1 = MAROON
        T2 = RED
        DX_COL = BLUE
        
        OG_func = MathTex("y = x^2").to_edge(UR).shift(LEFT)

        x = ValueTracker(0)
        e = ValueTracker(1)
        dx = ValueTracker(2)
        Finite_Difference = MathTex(r"\frac{dy}{dx} =", r"\frac{y_{i+1} - y_{i-1}}{2dx}").next_to(OG_func,DOWN*2)
        Finite_Difference[1][10:13].set_color(DX_COL)
        Finite_Difference[1][0:4].set_color(T1)
        Finite_Difference[1][5:9].set_color(T2)

        axes = (
            Axes(
                x_range = [0,10,1],
                x_length = 8,
                y_range = [0,100,20],
                y_length = 2.4,
                axis_config=CONFIG,
                ).to_edge(DL).set_color(GREY)
            )
        

        grid_labels_OG = axes.get_axis_labels(
            Tex("x").scale(0.7), Tex("y").scale(0.7)
        ).shift(DOWN/4).set_color(LIGHT_GREY)
        
        func = axes.plot(lambda x: (x**2), x_range = [0,10], color = YELLOW)
        OG_funcbox = SurroundingRectangle(axes, corner_radius=0.2,buff = 0.4)
        

        
        secant = always_redraw(
            lambda: axes.get_secant_slope_group(
                x = x.get_value(),
                graph = func,
                dx = dx.get_value(),
                dx_line_color = DX_COL,
                dy_line_color = DX_COL,
                dx_label = "2dx",
                secant_line_color=TEAL,
                secant_line_length = 3,
                )
            )
        secant[2].scale(0.1)
        Dot_x = always_redraw(
            lambda: Dot(color = T1).scale(1).move_to(axes.c2p(x.get_value(), 
                                                    func.underlying_function(x.get_value()))))
        
        Dot_dx = always_redraw(
            lambda: Dot(color = T2).scale(1).move_to(axes.c2p(x.get_value() + dx.get_value(), 
                                                    func.underlying_function(x.get_value() + dx.get_value()))))
        
        # Dot_OG = always_redraw(
        #     lambda: Dot(color = WHITE).scale(1.2).move_to(axes.c2p(e.get_value(), 
        #                                             func.underlying_function(e.get_value()))))
        
        group_OG = Group(axes,OG_funcbox,func,Dot_x,Dot_dx,grid_labels_OG)

        axes_derivative = (
            Axes(
                x_range = [0,10,1],
                x_length = 8,
                y_range = [0,20,5],
                y_length = 2.4,
                axis_config=CONFIG,
                ).to_edge(UL).set_color(GREY)
            )
        
        grid_labels_deriv = axes_derivative.get_axis_labels(
            Tex("x").scale(0.7), MathTex(r"\frac{dy}{dx}").scale(0.7).set_color(LIGHT_GREY)
        ).shift(DOWN/4)
        
        derivative_func = always_redraw(lambda: 
                                        axes_derivative.plot(lambda x: (2*x),
                                                             x_range = [1,e.get_value()], color = YELLOW))
        Dot_deriv = always_redraw(
            lambda: Dot(color = WHITE).scale(1.2).move_to(axes_derivative.c2p(e.get_value(), 
                                                    derivative_func.underlying_function(e.get_value()))))
        
        Derivative_funcbox = SurroundingRectangle(axes_derivative, corner_radius=0.2,buff=0.4)
        group_deriv = Group(axes_derivative,Derivative_funcbox,derivative_func,Dot_deriv,grid_labels_deriv)
        
        DL1 = always_redraw( lambda: Dot(color = WHITE, radius=0.00001).move_to(axes.c2p(e.get_value(),28.5)))
        DL2 = always_redraw( lambda: Dot(color = WHITE, radius=0.00001).move_to(axes.c2p(e.get_value(),0 )))
        #Deriv_Line = always_redraw(lambda: DashedLine(Dot_OG,Dot_deriv))
        Deriv_Line = always_redraw(lambda: DashedLine(DL1,DL2))

        tt1 = 0 
        y1_val = Variable(tt1, 'x', num_decimal_places=3)
        on_screen_var_y1 = Variable(function(tt1), MathTex("y_{i-1}"), num_decimal_places=2).next_to(Finite_Difference,DOWN*2)
        on_screen_var_y1.label.set_color(WHITE)
        on_screen_var_y1.value.set_color(T1)
        on_screen_var_y1.add_updater(lambda a: a.tracker.set_value(function(y1_val.tracker.get_value())))

        tt2 = 2 
        y2_val = Variable(tt1, 'x', num_decimal_places=3)
        on_screen_var_y2 = Variable(function(tt2), MathTex("y_{i+1}"), num_decimal_places=2).next_to(on_screen_var_y1,DOWN*2)
        on_screen_var_y2.label.set_color(WHITE)
        on_screen_var_y2.value.set_color(T2)
        on_screen_var_y2.add_updater(lambda b: b.tracker.set_value(  function(y2_val.tracker.get_value())))

        
        tt3 = 1 
        ans_val = Variable(tt1, 'x', num_decimal_places=3)
        on_screen_var_ans = Variable( deriv(tt3), MathTex(r"\frac{dy}{dx}"), num_decimal_places=2).next_to(on_screen_var_y2,DOWN*2)
        on_screen_var_ans.label.set_color(WHITE)
        on_screen_var_ans.value.set_color(WHITE)
        on_screen_var_ans.add_updater(lambda c: c.tracker.set_value(  deriv(ans_val.tracker.get_value())))

        Derivative_func = MathTex(r"\frac{dy}{dx} = 2x").next_to(on_screen_var_ans,DOWN*3)

        if player == True:
            # self.add(Deriv_Line,DL1,DL2)
            self.play(FadeOut(MathScene2Group_Fadeout),
                      ReplacementTransform(MathScene2Group_Transform[0],OG_func),
                      ReplacementTransform(MathScene2Group_Transform[1],Finite_Difference),
                      FadeIn(Derivative_func))
            self.wait(3)
            self.play(FadeIn(secant),
                      FadeIn(on_screen_var_y1,on_screen_var_y2,on_screen_var_ans),
                      FadeIn(group_OG,group_deriv),
                      FadeIn(y1_val.tracker.set_value(0),y2_val.tracker.set_value(2),ans_val.tracker.set_value(1)))
            # self.play(x.animate.set_value(8),e.animate.set_value(9),y1_val.animate.set_value(8),y2_val.animate.set_value(10))
            self.wait(4)
            for ii in range(8):
                self.play(x.animate.increment_value(1),y1_val.tracker.animate.set_value(ii+1),y2_val.tracker.animate.set_value(ii+3))
                self.play(e.animate.increment_value(1),ans_val.tracker.animate.set_value(ii+2))
                self.wait(1)
        else:
            None
    
    def construct(self):
        Matheqn2,txt2,divider = self.MathScene1(True)
        matreplacex,matreplacey,Central_Diff_Discrete = self.MathScene2(Matheqn2,divider,True)
        MathScene2Group_Fadeout = VGroup(divider,matreplacex,txt2)
        MathScene2Group_Transform = VGroup(matreplacey, Central_Diff_Discrete)
        self.ApplicationScene1(MathScene2Group_Fadeout,MathScene2Group_Transform,True)

        
        

        