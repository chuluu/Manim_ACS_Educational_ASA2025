# -*- coding: utf-8 -*-
"""
Created on Sat May 13 10:55:52 2023

@author: mbluu
"""

from manim import *

# class FEA_Derivation(Scene):
#     def construct(self):
#         pt1 = MathTex(
#             r"\frac{d^2 \xi}{dx^2} + \frac{1}{c^2}\frac{d^2 \xi}{dt^2} = 0"
#         )

#         self.play(Write(pt1))
#         self.wait(2)

#         pt1[0][0:7].set_color(YELLOW)
#         self.play(Transform(pt1,pt1))
        
#         pt2 = MathTex(
#             r"\frac{d^2 \xi}{dx^2} + \frac{1}{c^2}\frac{d^2 \xi}{dt^2} = 0"
#         ).next_to(pt1)
#         self.play(FadeIn(pt2))

def funMatLabMultip(f,t):
    """create an n x m array from 2 vectors of size n and m.
    Resulting rows are the multiplication of each element of the first vector for all the elements of the second vector
    f=np.array([2,4])
    t=np.array([1,2,3,4,5])
    [[ 2  4  6  8 10]
     [ 4  8 12 16 20]]
    """
    if t.size==t.shape[0]:
        k=f[0]*t
        for i in f[1:]:
            j=i*t
            k=np.vstack((k,j))
    else:
        raise Exception('arrays should 1D arrays')
    return k


class ColorGradients(ThreeDScene):
    def Create_Func(self):
        import numpy as np
        import matplotlib.pyplot as plt
        
        Nx = 101
        dx = 1
        
        f  = 10
        Nt = 1001
        dt = 1/(f*100)
        c  = 343
        
        
        x  = np.linspace(0,(Nx-1)*dx,Nx)
        t  = np.linspace(0,(Nt-1)*dt,Nt)
        
        C  = c*(dt/dx)
        
        U  = np.zeros((Nt,Nx))
        s1 = int(Nt/f)
        
        U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
        U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
        # U[:,0] = np.sin(2*np.pi*f*t)

        
        
        for jj in range(2,Nt):
            for nn in range(1,Nx-1):
                U1 = 2*U[jj-1,nn] - U[jj-2,nn]
                U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
                U[jj,nn] = U1 + C*C*U2
                U[:,-1] = U[:,-2]
        U = np.flip(U, 1)

        U_arr = []
        for ii in range(Nt):
            U_arr.append(funMatLabMultip(U[ii,:],np.ones(Nx))) 

        
        return U_arr
    
    def Func_3d(self,U, t, u, v):
        x = int(u)
        y = int(v)
        t = int(t)

        return U[t][x,y]
    
    def Func_3d_Plane(self,u,v):
        x = u
        y = v
        
        return 0
    
    def construct(self):
        U_arr = self.Create_Func()
        axes = ThreeDAxes(x_range=(0, 100, 1), y_range=(0, 100, 1), z_range=(-2, 2, 0.1), x_length=10, y_length = 5, z_length = 5)
        resolution_fa = 1 # Set for more point resolution on the wave plot
        time = ValueTracker(0)
        phi, theta, focal_distance, gamma, distance_to_origin = self.camera.get_value_trackers()

        surface_plane = always_redraw(lambda : Surface(
            lambda u, v: axes.c2p(u, v, self.Func_3d(U_arr,time.get_value(),u,v)),
            resolution=(resolution_fa, resolution_fa),
            u_range=[0, 100],
            v_range=[0, 100],
            should_make_jagged = True, 
        ))

        BC = Surface(
            lambda u, v: axes.c2p(self.Func_3d_Plane(u,v), u, v),
            resolution=(1, 1),
            u_range=[0, 100],
            v_range=[-2.2, 2.2],
            should_make_jagged = True, 
            checkerboard_colors  = False,
            fill_opacity  = 1,
            fill_color = GRAY
        )

        Fixed3d = MathTex(r'Fixed: \frac{d\xi}{dx} = 0').scale(1).shift(LEFT*5).shift(DOWN*3)
        self.add_fixed_in_frame_mobjects(Fixed3d) # <---- Add this line
        Free3d =  MathTex(r'Free: \xi = 0').scale(1).shift(RIGHT*5).shift(DOWN*3)
        self.add_fixed_in_frame_mobjects(Free3d) # <---- Add this line

        prismSmall = Prism(dimensions=[5,10,3],fill_color=LIGHT_GRAY,stroke_color=BLACK).rotate(PI / 2)
        prismSmall.set_style(fill_opacity=0.3)

        surface_plane.set_style(fill_opacity=1)
        surface_plane = always_redraw(lambda : surface_plane.set_fill_by_value(axes=axes, colorscale=[ (BLUE,-2), (BLUE,-1), (BLUE,-0.5), (GREEN, -0.1), (GREEN, 0), (GREEN, 0.1), (RED, 0.5), (RED, 1), (RED, 2)], axis=2))

        BigSpace_PlotGroup = Group(prismSmall,Fixed3d,Free3d,BC,surface_plane,axes).scale(0.1).move_to([-4,-4,0])
        BigSpace_PlotGroup.generate_target()
        BigSpace_PlotGroup.target.scale(10).move_to([0,0,0])
        
        self.set_camera_orientation(phi=65 * DEGREES, theta=-130 * DEGREES,distance=1)
        self.add(BigSpace_PlotGroup)
        self.play(MoveToTarget(BigSpace_PlotGroup))
        # self.play(ScaleInPlace(BigSpace_PlotGroup,2))
        #self.begin_ambient_camera_rotation(rate=0.2)
        prism2 = Prism(dimensions=[5,10,3],fill_color=WHITE,stroke_color=BLACK).rotate(PI / 2)

        self.wait(3)
        self.begin_ambient_camera_rotation(rate=10.0)
        self.play(ScaleInPlace(BigSpace_PlotGroup,0),runtime=2)
        self.play(ReplacementTransform(BigSpace_PlotGroup,prism2))
        self.stop_ambient_camera_rotation()
        self.wait(2)

        
        # self.begin_ambient_camera_rotation(rate=0.1)

        # self.play(time.animate.set_value(333),run_time=3,rate_func=linear)
        # self.stop_ambient_camera_rotation()
        
        # self.play(theta.animate.set_value(-90*DEGREES), phi.animate.set_value(90*DEGREES),run_time=0.5,rate_func=linear)
        # self.play(time.animate.set_value(666),run_time=3,rate_func=linear)

        # self.play(theta.animate.set_value(-90*DEGREES), phi.animate.set_value(0*DEGREES),run_time=0.5,rate_func=linear)

        # self.remove(BC)
               
        # BC = Surface(
        #     lambda u, v: axes.c2p(self.Func_3d_Plane(u,v), u, v),
        #     resolution=(1, 1),
        #     u_range=[0, 100],
        #     v_range=[0, 2.2],
        #     should_make_jagged = True, 
        #     checkerboard_colors  = False,
        #     fill_opacity  = 1,
        #     fill_color = GRAY
        # )
        # self.add(BC)
        # self.play(time.animate.set_value(999),run_time=4,rate_func=linear)

        # self.wait(2)
        
class FE_Difference_3D(ThreeDScene):
    def Create_Func(self):
        import numpy as np
        import matplotlib.pyplot as plt
        
        Nx = 101
        dx = 1
        
        f  = 10
        Nt = 1001
        dt = 1/(f*100)
        c  = 343
        
        
        x  = np.linspace(0,(Nx-1)*dx,Nx)
        t  = np.linspace(0,(Nt-1)*dt,Nt)
        
        C  = c*(dt/dx)
        
        U  = np.zeros((Nt,Nx))
        s1 = int(Nt/f)
        
        # U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
        # U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
        
        U[:,0] = np.sin(2*np.pi*f*t)
        
        for jj in range(2,Nt):
            for nn in range(1,Nx-1):
                U1 = 2*U[jj-1,nn] - U[jj-2,nn]
                U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
                U[jj,nn] = U1 + C*C*U2
                U[:,-1] = U[:,-2]
        U = np.flip(U, 1)

        U_arr = []
        for ii in range(Nt):
            U_arr.append(funMatLabMultip(U[ii,:],np.ones(Nx))) 

        
        return U_arr
    
    def Func_3d(self,U, t, u, v):
        x = int(u)
        y = int(v)
        t = int(t)

        return U[t][x,y]
    
    def Func_3d_Plane(self,u,v):
        x = u
        y = v
        
        return 0
    
    def construct(self):
        length = 20
        U_arr = self.Create_Func()
        axes = ThreeDAxes(x_range=(0, 100, 1), y_range=(0, 100, 1), z_range=(-2, 2, 0.1), x_length=length, y_length = 5, z_length = 5)
        resolution_fa = 100 # Set for more point resolution on the wave plot
        time = ValueTracker(0)
        phi, theta, focal_distance, gamma, distance_to_origin = self.camera.get_value_trackers()

        surface_plane = always_redraw(lambda : Surface(
            lambda u, v: axes.c2p(u, v, self.Func_3d(U_arr,time.get_value(),u,v)),
            resolution=(resolution_fa, resolution_fa),
            u_range=[0, 100],
            v_range=[0, 100],
            should_make_jagged = True, 
        ))

        BC = Surface(
            lambda u, v: axes.c2p(self.Func_3d_Plane(u,v), u, v),
            resolution=(1, 1),
            u_range=[0, 100],
            v_range=[-2.2, 2.2],
            should_make_jagged = True, 
            checkerboard_colors  = False,
            fill_opacity  = 1,
            fill_color = GRAY
        )

        Fixed3d = MathTex(r'Fixed: \frac{d\xi}{dx} = 0').scale(1).shift(LEFT*(int(length/2))).shift(DOWN*3)
        self.add_fixed_in_frame_mobjects(Fixed3d) # <---- Add this line
        Free3d =  MathTex(r'Free: \xi = 0').scale(1).shift(RIGHT*(int(length/2))).shift(DOWN*3)
        self.add_fixed_in_frame_mobjects(Free3d) # <---- Add this line

        prismSmall = Prism(dimensions=[5,length,3],fill_color=LIGHT_GRAY,stroke_color=BLACK).rotate(PI / 2)
        prismSmall.set_style(fill_opacity=0.3)

        surface_plane.set_style(fill_opacity=1)
        surface_plane = always_redraw(lambda : surface_plane.set_fill_by_value(axes=axes, colorscale=[ (BLUE,-2), (BLUE,-1), (BLUE,-0.5), (GREEN, -0.1), (GREEN, 0), (GREEN, 0.1), (RED, 0.5), (RED, 1), (RED, 2)], axis=2))
   
        BigSpace_PlotGroup = Group(prismSmall,Fixed3d,Free3d,BC,surface_plane,axes).scale(0.4).move_to([0,2,0])
        
        self.set_camera_orientation(phi=90 * DEGREES, theta=-90 * DEGREES,distance=1)
        #self.set_camera_orientation(phi=90 * DEGREES, theta=-90 * DEGREES,distance=1)
        self.add(BigSpace_PlotGroup)

        # self.play(ScaleInPlace(BigSpace_PlotGroup,2))
        
        self.wait(2)
        self.play(time.animate.set_value(500),run_time=8,rate_func=linear)
        
        # self.play(theta.animate.set_value(-90*DEGREES), phi.animate.set_value(90*DEGREES),run_time=0.5,rate_func=linear)
        # self.play(time.animate.set_value(666),run_time=2,rate_func=linear)

        self.play(theta.animate.set_value(-90*DEGREES), phi.animate.set_value(0*DEGREES),run_time=0.2,rate_func=linear)

        self.remove(BC)
               
        BC = Surface(
            lambda u, v: axes.c2p(self.Func_3d_Plane(u,v), u, v),
            resolution=(1, 1),
            u_range=[0, 100],
            v_range=[0, 2.2],
            should_make_jagged = True, 
            checkerboard_colors  = False,
            fill_opacity  = 1,
            fill_color = GRAY
        )
        self.add(BC)
        self.play(time.animate.set_value(999),run_time=8,rate_func=linear)

        # self.wait(2)

class FE_Difference(Scene):
    def Create_Func(self):
        import numpy as np
        import matplotlib.pyplot as plt
        
        Nx = 101
        dx = 1
        
        f  = 10
        Nt = 1001
        dt = 1/(f*100)
        c  = 343
        
        
        x  = np.linspace(0,(Nx-1)*dx,Nx)
        t  = np.linspace(0,(Nt-1)*dt,Nt)
        
        C  = c*(dt/dx)
        
        U  = np.zeros((Nt,Nx))
        s1 = int(Nt/f)
        
        U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
        U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
        
        
        
        for jj in range(2,Nt):
            for nn in range(1,Nx-1):
                U1 = 2*U[jj-1,nn] - U[jj-2,nn]
                U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
                U[jj,nn] = U1 + C*C*U2
                U[:,-1] = U[:,-2]
        U = np.flip(U, 1)

        return U

    def Fun_2D(self,U, t, u):
        x = int(u)
        t = int(t)

        return U[t,x]
    
    
    
    def Math_Text(self):
        pt1 = MathTex(
            r"\frac{d^2 \xi}{dx^2}"
        )
        
        pt2 = MathTex(
           r"\frac{1}{c^2}\frac{d^2 \xi}{dt^2}"
        )
        
        pt3 = MathTex(" = 0")
        variables_group = VGroup(pt1,pt2,pt3).arrange_submobjects().shift(UP)

        eqn1 = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn2 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn3 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", " = 0")
        eqn4 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "=", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", "")
        eqn5 = MathTex(r"\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}", "=", r" \left(\frac{\Delta x^2}{\Delta t^2c^2}\right)", r"\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}")
        text1 = MathTex("CFL")
        text1.set_color(YELLOW)
        text1.next_to(eqn5[2],DOWN)
        eqn5[2].set_color(YELLOW)
        group1 = VGroup(eqn5,text1)

        self.play(Write(eqn1))


        eqn1[0].set_color(YELLOW)
        eqn2[0].set_color(YELLOW)

        self.wait(2)

        self.play(Transform(eqn1,eqn1))

        self.play(TransformMatchingTex(eqn1,eqn2))
        self.wait(2)
        eqn2[0].set_color(WHITE)
        eqn2[2].set_color(YELLOW)
        eqn3[2].set_color(YELLOW)

        self.play(TransformMatchingTex(eqn2,eqn2))
        self.play(TransformMatchingTex(eqn2,eqn3))
        self.wait(2)
        self.play(TransformMatchingTex(eqn3,eqn4))
        self.wait(2)
        self.play(TransformMatchingTex(eqn4,group1))
        self.wait(2)
        self.play(group1.animate.shift(3*UP))
        
    def construct(self):
        
        self.Math_Text()
        
        U = self.Create_Func()
        resolution_fa = 25 # Set for more point resolution on the wave plot
        time = ValueTracker(100)
        axes = Axes(x_range = [0,100,10],
                     x_length = 20,
                     y_range = [-2,2,1],
                     y_length= 5)
        
        # y_label = axes.get_y_axis_label(MathTex(r"\xi"), edge=LEFT, direction=LEFT, buff=0.4)
        # x_label = axes.get_x_axis_label("x (position)")
        # grid_labels = VGroup(x_label, y_label)

        graph = always_redraw(lambda:
                              axes.plot(lambda x: self.Fun_2D(U,time.get_value(),x),
                                        x_range=[0,100],
                                        color = YELLOW,
                                        use_smoothing=False))
            
        text = Text("Gradient Color")
        text.set_color_by_gradient(BLUE,GREEN,RED).set_sheen_direction(UP)
        graph = always_redraw(lambda: graph.set_color_by_gradient(BLUE,GREEN,RED).set_sheen_direction(UP) )
        Total_2DGroup = Group(axes,graph).move_to([0,-2,0]).scale(0.4)
        self.play(Create(graph),Create(axes),Create(text))
        self.wait(2)
        self.play(time.animate.set_value(1000),run_time=10,rate_func=linear)
        

class FE_DifferenceTot(ThreeDScene):
    def Finite_Difference_Total_Code(self,signal = 'sine'):
        import numpy as np
        import matplotlib.pyplot as plt
        
        Nx = 101
        dx = 1
        
        f  = 10
        Nt = 1001
        dt = 1/(f*100)
        c  = 343
        
        
        x  = np.linspace(0,(Nx-1)*dx,Nx)
        t  = np.linspace(0,(Nt-1)*dt,Nt)
        
        C  = c*(dt/dx)
        
        U  = np.zeros((Nt,Nx))
        s1 = int(Nt/f)
        
        # U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
        # U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
        if signal == 'sine':
            U[:,0] = np.sin(2*np.pi*f*t)
            
        elif signal == 'wave':
            U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
            U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
            
        elif signal == 'pulse':
            U[0:s1,0] = np.ones(s1)
            U[0:s1,1] = np.ones(s1)
        else:
            print('choose proper')
        
        for jj in range(2,Nt):
            for nn in range(1,Nx-1):
                U1 = 2*U[jj-1,nn] - U[jj-2,nn]
                U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
                U[jj,nn] = U1 + C*C*U2
                U[:,-1] = U[:,-2]
        U = np.flip(U, 1)
        
        return U
    
    def Create_Func2D(self,U):
        return U
    
    def Create_Func3D(self,U):
        Nx = len(U[1,:])
        Nt = len(U)

        U_arr = []
        for ii in range(Nt):
            U_arr.append(funMatLabMultip(U[ii,:],np.ones(Nx))) 

        
        return U_arr
    
    def Fun_2D(self,U, t, u):
        x = int(u)
        t = int(t)

        return U[t,x]
    
    def Func_3d(self,U, t, u, v):
        x = int(u)
        y = int(v)
        t = int(t)

        return U[t][x,y]
    
    def Func_3d_Plane(self,u,v):
        x = u
        y = v
        
        return 0
    
    def Math_Text(self):
        pt1 = MathTex(r"\frac{d^2 \xi}{dx^2}")
        
        pt2 = MathTex(r"\frac{1}{c^2}\frac{d^2 \xi}{dt^2}")
        
        pt3 = MathTex(" = 0")
        variables_group = VGroup(pt1,pt2,pt3).arrange_submobjects().shift(UP)

        eqn1 = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn2 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn3 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", " = 0")
        eqn4 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "=", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", "")
        eqn5 = MathTex(r"\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}", "=", r" \left(\frac{\Delta x^2}{\Delta t^2c^2}\right)", r"\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}")
        eqn6 = MathTex(r"\xi_{j}^{n+1}", "=", r"\left(\frac{c \Delta t}{\Delta x}\right)^2", r"\left(\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}\right)", "+" , r" 2 \xi_{j}^{n} - \xi_{j}^{n-1}")

        text1 = MathTex("CFL")
        text1.set_color(GREEN)
        text1.next_to(eqn6[2],DOWN)
        eqn6[2].set_color(GREEN)
        group1 = VGroup(eqn6,text1)
        group1.scale(0.5).shift(UP)

        self.play(Write(eqn1))


        eqn1[0].set_color(YELLOW)
        eqn2[0].set_color(YELLOW)

        self.wait(2)

        self.play(Transform(eqn1,eqn1))

        self.play(TransformMatchingTex(eqn1,eqn2))
        self.wait(2)
        eqn2[0].set_color(WHITE)
        eqn2[2].set_color(YELLOW)
        eqn3[2].set_color(YELLOW)

        self.play(TransformMatchingTex(eqn2,eqn2))
        self.play(TransformMatchingTex(eqn2,eqn3))
        self.wait(2)
        self.play(TransformMatchingTex(eqn3,eqn4))
        self.wait(2)
        self.play(TransformMatchingTex(eqn4,eqn5))
        self.wait(2)
        self.play(TransformMatchingTex(eqn5,group1))

        #self.play(group1.animate.shift(3*UP))
        
    def construct(self):
        # 2D Setup
        length = 20
        runtime = 5
        resolution_fa = 1 # Set for more point resolution on the wave plot

        axes = Axes(x_range = [0,100,10],
                     x_length = length,
                     y_range = [-2,2,1],
                     y_length= 5)
        axes3d = ThreeDAxes(x_range=(0, 100, 10), y_range=(0, 100, 10), z_range=(-2, 2, 1), x_length=length, y_length = 5, z_length = 5)


        time = ValueTracker(0)

        U_OG = self.Finite_Difference_Total_Code()
        U_2D = self.Create_Func2D(U_OG)
        U_3D = self.Create_Func3D(U_OG)

        graph = always_redraw(lambda:
                              axes.plot(lambda x: self.Fun_2D(U_2D,time.get_value(),x),
                                        x_range=[0,100],
                                        color = GREEN,
                                        use_smoothing=False))
        Force2dArrow = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(axes,RIGHT)
        F2d = MathTex(r"F = sin(\omega t)",color=RED).next_to(Force2dArrow,RIGHT)
        Force2d = Group(Force2dArrow,F2d)
        Total_2DGroup = Group(axes,graph,Force2d).move_to([0,-2,0]).scale(0.4)


        
        # 3D Set up
        phi, theta, focal_distance, gamma, distance_to_origin = self.camera.get_value_trackers()

        surface_plane = always_redraw(lambda : Surface(
            lambda u, v: axes3d.c2p(u, v, self.Func_3d(U_3D,time.get_value(),u,v)),
            resolution=(resolution_fa, resolution_fa),
            u_range=[0, 100],
            v_range=[0, 100],
            should_make_jagged = True, 
        ))

        BC = Surface(
            lambda u, v: axes3d.c2p(self.Func_3d_Plane(u,v), u, v),
            resolution=(1, 1),
            u_range=[0, 100],
            v_range=[-2.2, 2.2],
            should_make_jagged = True, 
            checkerboard_colors  = False,
            fill_opacity  = 1,
            fill_color = GRAY
        )

        Fixed3d = MathTex(r'Fixed: \frac{d\xi}{dx} = 0').scale(1.2).shift(LEFT*(int(length/2))).shift(DOWN*5)
        Free3d =  MathTex(r'Free: \xi = 0').scale(1.2).shift(RIGHT*(int(length/2))).shift(DOWN*5)
        Force3dArrow = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(axes3d,RIGHT)
        F3d = MathTex(r"F = sin(\omega t)",color=RED).next_to(Force3dArrow,RIGHT)
        Force3d = Group(Force3dArrow,F3d)
        
        prismSmall = Prism(dimensions=[5,length,3],fill_color=LIGHT_GRAY,stroke_color=BLACK).rotate(PI / 2)
        prismSmall.set_style(fill_opacity=0.3)

        surface_plane.set_style(fill_opacity=1)
        surface_plane = always_redraw(lambda : surface_plane.set_fill_by_value(axes=axes3d, colorscale=[ (BLUE,-2), (BLUE,-1), (BLUE,-0.5), (GREEN, -0.1), (GREEN, 0), (GREEN, 0.1), (RED, 0.5), (RED, 1), (RED, 2)], axis=2))
   
        BigSpace_PlotGroup = Group(prismSmall,Fixed3d,Free3d,BC,surface_plane,axes3d,Force3d).scale(0.4).move_to([0,2,0])
        
        tt = 0
        on_screen_var = Variable(tt, Text("Time"), num_decimal_places=3).scale(0.5).shift(LEFT).shift(DOWN/2)
        time_var_text = Text("s").next_to(on_screen_var,RIGHT)
        on_screen_var.label.set_color(WHITE)
        on_screen_var.value.set_color(GREEN)
        var_tracker = on_screen_var.tracker
        
        
        self.Math_Text()
        
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES,distance=1)
        # self.set_camera_orientation(phi=90 * DEGREES, theta=-90 * DEGREES,distance=1)


        # self.add_fixed_in_frame_mobjects(Fixed3d) # <---- Add this line
        # self.add_fixed_in_frame_mobjects(Free3d) # <---- Add this line
        
        # self.add_fixed_in_frame_mobjects(Total_2DGroup) # <---- Add this line
        # self.add_fixed_in_frame_mobjects(graph) # <---- Add this line
        
        self.play(FadeIn(BigSpace_PlotGroup),FadeIn(Total_2DGroup),FadeIn(on_screen_var))
        self.wait(2)
        self.play(time.animate.set_value(1000),var_tracker.animate.set_value(runtime),run_time=runtime,rate_func=linear)

class FE_DifferenceTot_ChangeSignals(ThreeDScene):
    def Finite_Difference_Total_Code(self,signal = 0):
        import numpy as np
        import matplotlib.pyplot as plt
        print(signal)
        Nx = 101
        dx = 1
        
        f  = 10
        Nt = 1001
        dt = 1/(f*100)
        c  = 343
        
        
        x  = np.linspace(0,(Nx-1)*dx,Nx)
        t  = np.linspace(0,(Nt-1)*dt,Nt)
        
        C  = c*(dt/dx)
        
        U  = np.zeros((Nt,Nx))
        s1 = int(Nt/f)
        
        # U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
        # U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
        if signal == 'wave':
            U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
            U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
            
        elif signal == 'sine':
            U[:,0] = np.sin(2*np.pi*f*t)
            
            
        elif signal == 'pulse':
            U[0:s1,0] = np.ones(s1)
            U[0:s1,1] = np.ones(s1)
        else:
            print('choose proper')
        
        for jj in range(2,Nt):
            for nn in range(1,Nx-1):
                U1 = 2*U[jj-1,nn] - U[jj-2,nn]
                U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
                U[jj,nn] = U1 + C*C*U2
                U[:,-1] = U[:,-2]
        U = np.flip(U, 1)
        
        return U
    
    def Create_Func2D(self,U):
        return U
    
    def Create_Func3D(self,U):
        Nx = len(U[1,:])
        Nt = len(U)

        U_arr = []
        for ii in range(Nt):
            U_arr.append(funMatLabMultip(U[ii,:],np.ones(Nx))) 

        
        return U_arr
    
    def Func_2d(self,U, t, u):
        x = int(u)
        t = int(t)

        return U[t,x]
    
    def Func_3d(self,U, t, u, v):
        x = int(u)
        y = int(v)
        t = int(t)

        return U[t][x,y]
    
    def Func_3d_Plane(self,u,v):
        x = u
        y = v
        
        return -0.2
    
    def Math_Text(self):
        pt1 = MathTex(r"\frac{d^2 \xi}{dx^2}")
        
        pt2 = MathTex(r"\frac{1}{c^2}\frac{d^2 \xi}{dt^2}")
        
        pt3 = MathTex(" = 0")
        variables_group = VGroup(pt1,pt2,pt3).arrange_submobjects().shift(UP)

        eqn1 = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn2 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn3 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", " = 0")
        eqn4 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "=", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", "")
        eqn5 = MathTex(r"\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}", "=", r" \left(\frac{\Delta x^2}{\Delta t^2c^2}\right)", r"\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}")
        eqn6 = MathTex(r"\xi_{j}^{n+1}", "=", r"\left(\frac{c \Delta t}{\Delta x}\right)^2", r"\left(\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}\right)", "+" , r" 2 \xi_{j}^{n} - \xi_{j}^{n-1}")

        text1 = MathTex("CFL")
        text1.set_color(GREEN)
        text1.next_to(eqn6[2],DOWN)
        eqn6[2].set_color(GREEN)
        group1 = VGroup(eqn6,text1)
        group1.scale(0.5).shift(UP).shift(LEFT/2)

        self.play(Write(eqn1))


        eqn1[0].set_color(YELLOW)
        eqn2[0].set_color(YELLOW)

        self.wait(2)

        self.play(Transform(eqn1,eqn1))

        self.play(TransformMatchingTex(eqn1,eqn2))
        self.wait(2)
        eqn2[0].set_color(WHITE)
        eqn2[2].set_color(YELLOW)
        eqn3[2].set_color(YELLOW)

        self.play(TransformMatchingTex(eqn2,eqn2))
        self.play(TransformMatchingTex(eqn2,eqn3))
        self.wait(2)
        self.play(TransformMatchingTex(eqn3,eqn4))
        self.wait(2)
        self.play(TransformMatchingTex(eqn4,eqn5))
        self.wait(2)
        self.play(TransformMatchingTex(eqn5,group1))

        #self.play(group1.animate.shift(3*UP))

    
    def Create_Plots(self,resolution_fa,length,signal):
        time = ValueTracker(0)
        U_OG = self.Finite_Difference_Total_Code(signal = signal)
        U_2D = self.Create_Func2D(U_OG)
        U_3D = self.Create_Func3D(U_OG)
        
        
        axes = Axes(x_range = [0,100,10],
                         x_length = length,
                         y_range = [-2,2,1],
                         y_length= 5,
                         axis_config={"include_numbers": True,
                                      "include_tip": False})
            
        axes3d = ThreeDAxes(x_range=(0, 100, 10), 
                            y_range=(0, 100, 10), 
                            z_range=(-2, 2, 1), 
                            x_length=length, 
                            y_length = 5, 
                            z_length = 5,
                            axis_config={"include_numbers": True,
                                         "include_tip": False})
    
        graph = always_redraw(lambda: axes.plot(lambda x: self.Func_2d(U_2D,time.get_value(),x),
                                            x_range=[0,100],
                                            color = GREEN,
                                            use_smoothing=False))
    

        # 3D Set up    
        surface_plane = always_redraw(lambda : Surface(
            lambda u, v: axes3d.c2p(u, v, self.Func_3d(U_3D,time.get_value(),u,v)),
            resolution=(resolution_fa, resolution_fa),
            u_range=[0, 100],
            v_range=[0, 100],
            should_make_jagged = True, 
        ))

        surface_plane.set_style(fill_opacity=1)
        surface_plane = always_redraw(lambda : surface_plane.set_fill_by_value(axes=axes3d, colorscale=[ (BLUE,-2), (BLUE,-1), (BLUE,-0.5), (GREEN, -0.1), (GREEN, 0), (GREEN, 0.1), (RED, 0.5), (RED, 1), (RED, 2)], axis=2))
        
        # Setup 2D Variables
        Force2dArrow = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(axes,RIGHT*1.1)
        Force3dArrow = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(axes3d,RIGHT*1.1)

        if signal == 'sine':
            F2d = MathTex(r"F = sin(\omega t)",color=RED).next_to(Force2dArrow,RIGHT)
            F3d = MathTex(r"F = sin(\omega t)",color=RED).next_to(Force3dArrow,RIGHT)
        elif signal == 'wave':
            F2d = MathTex(r"F = sin(\omega t)",color=RED).next_to(Force2dArrow,RIGHT)
            F3d = MathTex(r"F = sin(\omega t)",color=RED).next_to(Force3dArrow,RIGHT)
        elif signal == 'pulse':
            F2d = MathTex(r"F = rect(t)",color=RED).next_to(Force2dArrow,RIGHT)
            F3d = MathTex(r"F = rect(t)",color=RED).next_to(Force3dArrow,RIGHT)    
        
        Force2d = Group(Force2dArrow,F2d)
        Force3d = Group(Force3dArrow,F3d)
        
        
        Fixed3d = MathTex(r'Fixed: \frac{d\xi}{dx} = 0').scale(1.2).shift(LEFT*(int(length/2))).shift(DOWN*5)
        Free3d =  MathTex(r'Free: \xi = 0').scale(1.2).shift(RIGHT*(int(length/2))).shift(DOWN*5)
        prismSmall = Prism(dimensions=[5,length,3],fill_color=LIGHT_GRAY,stroke_color=BLACK).rotate(PI / 2)
        prismSmall.set_style(fill_opacity=0.3)
        
        BC = Surface(
            lambda u, v: axes3d.c2p(self.Func_3d_Plane(u,v), u, v),
            resolution=(1, 1),
            u_range=[0, 100],
            v_range=[0, 2.2],
            should_make_jagged = True, 
            checkerboard_colors  = False,
            fill_opacity  = 1,
            fill_color = GRAY
        )
        Extra_3D_Decorators = Group(BC,Fixed3d,Free3d,prismSmall)

         
        Total_2DGroup = Group(axes,Force2d,graph).move_to([0,-2,0]).scale(0.4)
        Total_3DGroup = Group(axes3d,Force3d,surface_plane,Extra_3D_Decorators).move_to([0,2,0]).scale(0.4)

        return [axes,axes3d,graph,surface_plane,time,Force2d,Force3d,Total_2DGroup,Total_3DGroup]
                
    def construct(self):
        # Setup Variables
        length = 20
        resolution_fa = 100 # Set for more point resolution on the wave plot
        runtime = 13
        sig = ['wave','sine']

        axes2D = []
        axes3D = []
    
        graph2D = []
        surface3D = []
        time_var = []
        
        Group_2D = []
        Group_3D = []
        Group_Decorators = []
        Force2D = []
        Force3D = []

        # phi, theta, focal_distance, gamma, distance_to_origin = self.camera.get_value_trackers()

        
        # Setup Axes
        for ii in range(len(sig)):  
            [axes,axes3d,graph,surface_plane,time,force2d,force3d,Total_2DGroup,Total_3DGroup] = self.Create_Plots(resolution_fa,length,sig[ii])
            axes2D.append(axes)
            axes3D.append(axes3d)
            graph2D.append(graph)
            surface3D.append(surface_plane)
            time_var.append(time)
            
            # Setup 3D Variables     
            Force2D.append(force2d)
            Force3D.append(force3d)

            # Setup 2D Variables
            phi, theta, focal_distance, gamma, distance_to_origin = self.camera.get_value_trackers()
            # Total_2DGroup = Group(axes2D[ii],Force2D[ii],graph2D[ii]).move_to([0,-2,0]).scale(0.4)
            # Total_3DGroup = Group(BC,axes3D[ii],Force3D[ii],surface3D[ii]).move_to([0,2,0]).scale(0.4)
            
            
            Group_2D.append(Total_2DGroup)
            Group_3D.append(Total_3DGroup)

        # Animate Pre-Math
        self.Math_Text()
  
        self.play(FadeIn(Group_2D[0]),FadeIn(Group_3D[0]))
        for ii in range(len(sig)):
            print(ii)
            tt = 0
            on_screen_var = Variable(tt, Text("Time"), num_decimal_places=3).scale(0.5).shift(LEFT*1.5).shift(DOWN/2)
            time_var_text = Text("s").next_to(on_screen_var,RIGHT)
            on_screen_var.label.set_color(WHITE)
            on_screen_var.value.set_color(GREEN)
            var_tracker = on_screen_var.tracker
            
            
            self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES,distance=1)
            # self.set_camera_orientation(phi=90 * DEGREES, theta=-90 * DEGREES,distance=1)
            if ii == 0:
                self.play(FadeIn(on_screen_var))
            else:
                # self.play(FadeIn(graph2D[ii]),FadeIn(surface3D[ii]),FadeIn(Force2D[ii]),FadeIn(Force3D[ii]),FadeIn(on_screen_var))
                self.play(FadeIn(Group_2D[ii]),FadeIn(Group_3D[ii]),FadeIn(on_screen_var))

            
            self.wait(2)
            self.play(time_var[ii].animate.set_value(1000),var_tracker.animate.set_value(runtime),run_time=runtime,rate_func=linear)
            #self.play(FadeOut(graph2D[ii]),FadeOut(surface3D[ii]),FadeOut(Force2D[ii]),FadeOut(Force3D[ii]),FadeOut(on_screen_var))
            self.play(FadeOut(Group_2D[ii]),FadeOut(Group_3D[ii]),FadeOut(on_screen_var))

class FE_DifferenceTot_KeyFrame(ThreeDScene):
    def Finite_Difference_Total_Code(self):
        import numpy as np
        import matplotlib.pyplot as plt
        
        Nx = 101
        dx = 1
        
        f  = 10
        Nt = 1001
        dt = 1/(f*100)
        c  = 343
        
        
        x  = np.linspace(0,(Nx-1)*dx,Nx)
        t  = np.linspace(0,(Nt-1)*dt,Nt)
        
        C  = c*(dt/dx)
        
        U  = np.zeros((Nt,Nx))
        s1 = int(Nt/f)
        
        # U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
        # U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
        U[:,0] = np.sin(2*np.pi*f*t)

        
        
        for jj in range(2,Nt):
            for nn in range(1,Nx-1):
                U1 = 2*U[jj-1,nn] - U[jj-2,nn]
                U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
                U[jj,nn] = U1 + C*C*U2
                U[:,-1] = U[:,-2]
        U = np.flip(U, 1)
        
        return U
    
    def Create_Func2D(self):
        U = self.Finite_Difference_Total_Code()

        return U
    
    def Create_Func3D(self):
        U = self.Finite_Difference_Total_Code()
        Nx = len(U[1,:])
        Nt = len(U)

        U_arr = []
        for ii in range(Nt):
            U_arr.append(funMatLabMultip(U[ii,:],np.ones(Nx))) 

        
        return U_arr
    
    def Fun_2D(self,U, t, u):
        x = int(u)
        t = int(t)

        return U[t,x]
    
    def Func_3d(self,U, t, u, v):
        x = int(u)
        y = int(v)
        t = int(t)

        return U[t][x,y]
    
    def Func_3d_Plane(self,u,v):
        x = u
        y = v
        
        return 0
    
    def Math_Text(self):
        pt1 = MathTex(r"\frac{d^2 \xi}{dx^2}")
        
        pt2 = MathTex(r"\frac{1}{c^2}\frac{d^2 \xi}{dt^2}")
        
        pt3 = MathTex(" = 0")
        variables_group = VGroup(pt1,pt2,pt3).arrange_submobjects().shift(UP)

        eqn1 = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn2 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn3 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", " = 0")
        eqn4 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "=", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", "")
        eqn5 = MathTex(r"\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}", "=", r" \left(\frac{\Delta x^2}{\Delta t^2c^2}\right)", r"\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}")
        eqn6 = MathTex(r"\xi_{j}^{n+1}", "=", r"\left(\frac{c \Delta t}{\Delta x}\right)^2", r"\left(\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}\right)", "+" , r" 2 \xi_{j}^{n} - \xi_{j}^{n-1}")

        text1 = MathTex("CFL")
        text1.set_color(GREEN)
        text1.next_to(eqn6[2],DOWN)
        eqn6[2].set_color(GREEN)
        group1 = VGroup(eqn6,text1)
        group1.scale(0.5).shift(UP)

        self.play(Write(eqn1))


        eqn1[0].set_color(YELLOW)
        eqn2[0].set_color(YELLOW)

        self.wait(2)

        self.play(Transform(eqn1,eqn1))

        self.play(TransformMatchingTex(eqn1,eqn2))
        self.wait(2)
        eqn2[0].set_color(WHITE)
        eqn2[2].set_color(YELLOW)
        eqn3[2].set_color(YELLOW)

        self.play(TransformMatchingTex(eqn2,eqn2))
        self.play(TransformMatchingTex(eqn2,eqn3))
        self.wait(2)
        self.play(TransformMatchingTex(eqn3,eqn4))
        self.wait(2)
        self.play(TransformMatchingTex(eqn4,eqn5))
        self.wait(2)
        self.play(TransformMatchingTex(eqn5,group1))

        #self.play(group1.animate.shift(3*UP))
        
    def construct(self):
        # 2D Setup
        U = self.Create_Func2D()
        time = ValueTracker(250)
        axes = Axes(x_range = [0,100,10],
                     x_length = 20,
                     y_range = [-2,2,1],
                     y_length= 5)

        graph = always_redraw(lambda:
                              axes.plot(lambda x: self.Fun_2D(U,time.get_value(),x),
                                        x_range=[0,100],
                                        color = GREEN,
                                        use_smoothing=False))
        Force2dArrow = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(axes,RIGHT)
        F2d = MathTex(r"F = sin(\omega t)").next_to(Force2dArrow,RIGHT)
        Force2d = Group(Force2dArrow,F2d)
        Total_2DGroup = Group(axes,graph,Force2d).move_to([0,-2,0]).scale(0.4)


        
        # 3D Set up
        length = 20
        U_arr = self.Create_Func3D()
        axes3d = ThreeDAxes(x_range=(0, 100, 1), y_range=(0, 100, 1), z_range=(-2, 2, 1), x_length=length, y_length = 5, z_length = 5)
        resolution_fa = 25 # Set for more point resolution on the wave plot
        phi, theta, focal_distance, gamma, distance_to_origin = self.camera.get_value_trackers()

        surface_plane = always_redraw(lambda : Surface(
            lambda u, v: axes3d.c2p(u, v, self.Func_3d(U_arr,time.get_value(),u,v)),
            resolution=(resolution_fa, resolution_fa),
            u_range=[0, 100],
            v_range=[0, 100],
            should_make_jagged = True, 
        ))

        BC = Surface(
            lambda u, v: axes3d.c2p(self.Func_3d_Plane(u,v), u, v),
            resolution=(1, 1),
            u_range=[0, 100],
            v_range=[-2.2, 2.2],
            should_make_jagged = True, 
            checkerboard_colors  = False,
            fill_opacity  = 1,
            fill_color = GRAY
        )

        Fixed3d = MathTex(r'Fixed: \frac{d\xi}{dx} = 0').scale(1.2).shift(LEFT*(int(length/2))).shift(DOWN*5)
        Free3d =  MathTex(r'Free: \xi = 0').scale(1.2).shift(RIGHT*(int(length/2))).shift(DOWN*5)
        Force3dArrow = Arrow(start=RIGHT, end=LEFT, color=RED).next_to(axes3d,RIGHT)
        F3d = MathTex(r"F = sin(\omega t)",color=RED).next_to(Force3dArrow,DOWN)
        Force3d = Group(Force3dArrow,F3d)
        
        
        
        
        prismSmall = Prism(dimensions=[5,length,3],fill_color=LIGHT_GRAY,stroke_color=BLACK).rotate(PI / 2)
        prismSmall.set_style(fill_opacity=0.3)

        surface_plane.set_style(fill_opacity=1)
        surface_plane = always_redraw(lambda : surface_plane.set_fill_by_value(axes=axes3d, colorscale=[ (BLUE,-2), (BLUE,-1), (BLUE,-0.5), (GREEN, -0.1), (GREEN, 0), (GREEN, 0.1), (RED, 0.5), (RED, 1), (RED, 2)], axis=2))
   
        BigSpace_PlotGroup = Group(prismSmall,Fixed3d,Free3d,BC,surface_plane,axes3d,Force3d).scale(0.4).move_to([0,2,0])
        
        runtime = 5
        tt = 0
        on_screen_var = Variable(tt, Text("Time"), num_decimal_places=3).scale(0.5).shift(LEFT).shift(DOWN/2)
        time_var_text = Text("s").next_to(on_screen_var,RIGHT)
        on_screen_var.label.set_color(WHITE)
        on_screen_var.value.set_color(GREEN)
        var_tracker = on_screen_var.tracker
        
        self.set_camera_orientation(phi=0 * DEGREES, theta=-90 * DEGREES,distance=1)
        self.add(BigSpace_PlotGroup,Total_2DGroup,on_screen_var)













#%% OLD 
class FE_DifferenceOLD(Scene):
    def Create_Func(self):
        import numpy as np
        import matplotlib.pyplot as plt
        
        Nx = 101
        dx = 1
        
        f  = 10
        Nt = 1001
        dt = 1/(f*100)
        c  = 343
        
        
        x  = np.linspace(0,(Nx-1)*dx,Nx)
        t  = np.linspace(0,(Nt-1)*dt,Nt)
        
        C  = c*(dt/dx)
        
        U  = np.zeros((Nt,Nx))
        s1 = int(Nt/f)
        
        U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
        U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
        
        
        
        for jj in range(2,Nt):
            for nn in range(1,Nx-1):
                U1 = 2*U[jj-1,nn] - U[jj-2,nn]
                U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
                U[jj,nn] = U1 + C*C*U2
                U[:,-1] = U[:,-2]
        U = np.flip(U, 1)

        return 4*U
    
    def U_func(self,idx,U):
        idx = int(idx)
        f = FunctionGraph(
            lambda x: U[idx,int(x)],
            x_range=[0,100]
            ).move_to([-20,0,1])
        print(idx)
        return f
            
    def Math_Text(self):
        pt1 = MathTex(
            r"\frac{d^2 \xi}{dx^2}"
        )
        
        pt2 = MathTex(
           r"\frac{1}{c^2}\frac{d^2 \xi}{dt^2}"
        )
        
        pt3 = MathTex(" = 0")
        variables_group = VGroup(pt1,pt2,pt3).arrange_submobjects().shift(UP)

        eqn1 = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn2 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn3 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", " = 0")
        eqn4 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "=", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", "")
        eqn5 = MathTex(r"\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}", "=", r" \left(\frac{\Delta x^2}{\Delta t^2c^2}\right)", r"\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}")
        text1 = MathTex("CFL")
        text1.set_color(YELLOW)
        text1.next_to(eqn5[2],DOWN)
        eqn5[2].set_color(YELLOW)
        group1 = VGroup(eqn5,text1)

        self.play(Write(eqn1))


        eqn1[0].set_color(YELLOW)
        eqn2[0].set_color(YELLOW)

        self.wait(2)

        self.play(Transform(eqn1,eqn1))

        self.play(TransformMatchingTex(eqn1,eqn2))
        self.wait(2)
        eqn2[0].set_color(WHITE)
        eqn2[2].set_color(YELLOW)
        eqn3[2].set_color(YELLOW)

        self.play(TransformMatchingTex(eqn2,eqn2))
        self.play(TransformMatchingTex(eqn2,eqn3))
        self.wait(2)
        self.play(TransformMatchingTex(eqn3,eqn4))
        self.wait(2)
        self.play(TransformMatchingTex(eqn4,group1))
        self.wait(2)
        self.play(group1.animate.shift(3*UP))

        
    def construct(self):
        U = self.Create_Func()
        U_function = self.U_func(0,U).scale(0.12).center()
        time = ValueTracker(0)  
        axes = Axes(x_range = [0,100,5], y_range = [-1,1])
        
        y_label = axes.get_y_axis_label(MathTex(r"\xi"), edge=LEFT, direction=LEFT, buff=0.4)
        x_label = axes.get_x_axis_label("x (position)")
        grid_labels = VGroup(x_label, y_label)
        def update_U(U_function):
            U_function.become(
                self.U_func(time.get_value(),U).scale(0.12).center()
                )
            return U_function
 
        U_function.add_updater(update_U)
        graph = U_function

        pt1 = MathTex(
            r"\frac{d^2 \xi}{dx^2}"
        )
        
        pt2 = MathTex(
           r"\frac{1}{c^2}\frac{d^2 \xi}{dt^2}"
        )
        
        pt3 = MathTex(" = 0")
        variables_group = VGroup(pt1,pt2,pt3).arrange_submobjects().shift(UP)

        eqn1 = MathTex(r"\frac{\partial^2 \xi}{\partial x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn2 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\partial^2 \xi}{\partial t^2}", " = 0")
        eqn3 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "-", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", " = 0")
        eqn4 = MathTex(r"\frac{\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}}{\Delta x^2}", "=", r"\frac{1}{c^2}\frac{\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}}{\Delta t^2}", "")
        eqn5 = MathTex(r"\xi_{j+1}^{n} - 2 \xi_{j}^{n} + \xi_{j-1}^{n}", "=", r" \left(\frac{\Delta x^2}{\Delta t^2c^2}\right)", r"\xi_{j}^{n+1} - 2 \xi_{j}^{n} + \xi_{j}^{n-1}")
        text1 = MathTex("CFL")
        text1.set_color(YELLOW)
        text1.next_to(eqn5[2],DOWN)
        eqn5[2].set_color(YELLOW)
        group1 = VGroup(eqn5,text1)

        self.play(Write(eqn1))


        eqn1[0].set_color(YELLOW)
        eqn2[0].set_color(YELLOW)

        self.wait(2)

        self.play(Transform(eqn1,eqn1))

        self.play(TransformMatchingTex(eqn1,eqn2))
        self.wait(2)
        eqn2[0].set_color(WHITE)
        eqn2[2].set_color(YELLOW)
        eqn3[2].set_color(YELLOW)

        self.play(TransformMatchingTex(eqn2,eqn2))
        self.play(TransformMatchingTex(eqn2,eqn3))
        self.wait(2)
        self.play(TransformMatchingTex(eqn3,eqn4))
        self.wait(2)
        self.play(TransformMatchingTex(eqn4,group1))
        self.wait(2)
        self.play(group1.animate.shift(3*UP))
        
        self.play(Create(axes),Create(graph),Create(grid_labels))
        self.wait(2)
        self.play(time.animate.set_value(1000),run_time=10,rate_func=linear)


import numpy as np
import matplotlib.pyplot as plt


def Create_Func():
    import numpy as np
    import matplotlib.pyplot as plt
    
    Nx = 101
    dx = 1
    
    f  = 10
    Nt = 1001
    dt = 1/(f*100)
    c  = 343
    
    
    x  = np.linspace(0,(Nx-1)*dx,Nx)
    t  = np.linspace(0,(Nt-1)*dt,Nt)
    
    C  = c*(dt/dx)
    
    U  = np.zeros((Nt,Nx))
    s1 = int(Nt/f)
    
    U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
    U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])
    
    
    
    for jj in range(2,Nt):
        for nn in range(1,Nx-1):
            U1 = 2*U[jj-1,nn] - U[jj-2,nn]
            U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
            U[jj,nn] = U1 + C*C*U2

    return U

def U_func(idx,U):
    idx = int(idx)
    return U[idx,:]



U = Create_Func()
U_function = lambda t: U_func(t,U)

        
        

#%%

Nx = 101
dx = 1

f  = 10
Nt = 1001
dt = 1/(f*100)
c  = 343


x  = np.linspace(0,(Nx-1)*dx,Nx)
t  = np.linspace(0,(Nt-1)*dt,Nt)

C  = c*(dt/dx)

U  = np.zeros((Nt,Nx))
s1 = int(Nt/f)

# U[len(t)-s1:len(t),0] = np.sin(2*np.pi*f*t[0:s1])
# U[len(t)-s1:len(t),1] = np.sin(2*np.pi*f*t[0:s1])
U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])

for jj in range(2,Nt):
    for nn in range(1,Nx-1):
        U1 = 2*U[jj-1,nn] - U[jj-2,nn]
        U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
        U[jj,nn] = U1 + C*C*U2

U = np.flip(U, 1)

times = [100,200,450]
for jj in range(len(times)):
    plt.plot(x,U[times[jj],:])



#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

def funMatLabMultip(f,t):
    """create an n x m array from 2 vectors of size n and m.
    Resulting rows are the multiplication of each element of the first vector for all the elements of the second vector
    f=np.array([2,4])
    t=np.array([1,2,3,4,5])
    [[ 2  4  6  8 10]
     [ 4  8 12 16 20]]
    """
    if t.size==t.shape[0]:
        k=f[0]*t
        for i in f[1:]:
            j=i*t
            k=np.vstack((k,j))
    else:
        raise Exception('arrays should 1D arrays')
    return k

Nx = 101
dx = 1

f  = 10
Nt = 1001
dt = 1/(f*100)
c  = 343


x  = np.linspace(0,(Nx-1)*dx,Nx)
t  = np.linspace(0,(Nt-1)*dt,Nt)

C  = c*(dt/dx)

U  = np.zeros((Nt,Nx))
s1 = int(Nt/f)

U[0:s1,0] = np.sin(2*np.pi*f*t[0:s1])
U[0:s1,1] = np.sin(2*np.pi*f*t[0:s1])



for jj in range(2,Nt):
    for nn in range(1,Nx-1):
        U1 = 2*U[jj-1,nn] - U[jj-2,nn]
        U2 = U[jj-1,nn-1] - 2*U[jj-1,nn] + U[jj-1,nn+1]
        U[jj,nn] = U1 + C*C*U2
        U[:,-1] = U[:,-2]
U = np.flip(U, 1)

U_arr = []
for ii in range(Nt):
    U_arr.append(4*funMatLabMultip(U[ii,:],np.ones(Nx))) 
# Plot the surface.
X, Y = np.meshgrid(x, x)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(X, Y, U_arr[25], cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)


imageArray = np.uint8(100*U_arr[50])        



#%%
def func(u, v):
    return np.array([np.cos(u) * np.cos(v), np.cos(u) * np.sin(v), u])


aasdas = func(1,10)